"""Run, resume and summarize the paired RGLL/tree simulation in Python."""

import argparse
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import csv
import itertools
import json
import os
import platform
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import torch

from lmm import Evaluator
from search import Search, predict_group, state_id


HERE = Path(__file__).resolve().parent
# Figure numbering in RGLL/export_paper_figures.R.
OLD_CASE = (11, 12, 13, 14, 15, 3, 16, 1, 2, 4, 5, 6, 7, 8, 9, 10)
TREE_METHODS = ("GrowThenPrune", "TreeAnnealing", "BeamSearch")
RGLL_METHODS = ("ORACLE", "CIRG", "KM", "BLM", "CPF", "OBS", "LM")
SETTING = ["case", "rgll_case", "model", "R_true", "var_a", "var_b"]
PAIR = SETTING + ["replication", "seed"]


@dataclass(frozen=True)
class Job:
    case: int
    r: int
    var_a: float
    var_b: float
    rep: int

    @property
    def key(self):
        return f"case{self.case}_R{self.r}_A{self.var_a:g}_B{self.var_b:g}_rep{self.rep}"

    def metadata(self):
        return dict(case=self.case, rgll_case=OLD_CASE[self.case - 1], model="RS" if self.var_b else "RI",
                    R_true=self.r, var_a=self.var_a, var_b=self.var_b, replication=self.rep,
                    seed=12345 + 1000003 * self.rep + 7919 * 2500)


class RGLL:
    def __init__(self, executable, root):
        self.executable = executable
        self.root = root

    def __enter__(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="imspe-tree-")
        self.file = Path(self.tmp.name) / "transfer"
        self.proc = subprocess.Popen([self.executable, "--vanilla", str(HERE / "rgll_bridge.R"), str(self.root)],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     text=True, encoding="utf-8", errors="replace", bufsize=1)
        try:
            self._await("READY")
        except BaseException:
            self.__exit__()
            raise
        return self

    def __exit__(self, *error):
        try:
            if self.proc.poll() is None:
                self.proc.stdin.write("quit\n")
                self.proc.stdin.flush()
                self.proc.wait()
        finally:
            self.proc.stdin.close()
            self.proc.stdout.close()
            self.tmp.cleanup()

    def _await(self, token):
        for line in self.proc.stdout:
            line = line.strip()
            if line == token:
                return
            if line.startswith("ERROR "):
                raise RuntimeError(line[6:])
            if line:
                print(f"RGLL: {line}", flush=True)
        raise RuntimeError(f"RGLL bridge exited with status {self.proc.poll()}")

    def call(self, mode, job, *args):
        rgll_case = OLD_CASE[job.case - 1]
        fields = [mode, rgll_case, job.r, job.var_a, job.var_b, job.rep, self.file.as_posix(), *args]
        self.proc.stdin.write("\t".join(map(str, fields)) + "\n")
        self.proc.stdin.flush()
        self._await("OK")
        return self.file

    def data(self, job):
        vector = np.fromfile(self.call("data", job), dtype="<f8")
        offsets = np.cumsum([7500 * 50, 2500 * 50, 2500 * 50, 7500, 2500, 7500, 2500])
        if len(vector) != offsets[-1] or not np.isfinite(vector).all():
            raise ValueError("Invalid RGLL data transfer")
        x, ref, test, y, y_test, group, group_test = np.split(vector, offsets[:-1])
        return (x.reshape(-1, 50), ref.reshape(-1, 50), test.reshape(-1, 50), y, y_test,
                group.astype(np.int32) - 1, group_test.astype(np.int32) - 1)

    def baseline(self, job, methods):
        file = self.call("live", job, ",".join(methods))
        with file.open(newline="", encoding="utf-8") as handle:
            rows = []
            for d in csv.DictReader(handle):
                number = lambda k: None if d[k] in {"NA", ""} else float(d[k])
                rows.append(dict(method=d["method"], mspe=number("MSPE"), inferred_K=number("selected_K"),
                                 elapsed_s=number("runtime"), ari_test=number("group_ARI"),
                                 imspe=number("IMSPE_selected"), em_iterations=number("iterations"),
                                 converged=None if d["converged"] == "NA" else d["converged"] == "TRUE",
                                 backend="R/C++ CPU", timing_source="local",
                                 status="failed" if d["error"] else "ok", error=d["error"]))
        return rows


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, allow_nan=False)
    temporary.replace(path)


def adjusted_rand(a, b):
    _, a = np.unique(a, return_inverse=True)
    _, b = np.unique(b, return_inverse=True)
    k = b.max() + 1
    table = np.bincount(a * k + b, minlength=(a.max() + 1) * k).reshape(-1, k)
    pairs = lambda x: (x * (x - 1) / 2).sum()
    shared = pairs(table)
    row, col = pairs(table.sum(1)), pairs(table.sum(0))
    chance = row * col / pairs(np.array([len(a)], dtype=np.int64))
    upper = (row + col) * 0.5
    return float((shared - chance) / (upper - chance)) if upper != chance else 1.0


def make_evaluator(x, y, ref, job, config, device, solver="torch"):
    return Evaluator(x, y, ref, random_slope=bool(job.var_b), device=device,
                     tolerance=config["em_tolerance"], iterations=config["em_iterations"], solver=solver)


def run_trees(d, job, config, device, solver, directory):
    x, ref, test, y, ytest, group, group_test = d
    init_seed = job.metadata()["seed"] + 1777
    metrics = []
    for index, method in enumerate(TREE_METHODS):
        t0 = perf_counter()
        evaluator = make_evaluator(x, y, ref, job, config, device, solver)
        search = Search(x, ref, evaluator, config, init_seed + 1001 * (index + 1))
        best, termination = search.run(method, init_seed)
        backend = "PyTorch " + device if solver == "torch" else f"PyTorch {device} + C++ CPU fit"
        row = dict(method=method, backend=backend, timing_source="local",
                   termination=termination, evaluations=search.evaluations, status="ok" if best else "failed")
        if best:
            test_labels = predict_group(best.tree, test)
            prediction = best.fit.predict(test, best.remap[test_labels])
        if device == "cuda":
            torch.cuda.synchronize()
        row["elapsed_s"] = perf_counter() - t0
        if best:
            train_labels = predict_group(best.tree, x)
            row.update(mspe=float(np.mean(np.square(prediction - ytest))), imspe=best.fit.imspe,
                       inferred_K=best.k, depth=best.depth, min_leaf_size=min(best.counts),
                       max_leaf_size=max(best.counts), ari_train=adjusted_rand(group, train_labels),
                       ari_test=adjusted_rand(group_test, test_labels),
                       converged=best.fit.converged, em_iterations=best.fit.iterations)
            atomic_json(directory / f"{method}_best.json",
                        dict(candidate_id=state_id(best.tree), tree=asdict(best.tree),
                             leaf_sizes=best.counts, fit_leaf_order=best.remap.tolist(),
                             beta=best.fit.beta.tolist(), effects=best.fit.effects.tolist(),
                             variance=list(map(float, best.fit.variance)), metrics=row))
            by_k = {}
            for h in search.history:
                if h.get("imspe") is not None:
                    k, score = h["k"], h["imspe"]
                    by_k[k] = min(by_k.get(k, np.inf), score)
            pd.DataFrame(sorted(by_k.items()), columns=["K", "best_IMSPE"]).to_csv(
                directory / f"{method}_imspe_vs_k.csv", index=False)
        pd.DataFrame(search.history).assign(method=method, **job.metadata(),
                                            search_seed=init_seed + 1001 * (index + 1)).to_csv(
            directory / f"{method}_history.csv", index=False)
        metrics.append(row)
    return metrics


def summarize(directory):
    records = []
    for job_dir in sorted(directory.glob("runs/*")):
        for filename in ("tree_metrics.json", "baseline_metrics.json"):
            path = job_dir / filename
            if path.exists():
                records.extend(json.loads(path.read_text(encoding="utf-8")))
    if not records:
        return
    d = pd.DataFrame(records)
    for column in ("mspe", "elapsed_s", "inferred_K", "imspe", "ari_train", "ari_test", "evaluations"):
        if column not in d:
            d[column] = np.nan
    d["delta_K_vs_true"] = d.inferred_K - d.R_true
    d.to_csv(directory / "replications.csv", index=False, encoding="utf-8-sig")
    id_cols = SETTING + ["method", "backend", "timing_source"]
    summary = d.groupby(id_cols, dropna=False).agg(
        n_runs=("status", "size"), n_mspe=("mspe", "count"), mspe_mean=("mspe", "mean"),
        mspe_sd=("mspe", "std"), seconds_mean=("elapsed_s", "mean"),
        seconds_median=("elapsed_s", "median"), seconds_sd=("elapsed_s", "std"),
        k_mean=("inferred_K", "mean"), k_min=("inferred_K", "min"), k_max=("inferred_K", "max"),
        delta_k_mean=("delta_K_vs_true", "mean"), ari_train_mean=("ari_train", "mean"),
        ari_test_mean=("ari_test", "mean"), imspe_mean=("imspe", "mean"),
        evaluations_mean=("evaluations", "mean")).reset_index()
    interval = 1.96 * summary.mspe_sd / np.sqrt(summary.n_mspe.clip(lower=1))
    summary["mspe_mc95_low"], summary["mspe_mc95_high"] = summary.mspe_mean - interval, summary.mspe_mean + interval
    summary.to_csv(directory / "summary.csv", index=False, encoding="utf-8-sig")
    tree = d[d.method.isin(TREE_METHODS)]
    baseline = d[d.method.isin(RGLL_METHODS) & (d.timing_source == "local")]
    if tree.empty or baseline.empty:
        return
    pairs = tree.merge(baseline, on=PAIR, suffixes=("_tree", "_baseline"))
    columns = PAIR + ["method_tree", "method_baseline", "timing_source_baseline", "mspe_tree", "mspe_baseline",
                      "delta_mspe", "inferred_K_tree", "inferred_K_baseline", "ari_test_tree", "ari_test_baseline",
                      "elapsed_s_tree", "elapsed_s_baseline", "local_seconds_ratio", "delta_K"]
    pairs["delta_mspe"] = pairs.mspe_tree - pairs.mspe_baseline
    pairs["delta_K"] = pairs.inferred_K_tree - pairs.inferred_K_baseline
    paired_local = (pairs.timing_source_baseline == "local") & (pairs.elapsed_s_baseline > 0)
    pairs["local_seconds_ratio"] = pairs.elapsed_s_tree / pairs.elapsed_s_baseline.where(paired_local)
    pairs[columns].to_csv(directory / "paired_comparisons.csv", index=False, encoding="utf-8-sig")
    by_case = pairs.groupby(SETTING + ["method_tree", "method_baseline", "timing_source_baseline"], dropna=False).agg(
        n_pairs=("delta_mspe", "count"), delta_mspe_mean=("delta_mspe", "mean"),
        delta_mspe_sd=("delta_mspe", "std"), mspe_tree_mean=("mspe_tree", "mean"),
        mspe_baseline_mean=("mspe_baseline", "mean"),
        seconds_tree_mean=("elapsed_s_tree", "mean"), seconds_baseline_mean=("elapsed_s_baseline", "mean"),
        local_seconds_ratio_median=("local_seconds_ratio", "median"), delta_K_mean=("delta_K", "mean"),
        ari_test_tree_mean=("ari_test_tree", "mean"), ari_test_baseline_mean=("ari_test_baseline", "mean"),
    ).reset_index()
    difference = 1.96 * by_case.delta_mspe_sd / np.sqrt(by_case.n_pairs.clip(lower=1))
    by_case["delta_mspe_mc95_low"] = by_case.delta_mspe_mean - difference
    by_case["delta_mspe_mc95_high"] = by_case.delta_mspe_mean + difference
    by_case.to_csv(directory / "paired_summary.csv", index=False, encoding="utf-8-sig")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rgll", type=Path, default=HERE.parent.parent / "RGLL")
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", shutil.which("Rscript")))
    parser.add_argument("--config", type=Path, default=HERE / "config.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--solver", choices=("native", "torch"), default="native")
    parser.add_argument("--cpu-threads", type=int, default=2)
    parser.add_argument("--phase", choices=("all", "baseline", "trees"), default="all")
    parser.add_argument("--baseline-workers", type=int, default=1)
    parser.add_argument("--live-methods", nargs="+", choices=RGLL_METHODS, default=RGLL_METHODS)
    parser.add_argument("--case", type=int, nargs="+", choices=range(1, len(OLD_CASE) + 1))
    parser.add_argument("--r", type=int, nargs="+", choices=(10, 20, 50))
    parser.add_argument("--model", choices=("RI", "RS"))
    parser.add_argument("--var-a", type=float, nargs="+", choices=(0.5, 2.25))
    parser.add_argument("--limit", type=int, help="Maximum selected jobs for a test batch")
    parser.add_argument("--replications", type=int, choices=range(1, 6))
    parser.add_argument("--summarize-only", action="store_true")
    return parser.parse_args()


def read_config(path):
    config = json.loads(path.read_text(encoding="utf-8"))
    base_name = config.pop("extends", None)
    if base_name:
        base = json.loads((path.parent / base_name).read_text(encoding="utf-8"))
        config = {**base, **config}
    if not config["cases"] or len(set(config["cases"])) != len(config["cases"]) or any(
            type(case) is not int or not 1 <= case <= len(OLD_CASE) for case in config["cases"]):
        raise ValueError("cases must be distinct RGLL figure numbers in 1–16")
    search = config["search"]
    counts = ("min_leaf", "initial_candidates", "candidates_per_parent", "max_evaluations",
              "grow_evaluations", "fit_batch", "beam_width", "anneal_proposals", "em_iterations")
    for name in counts:
        if type(search[name]) is not int or search[name] <= 0:
            raise ValueError(f"search.{name} must be a positive integer")
    if not 2 <= search["initial_candidates"] < search["grow_evaluations"] < search["max_evaluations"]:
        raise ValueError("Evaluation budgets must reserve room for initialization, growth and pruning")
    if search["min_leaf"] > 3750:
        raise ValueError("min_leaf leaves no room for a two-leaf initialization")
    for name in ("anneal_temperature_fraction", "em_tolerance"):
        if not np.isfinite(search[name]) or search[name] <= 0:
            raise ValueError(f"search.{name} must be positive and finite")
    if not 0 < search["anneal_cooling"] <= 1:
        raise ValueError("search.anneal_cooling must be between 0 (exclusive) and 1")
    probabilities = np.asarray(search["anneal_move_probabilities"], dtype=float)
    if probabilities.shape != (3,) or not np.all(probabilities >= 0) or not np.isclose(probabilities.sum(), 1):
        raise ValueError("Annealing needs three nonnegative move probabilities summing to one")
    if type(config["repetitions"]) is not int or not 1 <= config["repetitions"] <= 5:
        raise ValueError("Use between one and five replications")
    return config


def do_baselines(jobs, args, root):
    pending = deque(job for job in jobs if not (root / "runs" / job.key / "baseline_metrics.json").exists())
    workers = min(args.baseline_workers, len(pending))

    def worker():
        with RGLL(args.rscript, args.rgll.resolve()) as rgll:
            while True:
                try:
                    job = pending.popleft()
                except IndexError:
                    return
                directory = root / "runs" / job.key
                file = directory / "baseline_metrics.json"
                if file.exists():
                    continue
                directory.mkdir(parents=True, exist_ok=True)
                t0 = perf_counter()
                baseline = rgll.baseline(job, args.live_methods)
                for row in baseline:
                    row["baseline_workers"] = args.baseline_workers
                atomic_json(file, [{**job.metadata(), **row} for row in baseline])
                print(f"RGLL {job.key} in {perf_counter() - t0:.1f}s", flush=True)

    if workers:
        with ThreadPoolExecutor(workers) as pool:
            futures = [pool.submit(worker) for _ in range(workers)]
            for future in futures:
                future.result()


def main():
    args = arguments()
    root = args.output.resolve()
    if args.summarize_only:
        summarize(root)
        return
    if not args.rscript:
        raise SystemExit("Provide --rscript or set RSCRIPT to Rscript.exe")
    if min(args.cpu_threads, args.baseline_workers) <= 0:
        raise ValueError("Thread and worker counts must be positive")
    config = read_config(args.config)
    s = config["search"]
    device = args.device if args.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu")
    if args.phase != "baseline" and device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA backend requested but CUDA is unavailable")
    torch.set_num_threads(args.cpu_threads)
    # Load CUDA libraries before recording method runtime.
    if args.phase != "baseline" and device == "cuda":
        t = torch.eye(51, device=device, dtype=torch.float64).expand(8, 51, 51).contiguous()
        torch.linalg.solve_ex(t, t)
        torch.cuda.synchronize()
    jobs = itertools.product(args.case or config["cases"], args.r or config["generating_groups"],
                             args.var_a or config["intercept_variances"],
                             [0 if args.model == "RI" else 0.1] if args.model else config["slope_variances"],
                             range(1, 1 + (args.replications or config["repetitions"])))
    jobs = [Job(*j) for j in jobs][:args.limit]
    spec = dict(config=config, device=device, cpu_threads=args.cpu_threads, torch=str(torch.__version__),
                gpu=torch.cuda.get_device_name() if device == "cuda" else None, platform=platform.platform(),
                baseline="live", live_methods=list(args.live_methods), archive_and_live=[],
                baseline_workers=args.baseline_workers,
                rgll=str(args.rgll.resolve()))
    root.mkdir(parents=True, exist_ok=True)
    manifest = root / "run_config.json"
    if manifest.exists() and spec != json.loads(manifest.read_text(encoding="utf-8")):
        raise SystemExit("Run settings differ from the existing output. Use a different output directory.")
    atomic_json(manifest, spec)
    if args.phase != "baseline":
        solver_manifest = root / "tree_backend.json"
        if solver_manifest.exists() and json.loads(solver_manifest.read_text())["solver"] != args.solver:
            raise SystemExit("The existing output uses a different solver")
        if args.solver == "native":
            from native import kernel
            kernel()
        atomic_json(solver_manifest, {"solver": args.solver})
    try:
        if args.phase == "baseline":
            do_baselines(jobs, args, root)
            return
        with RGLL(args.rscript, args.rgll.resolve()) as rgll:
            for number, job in enumerate(jobs, 1):
                directory = root / "runs" / job.key
                if (directory / "tree_metrics.json").exists():
                    continue
                directory.mkdir(parents=True, exist_ok=True)
                t0 = perf_counter()
                data = rgll.data(job)
                trees = run_trees(data, job, s, device, args.solver, directory)
                atomic_json(directory / "tree_metrics.json", [{**job.metadata(), **row} for row in trees])
                print(f"[{number}/{len(jobs)}] {job.key} in {perf_counter() - t0:.1f}s  " +
                      "  ".join(f"{r['method']}={r.get('mspe', '-')}" for r in trees), flush=True)
        if args.phase == "all":
            do_baselines(jobs, args, root)
    finally:
        summarize(root)


if __name__ == "__main__":
    main()
