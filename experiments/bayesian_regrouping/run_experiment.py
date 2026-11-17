"""Run one five-repeat setting using the original RGLL paper data."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
from time import perf_counter

import numpy as np
import pandas as pd
import torch
from scipy.stats import t
from scipy.special import gammaln

from evaluation import (Assignment, conditional_imspe, discrete_summary, load_draws,
                        representative, scalar_diagnostic, summarize_numbers)
from latent_lmm import fit_chain, peak_memory_bytes, scale_x, write_npz

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[1]
DEFAULT_RGLL = WORKSPACE if (WORKSPACE / "ROG_method_comparison.R").is_file() else WORKSPACE / "RGLL"
sys.path.insert(0, str(HERE.parent / "imspe_tree"))
from run import RGLL, Job, adjusted_rand, atomic_json, make_evaluator

METHODS = ("mfm", "dpm", "ppmx")


def load_config(path, variant):
    config = json.loads(path.read_text(encoding="utf-8"))
    if variant:
        if variant not in config["sensitivity"]:
            raise ValueError(f"Unknown sensitivity variant: {variant}")
        item = config["sensitivity"][variant]
        for section, patch in item["spec"].items():
            config[section].update(patch)
        config["methods"] = [item["method"]]
    exp = config["experiment"]
    if (type(exp["case"]) is not int or exp["case"] not in range(1, 8)
            or type(exp["generating_R"]) is not int or exp["generating_R"] not in {10, 20, 50}
            or exp["var_a"] not in {0.5, 2.25} or exp["var_b"] not in {0, 0.1}
            or type(exp["repetitions"]) is not int or exp["repetitions"] != 5):
        raise ValueError("Configure one paper Case 1–7 grid setting and exactly five repetitions")
    sampler = config["sampler"]
    evaluation = config["evaluation"]
    for section, fields in ((sampler, ("chains", "workers", "draws", "thin", "checkpoint_every")),
                            (evaluation, ("pair_budget", "small_group_below", "em_iterations")),
                            (config["mfm"], ("component_support",))):
        if any(type(section[key]) is not int or section[key] <= 0 for key in fields):
            raise ValueError(f"Expected positive integers for {fields}")
    if sampler["chains"] < 2 or sampler["draws"] < 6:
        raise ValueError("Diagnostics require at least two chains with six draws per chain")
    for section, key in ((sampler, "warmup"), (evaluation, "conditional_draws")):
        if type(section[key]) is not int or section[key] < 0:
            raise ValueError(f"Expected a nonnegative integer for {key}")
    for name in ("residual", "intercept", "slope"):
        prior = config["lmm"][name]
        if not np.isfinite([prior["shape"], prior["scale"]]).all() or prior["shape"] <= 1 or prior["scale"] <= 0:
            raise ValueError(f"{name} inverse-gamma shape must exceed one and scale must be positive")
    positive = (config["lmm"]["beta_variance"], config["mfm"]["gamma"], evaluation["em_tolerance"],
                *config["dpm"].values(), config["ppmx"]["alpha"], *config["similarity"].values())
    if not all(np.isfinite(v) and v > 0 for v in positive) or not 0 < config["mfm"]["geometric_probability"] < 1:
        raise ValueError("Invalid prior, proposal, or evaluation parameters")
    methods = config["methods"]
    if not methods or len(set(methods)) != len(methods) or not set(methods) <= set(METHODS):
        raise ValueError("Select distinct methods from mfm, dpm and ppmx")
    return config


def evaluate(draws, x, y, reference, test, test_y, truth, method, cfg, job, directory, device):
    metrics = {**job.metadata(), "method": method}
    spec = cfg["evaluation"]
    assignment = Assignment(x, method, cfg)
    index, loss = representative(draws, spec["pair_budget"], job.metadata()["seed"] + 219)
    chosen = draws[index]
    ref_membership = assignment.membership(reference, chosen)
    test_membership = assignment.membership(test, chosen)
    ref_labels = ref_membership[:, :-1].argmax(axis=1)
    test_labels = test_membership[:, :-1].argmax(axis=1)
    expected, uncertainty, new_mass, new_group_trace = assignment.predict_new(test, draws)
    if not np.isfinite([expected, uncertainty, new_mass]).all():
        raise FloatingPointError("Invalid posterior prediction")
    metrics["test_MSPE_posterior_predictive"] = float(np.mean((test_y - expected)**2))
    metrics["mean_test_new_group_probability"] = float(new_mass.mean())
    metrics.update(K_representative=chosen.k, K_occ_mean=float(np.mean([d.k for d in draws])),
                   K_total_mean=float(np.mean([d.total for d in draws])) if method == "mfm" else None,
                   delta_K_vs_R=chosen.k - job.r, ari_train=adjusted_rand(truth, chosen.labels))

    predictions = dict(actual=test_y, posterior_mean=expected, posterior_sd=uncertainty,
                       new_group_probability=new_mass, representative_existing_group=test_labels)
    metrics.update(plug_in_hard_IMSPE=None, test_MSPE_hard_assignment=None)
    plugin_status = "computed"
    if chosen.k == 1:
        plugin_status = "one_occupied_group_variance_components_not_identifiable"
    else:
        try:
            evaluate_hard = make_evaluator(x, y, reference, job, spec, device, solver="torch")
            fitted, = evaluate_hard([(chosen.labels, ref_labels)])
            if fitted is None:
                raise ValueError("RGLL hard LMM fit returned no finite IMSPE")
            pred = fitted.predict(test, test_labels)
            value = float(np.mean((test_y - pred)**2))
            if not np.isfinite(value):
                raise ValueError("Non-finite hard prediction MSPE")
            predictions["hard_prediction"] = pred
            metrics.update(plug_in_hard_IMSPE=fitted.imspe, test_MSPE_hard_assignment=value)
            if not fitted.converged:
                plugin_status = "em_iteration_limit"
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            plugin_status = f"numerical_failure: {error}"

    pd.DataFrame(predictions).to_csv(directory / "test_predictions.csv", index=False)
    pd.DataFrame(new_group_trace).to_csv(directory / "new_group_trace.csv", index=False)
    write_npz(directory / "representative_partition.npz", train_labels=chosen.labels,
              reference_group_weights=ref_membership, test_group_weights=test_membership)
    atomic_json(directory / "partition_summary.json", dict(chain=chosen.chain, iteration=chosen.iteration,
                method="Dahl least squares over sampled unordered pairs (with replacement)",
                pair_budget=spec["pair_budget"], pair_loss=loss,
                representative_sizes=np.bincount(chosen.labels).tolist(),
                assignment="native_ppmx" if method == "ppmx" else "training_partition_post_fit_gate",
                hard_assignment="maximum probability among occupied groups"))

    count = min(spec["conditional_draws"], len(draws))
    conditional = []
    values = []
    if count:
        selected = np.linspace(0, len(draws) - 1, count, dtype=int) if count > 1 else [len(draws) // 2]
        for j in selected:
            d = draws[j]
            row = dict(chain=d.chain, iteration=d.iteration, k=d.k, score=None)
            try:
                group = assignment.membership(reference, d)[:, :-1].argmax(axis=1)
                row["score"] = conditional_imspe(x, reference, d.labels, group, d.variances, d.ve)
                values.append(row["score"])
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
                row["error"] = str(error)
            conditional.append(row)
    metrics["posterior_conditional_IMSPE_mean"] = float(np.mean(values)) if values else None
    metrics["plugin_status"] = plugin_status
    atomic_json(directory / "imspe_summary.json", dict(plug_in_hard_IMSPE=metrics["plug_in_hard_IMSPE"],
                hard_status=plugin_status, conditional_draws=conditional, draws_attempted=count,
                draws_successful=len(values), conditional_summary=summarize_numbers(values) if values else None))
    return metrics


def posterior_diagnostics(draws, chain_dirs, method, cfg, directory, x):
    if method == "ppmx":
        mean, scale = scale_x(x)
        standardized = (x - mean) / scale
        total_sq = np.square(standardized).sum()
        kappa = cfg["similarity"]["mean_precision"]
        variance = cfg["similarity"]["within_variance"]
    rows = []
    for d in draws:
        sizes = np.bincount(d.labels)
        row = dict(chain=d.chain, iteration=d.iteration, K_occ=d.k, var_e=d.ve, var_a=d.va,
                   min_size=int(sizes.min()), max_size=int(sizes.max()),
                   small_group_fraction=float(np.mean(sizes < cfg["evaluation"]["small_group_below"])),
                   intercept_rms=float(np.sqrt(np.mean(d.effects[:, 0]**2))))
        if d.vb > 0:
            row["var_b"] = d.vb
            row["slope_rms"] = float(np.sqrt(np.mean(d.effects[:, 1:]**2)))
        if method == "dpm":
            row["alpha"] = d.alpha
        if method == "mfm":
            row["K_total"] = d.total
        if method == "ppmx":
            sums = np.zeros((d.k, x.shape[1]))
            np.add.at(sums, d.labels, standardized)
            row["log_cohesion"] = float(d.k * np.log(d.alpha) + gammaln(sizes).sum())
            row["log_covariate_similarity"] = float(
                -0.5 * x.size * np.log(2 * np.pi * variance)
                + 0.5 * x.shape[1] * np.log(kappa / (kappa + sizes)).sum()
                - 0.5 / variance * (total_sq - (np.square(sums).sum(axis=1) / (kappa + sizes)).sum()))
        rows.append(row)
    trace = pd.DataFrame(rows)
    trace.to_csv(directory / "posterior_trace.csv", index=False)
    k_fields = [col for col in ("K_occ", "K_total") if col in trace]
    counts = []
    for field in k_fields:
        distribution = discrete_summary(trace[field].to_numpy())
        counts.extend(dict(quantity=field, k=k, posterior_probability=probability)
                      for k, probability in distribution["probability"].items())
    pd.DataFrame(counts).to_csv(directory / "posterior_K.csv", index=False)
    posterior = {field: discrete_summary(trace[field].to_numpy()) for field in k_fields}
    if method == "mfm":
        support = cfg["mfm"]["component_support"]
        mass = (1 - cfg["mfm"]["geometric_probability"]) ** np.arange(support)
        posterior["K_total_prior"] = dict(zip(range(1, support + 1), (mass / mass.sum()).tolist()))
        posterior["K_total_upper_boundary_mass"] = float(np.mean(trace.K_total == support))
    posterior.update(beta=summarize_numbers(np.array([d.beta for d in draws])),
                     variance={col: summarize_numbers(trace[col].to_numpy())
                               for col in ("var_e", "var_a", "var_b", "alpha") if col in trace},
                     effects={col: summarize_numbers(trace[col].to_numpy())
                              for col in ("intercept_rms", "slope_rms") if col in trace},
                     small_group_fraction=summarize_numbers(trace.small_group_fraction.to_numpy()))
    atomic_json(directory / "posterior_summary.json", posterior)
    scalar = {field: scalar_diagnostic(np.stack([g[field].to_numpy() for _, g in trace.groupby("chain")]))
              for field in ("K_occ", "K_total", "var_e", "var_a", "var_b", "alpha") if field in trace}
    checkpoints = []
    for chain_dir in chain_dirs:
        with np.load(chain_dir / "checkpoint.npz", allow_pickle=False) as f:
            iteration, alpha_accepted, rebuilds = map(int, f["iterations"])
            checkpoints.append(dict(sweeps=iteration, rank_one_rebuilds=rebuilds, elapsed=float(f["elapsed"].item()),
                                    peak_resident_bytes=int(f["peak_memory"].item()) or None,
                                    alpha_acceptance=alpha_accepted / iteration if method == "dpm" else None))
    atomic_json(directory / "sampler_diagnostics.json", dict(
        scalar=scalar, chains=checkpoints, retained_draws=len(draws),
        ppmx_similarity=cfg["similarity"] if method == "ppmx" else None,
        min_chain_ess=min((v["ess"] for v in scalar.values() if v["ess"] is not None), default=None)))
    return checkpoints


def run_method(data, method, cfg, job, directory, device, executor):
    started = perf_counter()
    directory.mkdir(parents=True, exist_ok=True)
    done = directory / "metrics.json"
    x, reference, test, y, test_y, truth, _ = data
    if len(x) != 7500 or x.shape[1] != 50 or len(reference) != 2500 or len(test) != 2500:
        raise ValueError("Unexpected dimensions from the original RGLL case generator")
    sampler = cfg["sampler"]
    dirs = [directory / f"chain_{i + 1}" for i in range(sampler["chains"])]
    # Seed drives chains; R and the simulation truth labels never enter fit_chain.
    seed = job.metadata()["seed"] + 9300 + METHODS.index(method) * 1000
    seeds = [seed + i * 97 for i in range(len(dirs))]
    digest = sha256()
    for array in (x, y, reference, test, test_y):
        digest.update(np.ascontiguousarray(array).tobytes())
    session = directory / "resolved_config.json"
    resolved = {"config": cfg, "case": job.metadata(), "method": method, "seeds": seeds,
                "data_fingerprint": digest.hexdigest(), "evaluation_device": device}
    if session.exists() and json.loads(session.read_text(encoding="utf-8")) != resolved:
        raise ValueError("Checkpoint configuration or input data changed")
    if done.exists():
        return json.loads(done.read_text(encoding="utf-8"))
    atomic_json(session, resolved)
    args = [(x, y, method, bool(job.var_b), cfg, seeds[i], slot, bool(i % 2)) for i, slot in enumerate(dirs)]
    if executor is None:
        for arg in args:
            fit_chain(*arg)
    else:
        futures = [executor.submit(fit_chain, *arg) for arg in args]
        for future in futures:
            future.result()
    draws = load_draws(dirs)
    if len(draws) != len(dirs) * sampler["draws"]:
        raise RuntimeError("Unexpected saved draw count")
    time_before_eval = perf_counter()
    metrics = evaluate(draws, x, y, reference, test, test_y, truth, method, cfg, job, directory, device)
    checkpoints = posterior_diagnostics(draws, dirs, method, cfg, directory, x)
    evaluation_time = perf_counter() - time_before_eval
    metrics.update(draws=len(draws), evaluation_wall_s=evaluation_time,
                   chain_process_seconds=sum(c["elapsed"] for c in checkpoints),
                   current_invocation_wall_s=perf_counter() - started, status="computed")
    atomic_json(directory / "runtime.json", dict(
        resumed_chain_times_included=True, **{k: metrics[k] for k in
        ("evaluation_wall_s", "chain_process_seconds", "current_invocation_wall_s")},
        peak_chain_resident_bytes=max((c["peak_resident_bytes"] or 0 for c in checkpoints), default=0) or None,
        rank_one_rebuilds=sum(c["rank_one_rebuilds"] for c in checkpoints),
        peak_evaluation_resident_bytes=peak_memory_bytes(), sampler_device="cpu", evaluation_device=device))
    atomic_json(done, metrics)
    return metrics


def write_comparisons(rows, output):
    d = pd.DataFrame(rows)
    d.to_csv(output / "replications.csv", index=False)
    table = []
    fields = ("plug_in_hard_IMSPE", "posterior_conditional_IMSPE_mean", "test_MSPE_hard_assignment",
              "test_MSPE_posterior_predictive", "K_occ_mean", "K_total_mean", "current_invocation_wall_s")

    def record(name, metric, values):
        values = pd.to_numeric(values, errors="coerce").dropna()
        n = len(values)
        mean = float(values.mean()) if n else None
        margin = float(t.ppf(0.975, n - 1) * values.sem()) if n > 1 else None
        table.append(dict(comparison=name, metric=metric, n=n, mean=mean,
                          low95=mean - margin if margin is not None else None,
                          high95=mean + margin if margin is not None else None))

    for method, group in d.groupby("method"):
        for metric in fields:
            if metric in group:
                record(method, metric, group[metric])
    for metric in fields[:4]:
        if metric in d:
            m = d.pivot(index="replication", columns="method", values=metric).apply(pd.to_numeric, errors="coerce")
            for i, method in enumerate(m.columns):
                for other in m.columns[:i]:
                    record(f"{method} - {other}", metric, m[method] - m[other])
    pd.DataFrame(table).to_csv(output / "comparison.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.pilot.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", shutil.which("Rscript")))
    parser.add_argument("--rgll", type=Path, default=DEFAULT_RGLL)
    parser.add_argument("--variant")
    parser.add_argument("--device", choices=("cpu", "cuda", "auto"), default="auto")
    args = parser.parse_args()
    if not args.rscript:
        parser.error("Set RSCRIPT or pass --rscript")
    required_rgll = ("ROG_method_comparison.R", "ROG_aligned_fixed.R", "lmm_fast_aligned.cpp")
    missing = [name for name in required_rgll if not (args.rgll / name).is_file()]
    if missing:
        parser.error(f"DGP_SPEC_INCOMPLETE: missing in {args.rgll}: {', '.join(missing)}")
    if args.device == "cuda" and not torch.cuda.is_available():
        parser.error("CUDA is unavailable")
    device = ("cuda" if torch.cuda.is_available() else "cpu") if args.device == "auto" else args.device
    cfg = load_config(args.config, args.variant)
    torch.set_num_threads(2)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    signature = output / "run_config.json"
    run_spec = dict(config=cfg, variant=args.variant, rgll=str(args.rgll.resolve()), evaluation_device=device)
    if signature.exists() and json.loads(signature.read_text(encoding="utf-8")) != run_spec:
        parser.error("Output contains a different run configuration")
    atomic_json(signature, run_spec)
    rows = []
    exp = cfg["experiment"]
    workers = min(cfg["sampler"]["workers"], cfg["sampler"]["chains"])
    # A serial run is useful when checking one process's resource use.
    from contextlib import nullcontext
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else nullcontext()
    with pool as executor, RGLL(args.rscript, args.rgll.resolve()) as rgll:
        for rep in range(1, exp["repetitions"] + 1):
            job = Job(exp["case"], exp["generating_R"], exp["var_a"], exp["var_b"], rep)
            try:
                data = rgll.data(job)
            except RuntimeError as exc:
                message = str(exc)
                if "not found" in message or "Unknown dist_x" in message:
                    raise RuntimeError(f"DGP_SPEC_INCOMPLETE: Case {job.case} generator: {message}") from exc
                raise
            for method in cfg["methods"]:
                directory = output / job.key / method
                try:
                    row = run_method(data, method, cfg, job, directory, device, executor)
                except (RuntimeError, ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                    row = {**job.metadata(), "method": method, "status": "failed", "error": str(exc)}
                    atomic_json(directory / "failure.json", row)
                rows.append(row)
                print(f"replication={rep} method={method} status={row['status']}", flush=True)
            write_comparisons(rows, output)
    if any(row["status"] == "failed" for row in rows):
        raise SystemExit("One or more fits failed. See failure.json and resume with the same arguments.")


if __name__ == "__main__":
    main()
