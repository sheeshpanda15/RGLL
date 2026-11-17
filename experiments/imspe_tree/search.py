"""Tree rules and the three adaptive IMSPE searches."""

from collections import defaultdict
from dataclasses import dataclass
from hashlib import blake2b
from math import exp
from time import perf_counter

import numpy as np


@dataclass(frozen=True)
class Split:
    feature: int
    threshold: float
    left: "Split | None" = None
    right: "Split | None" = None

    @classmethod
    def from_dict(cls, value):
        if value is None:
            return None
        return cls(value["feature"], value["threshold"],
                   cls.from_dict(value["left"]), cls.from_dict(value["right"]))


def walk(tree, path=()):
    yield path, tree
    if tree is not None:
        yield from walk(tree.left, (*path, 0))
        yield from walk(tree.right, (*path, 1))


def edit(tree, path, replacement):
    if not path:
        return replacement
    if path[0] == 0:
        return Split(tree.feature, tree.threshold, edit(tree.left, path[1:], replacement), tree.right)
    return Split(tree.feature, tree.threshold, tree.left, edit(tree.right, path[1:], replacement))


def parts(tree, x, rows=None, path=()):
    rows = np.arange(len(x)) if rows is None else rows
    yield path, rows
    if tree is not None:
        left = x[rows, tree.feature] <= tree.threshold
        yield from parts(tree.left, x, rows[left], (*path, 0))
        yield from parts(tree.right, x, rows[~left], (*path, 1))


def predict_group(tree, x):
    locations = dict(parts(tree, x))
    labels = np.empty(len(x), dtype=np.int32)
    leaves = [path for path, node in walk(tree) if node is None]
    for group, path in enumerate(leaves):
        labels[locations[path]] = group
    return labels


def canonical(train, ref, k):
    first = np.full(k, len(train))
    np.minimum.at(first, train, np.arange(len(train)))
    remap = first.argsort().argsort().astype(np.int32)
    return remap[train], remap[ref], remap


def state_id(tree):
    return blake2b(repr(tree).encode(), digest_size=12).hexdigest()


@dataclass
class Candidate:
    tree: Split
    parent: "Candidate | None"
    move: str
    k: int = 0
    depth: int = 0
    counts: list | None = None
    remap: np.ndarray | None = None
    fit: object | None = None
    key: bytes | None = None
    row: dict | None = None


class Search:
    def __init__(self, x, ref, evaluate, config, seed):
        self.x, self.ref, self.evaluate, self.config = x, ref, evaluate, config
        self.rng = np.random.default_rng(seed)
        self.history = []
        self.cache = {}
        self.best = None
        self.evaluations = 0

    def assess(self, candidates, stage, cap, step):
        pending = defaultdict(list)
        evaluated_now = {}
        recorded = []
        for c in candidates:
            leaves = [path for path, node in walk(c.tree) if node is None]
            c.k, c.depth = len(leaves), max(map(len, leaves))
            labels = predict_group(c.tree, self.x)
            c.counts = list(map(int, np.bincount(labels, minlength=c.k)))
            status = "pending"
            if min(c.counts) < self.config["min_leaf"]:
                status = "minimum_leaf_size"
            else:
                train, ref, c.remap = canonical(labels, predict_group(c.tree, self.ref), c.k)
                # Both partitions matter: equal training labels can induce
                # different reference integration regions and IMSPE values.
                c.key = blake2b(train.tobytes() + ref.tobytes(), digest_size=20).digest()
                if c.key not in self.cache and c.key not in evaluated_now:
                    if self.evaluations >= cap:
                        status = "evaluation_budget"
                    else:
                        pending[c.k].append((c, train, ref))
                        self.evaluations += 1
                        evaluated_now[c.key] = None
                else:
                    status = "cache_hit"
            c.row = dict(step=step, stage=stage, action=c.move, candidate_id=state_id(c.tree),
                         parent_id=state_id(c.parent.tree) if c.parent else None,
                         k=c.k, depth=c.depth, leaf_sizes=c.counts,
                         evaluations=self.evaluations, status=status, accepted=False,
                         partition_hash=c.key.hex() if c.key else None)
            self.history.append(c.row)
            recorded.append(c)
        for group in pending.values():
            for start in range(0, len(group), self.config["fit_batch"]):
                batch = group[start:start + self.config["fit_batch"]]
                t0 = perf_counter()
                try:
                    fits = self.evaluate([(train, ref) for _, train, ref in batch])
                except (ValueError, RuntimeError) as error:
                    fits = [None] * len(batch)
                    for c, _, _ in batch:
                        c.row["error"] = str(error)
                elapsed = perf_counter() - t0
                for (c, _, _), fit in zip(batch, fits, strict=True):
                    self.cache[c.key] = fit
                    c.row["batch_seconds"] = elapsed
                    c.row["batch_size"] = len(batch)
                    c.row["eval_seconds_share"] = elapsed / len(batch)
        valid = []
        for c in recorded:
            if c.row["status"] in {"pending", "cache_hit"}:
                c.fit = self.cache[c.key]
                c.row["fit_status"] = ("failed" if c.fit is None else
                                       "converged" if c.fit.converged else "em_iteration_limit")
                if c.row["status"] == "pending":
                    c.row["status"] = "evaluated"
                if c.fit is not None:
                    c.row["imspe"] = c.fit.imspe
                    c.row["score"] = c.fit.imspe
                    if self.best is None or c.fit.imspe < self.best.fit.imspe:
                        self.best = c
                    valid.append(c)
            c.row["best_so_far"] = self.best.fit.imspe if self.best else None
        return valid

    def propose(self, tree, move, count, rng, parent=None):
        nodes = list(walk(tree))
        leaves = [path for path, node in nodes if node is None]
        internal = [(path, node) for path, node in nodes if node is not None]
        if move == "prune_subtree":
            return [Candidate(edit(tree, path, None), parent, move)
                    for path, _ in internal if path]  # K=1 has no group variance estimate.

        paths = leaves if move == "split_leaf" else [path for path, _ in internal]
        locations = dict(parts(tree, self.x))
        minimum = self.config["min_leaf"]
        paths = [path for path in paths if len(locations[path]) >= 2 * minimum]
        if not paths:
            return []
        by_path = dict(nodes)
        out = {}
        for _ in range(10 * count):
            path = paths[int(rng.integers(len(paths)))]
            feature = int(rng.integers(self.x.shape[1]))
            values = self.x[locations[path], feature]
            rank = int(rng.integers(minimum, len(values) - minimum + 1))
            lower, upper = np.partition(values, [rank - 1, rank])[rank - 1:rank + 1]
            if lower >= upper:
                continue
            old = by_path[path]
            rule = Split(feature, float(lower + (upper - lower) / 2),
                         old.left if old else None, old.right if old else None)
            changed = edit(tree, path, rule)
            if changed != tree and changed not in out:
                out[changed] = Candidate(changed, parent, move)
            if len(out) == count:
                break
        return list(out.values())

    @staticmethod
    def retain(candidates, width):
        retained = []
        seen = set()
        for candidate in sorted(candidates, key=lambda c: c.fit.imspe):
            if candidate.key in seen:
                continue
            retained.append(candidate)
            seen.add(candidate.key)
            candidate.row["accepted"] = True
            if len(retained) == width:
                break
        return retained

    def run(self, method, initial_seed):
        budget = self.config["max_evaluations"]
        initial = self.propose(None, "split_leaf", self.config["initial_candidates"],
                               np.random.default_rng(initial_seed))
        valid = self.assess(initial, "initial", budget, 0)
        if not valid:
            return None, "no_feasible_initial_tree"
        if method == "GrowThenPrune":
            reason = self.grow_prune(valid, budget)
        elif method == "TreeAnnealing":
            reason = self.anneal(valid, budget)
        elif method == "BeamSearch":
            reason = self.beam(valid, budget)
        else:
            raise ValueError(f"Unknown search method {method}")
        return self.best, reason

    def grow_prune(self, initial, budget):
        current, = self.retain(initial, 1)
        count = self.config["candidates_per_parent"]
        grow_cap = self.config["grow_evaluations"]
        step = 0
        while self.evaluations < grow_cap:
            step += 1
            new = self.propose(current.tree, "split_leaf", count, self.rng, current)
            valid = self.assess(new, "grow", grow_cap, step)
            if not valid:
                break
            current, = self.retain(valid, 1)
        while self.evaluations < budget:
            step += 1
            new = self.propose(current.tree, "prune_subtree", count, self.rng, current)
            valid = self.assess(new, "prune", budget, step)
            if not valid:
                break
            current, = self.retain(valid, 1)
        return "evaluation_budget" if self.evaluations >= budget else "pruning_path_exhausted"

    def anneal(self, initial, budget):
        current, = self.retain(initial, 1)
        temperature = max(current.fit.imspe * self.config["anneal_temperature_fraction"], 1e-10)
        for step in range(1, self.config["anneal_proposals"] + 1):
            if self.evaluations >= budget:
                return "evaluation_budget"
            move = self.rng.choice(["split_leaf", "prune_subtree", "change_split"],
                                   p=self.config["anneal_move_probabilities"])
            new = self.propose(current.tree, move, 1, self.rng, current)
            if new:
                c = new[int(self.rng.integers(len(new)))]
                valid = self.assess([c], "anneal", budget, step)
                probability = exp(-max(0.0, c.fit.imspe - current.fit.imspe) / temperature) if valid else 0.0
                accepted = bool(valid) and self.rng.random() < probability
                if accepted:
                    current = c
                c.row.update(temperature=temperature, probability=probability, accepted=accepted,
                             current_score=current.fit.imspe)
            else:
                self.history.append(dict(step=step, stage="anneal", action=move, status="no_available_move",
                                         temperature=temperature, probability=0.0, accepted=False,
                                         parent_id=state_id(current.tree), current_score=current.fit.imspe,
                                         best_so_far=self.best.fit.imspe, evaluations=self.evaluations))
            temperature *= self.config["anneal_cooling"]
        return "proposal_budget"

    def beam(self, initial, budget):
        width = self.config["beam_width"]
        frontier = self.retain(initial, width)
        step = 0
        while self.evaluations < budget:
            step += 1
            per_parent = [self.propose(parent.tree, "split_leaf", self.config["candidates_per_parent"],
                                      self.rng, parent) for parent in frontier]
            # Interleave parents so the budget does not favor the first tree.
            new = [items[j] for j in range(max(map(len, per_parent))) for items in per_parent if j < len(items)]
            valid = self.assess(new, "beam", budget, step)
            if not valid:
                break
            frontier = self.retain(valid, width)
        return "evaluation_budget" if self.evaluations >= budget else "frontier_exhausted"
