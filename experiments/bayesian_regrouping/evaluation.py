"""Label-invariant summaries, posterior prediction and hard Eq. (10)."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from latent_lmm import design, scale_x
from partition import PartitionLaw, predict_membership


@dataclass(frozen=True)
class Draw:
    labels: np.ndarray
    beta: np.ndarray
    effects: np.ndarray
    ve: float
    va: float
    vb: float
    alpha: float
    total: int
    iteration: int
    chain: int

    @property
    def k(self):
        return len(self.effects)

    @property
    def variances(self):
        return np.r_[self.va, np.repeat(self.vb, self.effects.shape[1] - 1)]


def load_draws(directories):
    draws = []
    for chain, directory in enumerate(directories, start=1):
        for path in sorted(directory.glob("draw_*.npz")):
            with np.load(path, allow_pickle=False) as f:
                draws.append(Draw(f["labels"], f["beta"], f["effects"], *map(float, f["hyper"]),
                                  int(f["total_components"].item()), int(f["iteration"].item()), chain))
    return draws


def representative(draws, pair_budget, seed):
    """Dahl decision restricted to saved draws, using sampled co-clustering pairs."""
    n = len(draws[0].labels)
    rng = np.random.default_rng(seed)
    a = rng.integers(n, size=pair_budget)
    b = (a + rng.integers(1, n, size=pair_budget)) % n
    average = np.zeros(pair_budget)
    for draw in draws:
        average += draw.labels[a] == draw.labels[b]
    average /= len(draws)
    distance = [np.square((draw.labels[a] == draw.labels[b]) - average).mean() for draw in draws]
    index = int(np.argmin(distance))
    return index, float(distance[index])


class Assignment:
    """Per-draw group probabilities conditional on new X, never conditional on new y."""

    def __init__(self, train_x, method, config):
        self.center, self.scale = scale_x(train_x)
        self.train = (train_x - self.center) / self.scale
        self.law = PartitionLaw(method, len(train_x), config[method])
        self.similarity = config["similarity"]

    def membership(self, x_new, draw):
        return predict_membership((x_new - self.center) / self.scale, self.train, draw.labels,
                                  self.law, draw.alpha, self.similarity)

    def predict_new(self, x_new, draws):
        f = design(x_new)
        first = np.zeros(len(f))
        second = np.zeros(len(f))
        new_mass = np.zeros(len(f))
        new_group_trace = []
        for draw in draws:
            w = self.membership(x_new, draw)
            z = f[:, :draw.effects.shape[1]]
            fixed = f @ draw.beta
            group_mean = fixed[:, None] + z @ draw.effects.T
            new_variance = np.square(z) @ draw.variances
            first += np.einsum("ij,ij->i", w[:, :-1], group_mean) + w[:, -1] * fixed
            second += (np.einsum("ij,ij->i", w[:, :-1], group_mean**2)
                       + w[:, -1] * (fixed**2 + new_variance) + draw.ve)
            new_mass += w[:, -1]
            new_group_trace.append(dict(chain=draw.chain, iteration=draw.iteration,
                                        mean_new_group_probability=float(w[:, -1].mean()) if len(w) else None))
        first /= len(draws)
        return first, np.sqrt(np.maximum(second / len(draws) - first**2, 0)), new_mass / len(draws), new_group_trace


def conditional_imspe(train_x, reference_x, train_labels, reference_labels, variances, ve):
    """Eq. (10) at a supplied variance draw, conditional on hard existing groups."""
    f, f_ref = design(train_x), design(reference_x)
    q = f.shape[1]
    rank = len(variances)
    k = int(train_labels.max()) + 1
    information = f.T @ f / ve
    fixed_weight = np.zeros((q, q))
    random_term = 0.0
    for group in range(k):
        fg = f[train_labels == group]
        fr = f_ref[reference_labels == group]
        w = fr.T @ fr / len(f_ref)
        z = fg[:, :rank]
        b = fg.T @ z / ve
        factor = cho_factor(np.diag(1 / variances) + z.T @ z / ve, lower=True)
        influence = cho_solve(factor, b.T).T
        information -= influence @ b.T
        transform = np.eye(q)
        transform[:, :rank] -= influence
        fixed_weight += transform @ w @ transform.T
        random_term += np.trace(cho_solve(factor, w[:rank, :rank]))
    fixed_term = np.trace(cho_solve(cho_factor((information + information.T) / 2, lower=True), fixed_weight))
    score = float(ve + fixed_term + random_term)
    if not np.isfinite(score) or score <= 0:
        raise FloatingPointError("Invalid hard IMSPE at supplied variance components")
    return score


def summarize_numbers(values):
    quantiles = np.quantile(values, [0.025, 0.5, 0.975], axis=0)
    return dict(mean=np.asarray(np.mean(values, axis=0)).tolist(), median=quantiles[1].tolist(),
                low95=quantiles[0].tolist(), high95=quantiles[2].tolist())


def discrete_summary(values):
    unique, counts = np.unique(values, return_counts=True)
    return {**summarize_numbers(values), "mode": int(unique[np.argmax(counts)]),
            "probability": dict(zip(map(int, unique), (counts / counts.sum()).tolist()))}


def scalar_diagnostic(chains):
    """Classical split R-hat and an initial-positive-sequence ESS estimate."""
    chains = np.asarray(chains, dtype=float)
    half = chains.shape[1] // 2
    if chains.shape[0] < 2 or half < 3:
        return dict(split_rhat=None, ess=None, status="insufficient_draws")
    x = np.concatenate((chains[:, :half], chains[:, -half:]), axis=0)
    within = x.var(axis=1, ddof=1).mean()
    if within <= 0:
        return dict(split_rhat=None, ess=None, status="constant_or_stuck_chains")
    var_plus = (half - 1) / half * within + x.mean(axis=1).var(ddof=1)
    centered = x - x.mean(axis=1, keepdims=True)
    rho = [1.0] + [1 - (within - (centered[:, :-lag] * centered[:, lag:]).sum(axis=1).mean() / half) / var_plus
                   for lag in range(1, half)]
    pairs = []
    for i in range(0, len(rho) - 1, 2):
        total = rho[i] + rho[i + 1]
        if total <= 0:
            break
        pairs.append(min(total, pairs[-1]) if pairs else total)
    ess = x.size / max(1.0, -1 + 2 * sum(pairs))
    rhat = float(max(1.0, np.sqrt(var_plus / within)))
    return dict(split_rhat=rhat, ess=float(ess), status="review" if rhat > 1.01 else "computed")
