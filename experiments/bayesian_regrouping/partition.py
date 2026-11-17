"""Collapsed partition weights and covariate-only predictive assignment."""

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln, logsumexp


@dataclass(frozen=True)
class PartitionLaw:
    method: str
    n: int
    spec: dict

    def __post_init__(self):
        if self.method not in {"mfm", "dpm", "ppmx"}:
            raise ValueError(f"Unknown model: {self.method}")
        if self.method == "mfm":
            support = np.arange(1, self.spec["component_support"] + 1)
            gamma = self.spec["gamma"]
            probability = self.spec["geometric_probability"]
            log_prior = (support - 1) * np.log1p(-probability)
            log_prior -= logsumexp(log_prior)
            possible_t = np.arange(support.size + 1)[:, None]
            log_falling = gammaln(support + 1) - gammaln(support - possible_t + 1)
            coefficients = log_prior + log_falling
            terms = {size: coefficients + gammaln(gamma * support) - gammaln(gamma * support + size)
                     for size in (self.n, self.n + 1)}
            object.__setattr__(self, "support", support)
            object.__setattr__(self, "terms", terms)
            object.__setattr__(self, "log_v", {size: logsumexp(values, axis=1) for size, values in terms.items()})

    def weights(self, counts, alpha, predictive=False):
        if self.method != "mfm":
            return np.r_[np.log(counts), np.log(alpha)]
        gamma = self.spec["gamma"]
        k = len(counts)
        log_v = self.log_v[self.n + int(predictive)]
        new = np.log(gamma) + log_v[k + 1] - log_v[k] if k < len(self.support) else -np.inf
        return np.r_[np.log(counts + gamma), new]

    def total_distribution(self, occupied):
        log_p = self.terms[self.n][occupied]
        return self.support, np.exp(log_p - logsumexp(log_p))

    def sample_alpha(self, alpha, occupied, rng):
        if self.method != "dpm":
            return alpha, False
        shape = self.spec["alpha_prior_shape"]
        rate = self.spec["alpha_prior_rate"]

        def log_target(log_alpha):
            value = np.exp(log_alpha)
            # CRP EPPF times a gamma prior and the log-scale Jacobian.
            return (occupied + shape) * log_alpha + gammaln(value) - gammaln(value + self.n) - rate * value

        current = np.log(alpha)
        proposal = current + rng.normal(scale=self.spec["alpha_log_step"])
        accept = np.log(rng.random()) < log_target(proposal) - log_target(current)
        return (float(np.exp(proposal)), True) if accept else (alpha, False)


def similarity_parameters(counts, sums, spec):
    """Predictive x density from a diagonal Gaussian with integrated group means."""
    kappa = spec["mean_precision"]
    precision = kappa + counts[:, None]
    mean = sums / precision  # Training-standardized X, with zero prior mean.
    variance = spec["within_variance"] * (1 + 1 / precision)
    return mean, variance


def covariate_log_density(x, counts, sums, spec):
    counts = np.r_[counts, 0]
    sums = np.vstack((sums, np.zeros(sums.shape[1])))
    means, variance = similarity_parameters(counts, sums, spec)
    squared = np.square(np.asarray(x)[..., None, :] - means).sum(axis=-1)
    return -0.5 * (sums.shape[1] * np.log(2 * np.pi * variance[:, 0]) + squared / variance[:, 0])


def predict_membership(x_standardized, train_standardized, labels, law, alpha, similarity):
    """Per-draw probabilities, last column is the new-group event."""
    k = int(labels.max()) + 1
    counts = np.bincount(labels, minlength=k)
    sums = np.zeros((k, train_standardized.shape[1]))
    np.add.at(sums, labels, train_standardized)
    weights = law.weights(counts, alpha, predictive=True)
    # Native PPMx predictive; the same density is a post-fit gate for DPM/MFM.
    result = []
    batch = min(128, max(1, 2_000_000 // ((k + 1) * sums.shape[1])))
    for start in range(0, len(x_standardized), batch):
        x = x_standardized[start:start + batch]
        score = weights + covariate_log_density(x, counts, sums, similarity)
        result.append(np.exp(score - logsumexp(score, axis=1, keepdims=True)))
    return np.concatenate(result, axis=0) if result else np.empty((0, k + 1))
