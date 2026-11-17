"""Partially collapsed single-site Gibbs for the shared RI/RS model."""

import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve, cholesky, solve_triangular
from scipy.special import logsumexp

from partition import PartitionLaw, covariate_log_density


def design(x):
    return np.column_stack((np.ones(len(x)), x))


def scale_x(x):
    mean = x.mean(axis=0)
    sd = x.std(axis=0)
    return mean, np.where(sd > 0, sd, 1.0)


def gaussian_conditional(covariates, residual, variances, ve):
    r = len(variances)
    precision = np.diag(1 / variances) + covariates.T @ covariates / ve
    factor = cho_factor(precision, lower=True)
    cov = cho_solve(factor, np.eye(r))
    mean = cho_solve(factor, covariates.T @ residual / ve)
    return mean, cov


def update_gaussian(mean, cov, g, residual, ve, sign):
    """Add/remove one normal likelihood factor via a rank-one update."""
    projected = cov @ g
    denominator = ve + sign * (g @ projected)
    if not np.isfinite(denominator) or denominator <= 1e-12 * ve:
        raise np.linalg.LinAlgError("Unstable Gaussian rank-one update")
    gain = sign * projected / denominator
    return mean + gain * (residual - g @ mean), cov - np.outer(gain, projected)


def sample_normal_precision(precision, information, rng):
    factor = cholesky(precision, lower=True)
    mean = cho_solve((factor, True), information)
    return mean + solve_triangular(factor.T, rng.normal(size=len(mean)), lower=False)


def inverse_gamma(spec, n, squares, rng):
    return (spec["scale"] + 0.5 * squares) / rng.gamma(spec["shape"] + 0.5 * n)


def write_npz(path, **values):
    temporary = path.with_suffix(".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **values)
    temporary.replace(path)


def peak_memory_bytes():
    if sys.platform != "win32":
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(rss if sys.platform == "darwin" else 1024 * rss)

    import ctypes
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in ("peak_ws", "ws", "peak_paged", "paged", "peak_nonpaged",
                                                "nonpaged", "pagefile", "peak_pagefile")]

    kernel = ctypes.windll.kernel32
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
    data = Counters()
    data.cb = ctypes.sizeof(data)
    return int(data.peak_ws) if query(kernel.GetCurrentProcess(), ctypes.byref(data), data.cb) else None


def gibbs_sweep(f, z, x_std, y, labels, beta, variances, ve, law, alpha, similarity, rng):
    """Update g after marginalizing u, conditional on beta, D and residual variance."""
    n, r = z.shape
    residual = y - f @ beta
    base = np.diag(variances)
    counts = np.bincount(labels, minlength=n)
    active = set(np.flatnonzero(counts).tolist())
    free = sorted(set(range(n)).difference(active), reverse=True)
    means = np.empty((n, r))
    covs = np.empty((n, r, r))
    sums = np.zeros((n, x_std.shape[1])) if law.method == "ppmx" else None
    rebuilds = 0
    for group in active:
        mask = labels == group
        means[group], covs[group] = gaussian_conditional(z[mask], residual[mask], variances, ve)
        if sums is not None:
            sums[group] = x_std[mask].sum(axis=0)

    for i in rng.permutation(n):
        old = labels[i]
        g = z[i]
        counts[old] -= 1
        labels[i] = -1
        if counts[old]:
            try:
                means[old], covs[old] = update_gaussian(means[old], covs[old], g, residual[i], ve, -1)
            except np.linalg.LinAlgError:
                # Same conditional, rebuilt by Cholesky when downdating is ill-conditioned.
                mask = labels == old
                means[old], covs[old] = gaussian_conditional(z[mask], residual[mask], variances, ve)
                rebuilds += 1
            if sums is not None:
                sums[old] -= x_std[i]
        else:
            active.remove(old)
            free.append(old)
        ids = np.array(sorted(active), dtype=np.intp)
        expected = np.r_[means[ids] @ g, 0.0]
        variance = np.r_[ve + np.einsum("i,kij,j->k", g, covs[ids], g), ve + np.dot(variances, g**2)]
        if np.any(variance <= 0) or not np.isfinite(variance).all():
            raise FloatingPointError("Invalid collapsed response predictive variance")
        weights = law.weights(counts[ids], alpha) - 0.5 * (np.log(variance) + (residual[i] - expected)**2 / variance)
        if sums is not None:
            weights += covariate_log_density(x_std[i], counts[ids], sums[ids], similarity)
        probs = np.exp(weights - logsumexp(weights))
        chosen = rng.choice(len(probs), p=probs)
        if chosen == len(ids):
            group = free.pop()
            active.add(group)
            means[group], covs[group] = 0.0, base
            if sums is not None:
                sums[group] = 0.0
        else:
            group = ids[chosen]
        means[group], covs[group] = update_gaussian(means[group], covs[group], g, residual[i], ve, 1)
        counts[group] += 1
        labels[i] = group
        if sums is not None:
            sums[group] += x_std[i]

    occupied, labels = np.unique(labels, return_inverse=True)
    effects = np.empty((len(occupied), r))
    for j in range(len(occupied)):
        # Recompute the effect conditional to avoid accumulating roundoff in draws.
        mask = labels == j
        gj = z[mask]
        precision = np.diag(1 / variances) + gj.T @ gj / ve
        effects[j] = sample_normal_precision(precision, gj.T @ residual[mask] / ve, rng)
    return labels.astype(np.int32), effects, rebuilds


def fit_chain(x, y, method, random_slope, config, seed, directory, prior_start=False):
    """Run/resume one independent chain. Checkpoints contain the full RNG state."""
    started = perf_counter()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    f = design(x)
    z = f if random_slope else f[:, :1]
    mean, scale = scale_x(x)
    x_std = (x - mean) / scale
    law = PartitionLaw(method, len(y), config[method])
    sampler, prior = config["sampler"], config["lmm"]
    checkpoint = directory / "checkpoint.npz"
    rng = np.random.default_rng(seed)
    labels = np.zeros(len(y), dtype=np.int32)
    ve = prior["residual"]["scale"] / (prior["residual"]["shape"] - 1)
    va = prior["intercept"]["scale"] / (prior["intercept"]["shape"] - 1)
    vb = prior["slope"]["scale"] / (prior["slope"]["shape"] - 1) if random_slope else 0.0
    alpha = config[method].get("alpha", 1.0)
    ftf = f.T @ f
    beta_precision = np.eye(f.shape[1]) / prior["beta_variance"]
    beta = cho_solve(cho_factor(beta_precision + ftf / ve, lower=True), f.T @ y / ve)
    if prior_start and not checkpoint.exists():
        if method == "mfm":
            mass = (1 - config[method]["geometric_probability"]) ** (law.support - 1)
            components = rng.choice(law.support, p=mass / mass.sum())
            weights = rng.dirichlet(np.full(components, config[method]["gamma"]))
            labels = np.unique(rng.choice(components, size=len(y), p=weights), return_inverse=True)[1].astype(np.int32)
        else:
            counts = [1]
            for i in range(1, len(y)):
                group = rng.choice(len(counts) + 1, p=np.r_[counts, alpha] / (i + alpha))
                labels[i] = group
                if group == len(counts):
                    counts.append(0)
                counts[group] += 1
    completed, accepted, rebuild_count, elapsed, memory = 0, 0, 0, 0.0, 0
    if checkpoint.exists():
        with np.load(checkpoint, allow_pickle=False) as saved:
            labels, beta = saved["labels"], saved["beta"]
            ve, va, vb, alpha = map(float, saved["hyper"])
            completed, accepted, rebuild_count = map(int, saved["iterations"])
            elapsed = float(saved["elapsed"].item())
            memory = int(saved["peak_memory"].item())
            rng.bit_generator.state = json.loads(saved["rng"].item())

    total = sampler["warmup"] + sampler["draws"] * sampler["thin"]
    for iteration in range(completed + 1, total + 1):
        variances = np.r_[va, np.repeat(vb, z.shape[1] - 1)]
        labels, effects, rebuilds = gibbs_sweep(f, z, x_std, y, labels, beta, variances, ve, law, alpha,
                                               config["similarity"], rng)
        rebuild_count += rebuilds
        random = np.einsum("ij,ij->i", z, effects[labels])
        beta = sample_normal_precision(beta_precision + ftf / ve, f.T @ (y - random) / ve, rng)
        residual = y - f @ beta - random
        ve = inverse_gamma(prior["residual"], len(y), residual @ residual, rng)
        va = inverse_gamma(prior["intercept"], len(effects), np.square(effects[:, 0]).sum(), rng)
        if random_slope:
            vb = inverse_gamma(prior["slope"], effects[:, 1:].size, np.square(effects[:, 1:]).sum(), rng)
        alpha, changed = law.sample_alpha(alpha, len(effects), rng)
        accepted += changed
        if not np.isfinite([ve, va, vb, alpha]).all() or min(ve, va, alpha) <= 0 or (random_slope and vb <= 0):
            raise FloatingPointError("Invalid variance or concentration draw")
        if not np.isfinite(beta).all() or not np.isfinite(effects).all():
            raise FloatingPointError("Invalid regression coefficient draw")
        kept = iteration > sampler["warmup"] and (iteration - sampler["warmup"]) % sampler["thin"] == 0
        if kept:
            total_components = -1
            if method == "mfm":
                support, probability = law.total_distribution(len(effects))
                total_components = int(rng.choice(support, p=probability))
            write_npz(directory / f"draw_{iteration:06d}.npz", labels=labels, beta=beta, effects=effects,
                      hyper=[ve, va, vb, alpha], total_components=total_components, iteration=iteration)
        if kept or iteration % sampler["checkpoint_every"] == 0 or iteration == total:
            memory = max(memory, peak_memory_bytes() or 0)
            write_npz(checkpoint, labels=labels, beta=beta, hyper=[ve, va, vb, alpha],
                      iterations=[iteration, accepted, rebuild_count], rng=json.dumps(rng.bit_generator.state),
                      elapsed=elapsed + perf_counter() - started, peak_memory=memory)
    return directory
