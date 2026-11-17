"""RGLL ML-EM and Eq. (10), with batched sufficient statistics."""

from dataclasses import dataclass

import numpy as np
import torch


EPS = 1e-8  # RGLL lmm_fast_aligned.cpp


@dataclass
class Fit:
    imspe: float
    beta: np.ndarray
    effects: np.ndarray
    variance: list[float]
    iterations: int
    converged: bool

    def predict(self, x: np.ndarray, labels: np.ndarray) -> np.ndarray:
        design = np.column_stack((np.ones(len(x)), x))
        random_design = design[:, : self.effects.shape[1]]
        return design @ self.beta + np.einsum("ij,ij->i", random_design, self.effects[labels])


class Evaluator:
    def __init__(self, x, y, ref, *, random_slope, device, tolerance, iterations, solver="torch"):
        opts = {"device": device, "dtype": torch.float64}
        self.train = torch.as_tensor(np.column_stack((np.ones(len(x)), x, y)), **opts)
        self.reference = torch.as_tensor(np.column_stack((np.ones(len(ref)), ref)), **opts)
        self.total = self.train.T @ self.train
        self.w = self.reference.T @ self.reference / len(ref)
        self.q = self.w.shape[0]
        self.rank = self.q if random_slope else 1
        self.eye = torch.eye(self.q, **opts)
        self.tolerance = tolerance
        self.max_iterations = iterations
        self.device = device
        self.solver = solver
        self.beta0, info = torch.linalg.solve_ex(self.total[:-1, :-1], self.total[:-1, -1])
        if int(info) or not torch.isfinite(self.beta0).all():
            raise ValueError("Fixed design is rank deficient")
        self.ve0 = ((self.train[:, -1] - self.train[:, :-1] @ self.beta0).square().mean()).clamp_min(EPS)
        if solver == "native":
            self.native_initial = (self.beta0.cpu().contiguous().numpy(), float(self.ve0))

    def _moments(self, data, groups, k):
        # BLAS forms moments on the accelerator. Grouping uses deterministic
        # sample order and hard leaf labels, including empty reference regions.
        row_ids = [np.flatnonzero(labels == j) for labels in groups for j in range(k)]
        lengths = np.array([len(ids) for ids in row_ids])
        index = np.zeros((len(row_ids), lengths.max()), dtype=np.int64)
        for row, ids in enumerate(row_ids):
            index[row, :len(ids)] = ids
        index = torch.as_tensor(index, device=self.device)
        valid = torch.arange(index.shape[1], device=self.device) < torch.as_tensor(lengths, device=self.device)[:, None]
        rows = data[index] * valid[..., None]
        return (rows.mT @ rows).reshape(len(groups), k, data.shape[-1], data.shape[-1])

    def __call__(self, assignments):
        # Same-K batches avoid artificial groups in the likelihood updates.
        k = int(assignments[0][0].max()) + 1
        gram = self._moments(self.train, [a[0] for a in assignments], k)
        w_k = self._moments(self.reference, [a[1] for a in assignments], k) / len(self.reference)
        if self.solver == "native":
            from native import fit_batch
            results = fit_batch(gram.cpu().contiguous().numpy(), w_k.cpu().contiguous().numpy(),
                                *self.native_initial, len(self.train), self.rank,
                                self.tolerance, self.max_iterations)
            fits = []
            for result in results:
                if result is None:
                    fits.append(None)
                else:
                    beta, effects, (score, a, slope, e, steps, done) = result
                    fits.append(Fit(float(score), beta, effects, [a, slope, e], int(steps), bool(done)))
            return fits
        return self.fit(gram, w_k)

    def fit(self, gram, w_k):
        # A=F'F, h=F'y. The random design uses either the intercept alone
        # or the full fixed design, so its moments are submatrices of A.
        b, k = gram.shape[:2]
        r, q = self.rank, self.q
        x = gram[..., :q, :q].contiguous()
        xy = gram[..., :q, q, None].contiguous()
        yy = gram[..., q, q, None, None]
        xtz, ztz, zty = x[..., :r], x[..., :r, :r], xy[..., :r, :]
        beta = self.beta0.expand(b, q)[..., None]
        ve = self.ve0.expand(b)
        va = (0.10 * ve).clamp_min(0.1)
        vb = (0.01 * ve).clamp_min(0.01) if r > 1 else torch.zeros_like(va)
        done = torch.zeros(b, dtype=torch.bool, device=self.device)
        used = torch.zeros(b, dtype=torch.int32, device=self.device)
        bad = torch.zeros_like(done)

        def posterior(ve, va, vb):
            if r == 1:
                return (ve[:, None, None, None] / va[:, None, None, None] + ztz).reciprocal()
            diag = torch.cat((ve[:, None] / va[:, None], (ve / vb)[:, None].expand(-1, r - 1)), 1)
            # solve_ex keeps CUDA kernels asynchronous. The right hand side is
            # an identity because both EM and IMSPE use the posterior diagonal.
            h, status = torch.linalg.solve_ex(
                ztz + torch.diag_embed(diag)[:, None], self.eye.expand(b, k, q, q)
            )
            nonlocal bad
            bad = bad | (status != 0).any(dim=1)
            return h

        def information(h):
            correction = xtz @ h
            m = (x - correction @ xtz.mT).sum(1)
            return (m + m.mT) * 0.5, correction

        h = posterior(ve, va, vb)
        for step in range(1, self.max_iterations + 1):
            br = beta[:, None]
            rg = zty - xtz.mT @ br
            u = h @ rg
            posterior_cov = ve[:, None, None, None] * h
            e_diag = posterior_cov.diagonal(dim1=-2, dim2=-1) + u[..., 0].square()
            va_new = e_diag[..., 0].mean(1).clamp_min(EPS)
            vb_new = e_diag[..., 1:].mean((1, 2)).clamp_min(EPS) if r > 1 else vb
            rss = yy - 2 * (br.mT @ xy) + (br.mT @ x @ br)
            rss = rss - 2 * (u.mT @ rg) + (u.mT @ ztz @ u)
            trace = (posterior_cov * ztz).sum((-1, -2))
            ve_new = ((rss[..., 0, 0] + trace).sum(1) / len(self.train)).clamp_min(EPS)

            # Freeze individually converged candidates while the batch continues.
            va_new = torch.where(done, va, va_new)
            vb_new = torch.where(done, vb, vb_new)
            ve_new = torch.where(done, ve, ve_new)
            h = posterior(ve_new, va_new, vb_new)
            m, correction = information(h)
            rhs = (xy - correction @ zty).sum(1)
            beta_new, status = torch.linalg.solve_ex(m, rhs)
            bad = bad | (status != 0)
            delta = torch.stack((
                (va_new - va).abs() / (1 + va.abs()),
                (vb_new - vb).abs() / (1 + vb.abs()),
                (ve_new - ve).abs() / (1 + ve.abs()),
                (beta_new - beta).norm(dim=1)[:, 0] / (1 + beta.norm(dim=1)[:, 0]),
            )).amax(0)
            beta = torch.where(done[:, None, None], beta, beta_new)
            ve, va, vb = ve_new, va_new, vb_new
            used += ~done
            done = done | (delta < self.tolerance)
            if step % 5 == 0 and bool((done | bad).all()):
                break

        # Schur complement form of Eq. (10): the fixed and cross terms
        # combine into E[(f - FZ H^-1 z)' C11 (f - FZ H^-1 z)].
        m, correction = information(h)
        c11, status = torch.linalg.solve_ex(m, self.eye.expand(b, q, q))
        bad = bad | (status != 0)
        c11 = c11 * ve[:, None, None]
        design_map = self.eye - torch.nn.functional.pad(correction, (0, q - r))
        fixed = (c11[:, None] * (design_map @ w_k @ design_map.mT)).sum((1, 2, 3))
        random = (h * w_k[..., :r, :r]).sum((1, 2, 3)) * ve
        imspe = ve + fixed + random
        u = (h @ (zty - xtz.mT @ beta[:, None]))[..., 0]
        beta = beta[..., 0]
        diagnostic = torch.stack((imspe, va, vb, ve, used, done, bad), dim=1).cpu().numpy()
        betas, effects = beta.cpu().numpy(), u.cpu().numpy()
        out = []
        for i, (score, a, slope, e, steps, converged, failed) in enumerate(diagnostic):
            if failed or not np.isfinite([score, a, slope, e]).all() or score <= 0:
                out.append(None)
            else:
                out.append(Fit(float(score), betas[i], effects[i], [a, slope, e], int(steps), bool(converged)))
        return out
