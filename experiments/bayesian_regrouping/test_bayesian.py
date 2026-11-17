"""Offline numerical and sampler acceptance tests for the Bayesian pilot."""

from dataclasses import replace
import inspect
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
from scipy.stats import multivariate_normal, norm

import latent_lmm
from evaluation import Assignment, Draw, conditional_imspe, load_draws, representative
from latent_lmm import fit_chain, gaussian_conditional, update_gaussian
from partition import PartitionLaw, covariate_log_density
from run_experiment import DEFAULT_RGLL, HERE, Job, RGLL, evaluate, make_evaluator


def pilot():
    return json.loads((HERE / "config.pilot.json").read_text(encoding="utf-8"))


class PartitionTests(unittest.TestCase):
    def test_crp_existing_and_new_mass(self):
        law = PartitionLaw("dpm", 5, {})
        np.testing.assert_allclose(np.exp(law.weights(np.array([2, 2]), 1.5)), [2, 2, 1.5])

    def test_mfm_predictive_matches_count_mixture(self):
        law = PartitionLaw("mfm", 6, dict(component_support=8, geometric_probability=0.3, gamma=1.7))
        support, posterior = law.total_distribution(2)
        assert posterior[0] == 0 and posterior[2] > 0
        conditional = 1.7 * (support - 2) / (6 + 1.7 * support)
        expected = conditional @ posterior
        weights = np.exp(law.weights(np.array([4, 2]), 1.0, predictive=True))
        np.testing.assert_allclose(weights[-1] / weights.sum(), expected, rtol=1e-12)
        at_limit = law.weights(np.ones(8), 1.0, predictive=True)
        self.assertEqual(at_limit[-1], -np.inf)

    def test_ppmx_uses_integrated_predictive_and_new_group(self):
        spec = dict(mean_precision=2.0, within_variance=1.5)
        log_density = covariate_log_density([1.1], np.array([3]), np.array([[3.0]]), spec)
        expected = norm.logpdf(1.1, loc=[3 / 5, 0], scale=np.sqrt(1.5 * np.array([1 + 1 / 5, 1 + 1 / 2])))
        np.testing.assert_allclose(log_density, expected)
        train = np.array([[-2.1], [-2.0], [1.9], [2.0]])
        assign = Assignment(train, "ppmx", pilot())
        draw = Draw(np.array([0, 0, 1, 1]), np.zeros(2), np.array([[-1.0], [1.0]]), 1, 1, 0, 1, -1, 1, 1)
        w = assign.membership(np.array([[-2.0], [2.0]]), draw)
        self.assertGreater(w[0, 0], w[0, 1])
        self.assertGreater(w[1, 1], w[1, 0])
        self.assertTrue((w[:, -1] > 0).all())
        np.testing.assert_allclose(w.sum(axis=1), 1)


class LMMTests(unittest.TestCase):
    def test_shared_group_joint_likelihood(self):
        z = np.array([[1.0, -0.4], [1.0, 1.2]])
        y = np.array([0.1, 1.1])
        variances = np.array([1.5, 0.3])
        ve = 0.8
        full = (z * variances) @ z.T + ve * np.eye(len(y))
        mean, covariance = gaussian_conditional(z[:1], y[:1], variances, ve)
        second = norm.logpdf(y[1], z[1] @ mean, np.sqrt(ve + z[1] @ covariance @ z[1]))
        np.testing.assert_allclose(norm.logpdf(y[0], 0, np.sqrt(full[0, 0])) + second,
                                   multivariate_normal.logpdf(y, cov=full), rtol=1e-12)
        added = update_gaussian(mean, covariance, z[1], y[1], ve, 1)
        expected = gaussian_conditional(z, y, variances, ve)
        for a, b in zip(added, expected):
            np.testing.assert_allclose(a, b, rtol=1e-12)
        removed = update_gaussian(*added, z[1], y[1], ve, -1)
        np.testing.assert_allclose(removed[0], mean)
        np.testing.assert_allclose(removed[1], covariance)

    def test_eq10_matches_reused_rgll_evaluator(self):
        rng = np.random.default_rng(20)
        x, ref = rng.normal(size=(96, 2)), rng.normal(size=(64, 2))
        labels = np.repeat(np.arange(3), 32)
        target = np.arange(len(ref)) % 3
        for vb in (0, 0.1):
            job = Job(1, 10, 0.5, vb, 1)
            y = x[:, 0] + labels + (labels * x[:, 1] if vb else 0) + rng.normal(size=len(x))
            original, = make_evaluator(x, y, ref, job, pilot()["evaluation"], "cpu", solver="torch")([(labels, target)])
            ve, va, slope = original.variance[2], original.variance[0], original.variance[1]
            variance = np.array([va, slope, slope] if vb else [va])
            np.testing.assert_allclose(conditional_imspe(x, ref, labels, target, variance, ve), original.imspe,
                                       rtol=1e-6, atol=1e-8)

    def test_partition_permutation_and_prediction(self):
        x = np.array([[-2.0], [-1.0], [1.0], [2.0]])
        labels = np.array([0, 0, 1, 1])
        draw = Draw(labels, np.array([1., 0.5]), np.array([[-1.0], [1.0]]), 1, 1, 0, 1, -1, 1, 1)
        relabeled = replace(draw, labels=1 - labels, effects=draw.effects[::-1])
        first = representative([draw, draw], 60, 7)
        self.assertEqual(first, representative([draw, relabeled], 60, 7))
        for method in ("mfm", "dpm", "ppmx"):
            gate = Assignment(x, method, pilot())
            pred1 = gate.predict_new(x, [draw])
            pred2 = gate.predict_new(x, [relabeled])
            for a, b in zip(pred1[:3], pred2[:3]):
                np.testing.assert_allclose(a, b)
        score = conditional_imspe(x, x, labels, labels, draw.variances, draw.ve)
        np.testing.assert_allclose(score, conditional_imspe(x, x, 1 - labels, 1 - labels, draw.variances, draw.ve))

    def test_prediction_scoring_cannot_change_the_fitted_partition(self):
        rng = np.random.default_rng(42)
        x = rng.normal(size=(36, 1))
        labels = (x[:, 0] > 0).astype(int)
        y = labels + x[:, 0] + rng.normal(size=len(x))
        draw = Draw(labels.copy(), np.array([1., 0.5]), np.array([[-1.0], [1.0]]), 1, 1, 0, 1, -1, 1, 1)
        job = Job(1, 10, 0.5, 0, 1)
        with TemporaryDirectory() as first, TemporaryDirectory() as second:
            cfg = pilot()
            a = evaluate([draw], x, y, x, x, np.zeros(len(x)), labels, "dpm", cfg, job, Path(first), "cpu")
            b = evaluate([draw], x, y, x, x, np.repeat(20.0, len(x)), labels, "dpm", cfg, job, Path(second), "cpu")
        self.assertIsNotNone(a["plug_in_hard_IMSPE"])
        self.assertEqual(a["plug_in_hard_IMSPE"], b["plug_in_hard_IMSPE"])
        self.assertEqual(a["K_representative"], b["K_representative"])
        self.assertNotEqual(a["test_MSPE_posterior_predictive"], b["test_MSPE_posterior_predictive"])
        np.testing.assert_array_equal(draw.labels, labels)
        self.assertTrue({"R", "truth", "test_y"}.isdisjoint(inspect.signature(fit_chain).parameters))


class SamplerTests(unittest.TestCase):
    def test_draws_and_resume(self):
        rng = np.random.default_rng(15)
        x = rng.normal(size=(24, 2))
        y = rng.normal(size=len(x))
        cfg = pilot()
        cfg["sampler"].update(warmup=4, draws=16, thin=1, checkpoint_every=2)
        for method in ("mfm", "dpm", "ppmx"):
            with self.subTest(method=method), TemporaryDirectory() as temp:
                root = Path(temp)
                full, resumed = root / "full", root / "resumed"
                fit_chain(x, y, method, False, cfg, 319, full)
                write = latent_lmm.write_npz

                def interrupt(path, **values):
                    write(path, **values)
                    if path.name == "checkpoint.npz" and values["iterations"][0] == 6:
                        raise InterruptedError("test interrupt after completing checkpoint")

                with patch.object(latent_lmm, "write_npz", side_effect=interrupt):
                    with self.assertRaises(InterruptedError):
                        fit_chain(x, y, method, False, cfg, 319, resumed)
                fit_chain(x, y, method, False, cfg, 319, resumed)
                draws = load_draws([full])
                self.assertEqual(len(draws), 16)
                self.assertGreater(len({d.k for d in draws}), 1)
                for a, b in zip(draws, load_draws([resumed]), strict=True):
                    np.testing.assert_array_equal(a.labels, b.labels)
                    np.testing.assert_array_equal(a.effects, b.effects)
                    np.testing.assert_array_equal(a.variances, b.variances)
                    self.assertEqual(a.total, b.total)
                    if method == "mfm":
                        self.assertGreaterEqual(a.total, a.k)


@unittest.skipUnless(os.environ.get("RSCRIPT"), "Set RSCRIPT to opt into the full RGLL generator checks")
class GeneratorTests(unittest.TestCase):
    def test_paper_case_adapter(self):
        with RGLL(os.environ["RSCRIPT"], DEFAULT_RGLL) as rgll:
            for case in range(1, 8):
                x, reference, test, y, test_y, groups, _ = rgll.data(Job(case, 10, 0.5, 0, 1))
                self.assertEqual((x.shape, reference.shape, test.shape), ((7500, 50), (2500, 50), (2500, 50)))
                self.assertEqual((len(y), len(test_y), len(np.unique(groups))), (7500, 2500, 10))
                self.assertFalse(np.array_equal(reference, test))


if __name__ == "__main__":
    unittest.main(verbosity=2)
