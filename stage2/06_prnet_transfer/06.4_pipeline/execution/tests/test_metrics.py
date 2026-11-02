"""Metric correctness: the frozen direction-AUPRC definition and its vectorised twin.

The vectorised AP matrix must agree with ``sklearn.metrics.average_precision_score``
— including exact ties — because checkpoint selection depends on it every epoch.
"""

from __future__ import annotations

import unittest

import numpy as np
from sklearn.metrics import average_precision_score

import _support  # noqa: F401  (path setup)
from ptlib import metrics


class TrueSetTests(unittest.TestCase):
    def test_top_bottom_are_disjoint_and_signed(self):
        rng = np.random.default_rng(0)
        delta = rng.normal(size=500)
        up, down = metrics.true_set_indices(delta, 200)
        self.assertEqual(up.size, 200)
        self.assertEqual(down.size, 200)
        self.assertTrue((delta[up] > 0).all())
        self.assertTrue((delta[down] < 0).all())
        self.assertEqual(np.intersect1d(up, down).size, 0)
        # top-k by delta value
        self.assertEqual(up[0], int(np.argmax(delta)))

    def test_fallback_when_fewer_nonzero_genes_than_k(self):
        delta = np.array([0.0, 1.0, 0.0, -2.0, 0.0, 0.5])
        up, down = metrics.true_set_indices(delta, 5)
        self.assertEqual(sorted(up.tolist()), [1, 5])
        self.assertEqual(sorted(down.tolist()), [3])

    def test_matrix_masks_match_scalar_path(self):
        rng = np.random.default_rng(1)
        delta = rng.normal(size=(7, 300))
        delta[0, :] = 0.0  # degenerate condition
        delta[1, :] = 0.0
        delta[1, :3] = [1.0, -1.0, 0.0]
        up, down = metrics.true_set_masks(delta, 20)
        for row in range(delta.shape[0]):
            up_ref, down_ref = metrics.true_set_indices(delta[row], 20)
            np.testing.assert_array_equal(np.sort(np.flatnonzero(up[row])), np.sort(up_ref))
            np.testing.assert_array_equal(np.sort(np.flatnonzero(down[row])), np.sort(down_ref))


class AveragePrecisionTests(unittest.TestCase):
    def test_matches_sklearn_random_scores(self):
        rng = np.random.default_rng(2)
        scores = rng.normal(size=(9, 1000))
        labels = np.zeros((9, 1000), dtype=np.int8)
        for row in range(9):
            idx = rng.choice(1000, size=40, replace=False)
            labels[row, idx] = 1
        ours = metrics.average_precision_matrix(scores, labels)
        ref = np.array(
            [average_precision_score(labels[i], scores[i]) for i in range(9)]
        )
        np.testing.assert_allclose(ours, ref, rtol=1e-12, atol=1e-12)

    def test_matches_sklearn_with_exact_ties(self):
        rng = np.random.default_rng(3)
        scores = np.round(rng.normal(size=(6, 500)), 1)  # many ties
        labels = (rng.random((6, 500)) < 0.1).astype(np.int8)
        ours = metrics.average_precision_matrix(scores, labels)
        ref = np.array(
            [average_precision_score(labels[i], scores[i]) for i in range(6)]
        )
        np.testing.assert_allclose(ours, ref, rtol=1e-12, atol=1e-12)

    def test_constant_scores_give_positive_rate(self):
        scores = np.zeros((1, 400))
        labels = np.zeros((1, 400), dtype=np.int8)
        labels[0, :50] = 1
        ap = metrics.average_precision_matrix(scores, labels)
        self.assertAlmostEqual(float(ap[0]), 50 / 400, places=12)

    def test_no_positive_is_nan(self):
        scores = np.zeros((1, 10))
        labels = np.zeros((1, 10), dtype=np.int8)
        self.assertTrue(np.isnan(metrics.average_precision_matrix(scores, labels)[0]))


class DirectionAuprcTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(4)
        self.true = rng.normal(size=(5, 1000))
        self.true[0, :] *= 3  # one condition with a stronger signal

    def test_perfect_prediction_scores_one(self):
        out = metrics.direction_auprc_matrix(self.true, self.true.copy(), k=50)
        np.testing.assert_allclose(out["auprc"], 1.0, atol=1e-12)
        self.assertEqual(int(out["n_up"].min()), 50)

    def test_constant_prediction_is_near_base_rate(self):
        pred = np.zeros_like(self.true)
        out = metrics.direction_auprc_matrix(self.true, pred, k=50)
        # both directions are all-tied: AP = positives/genes
        np.testing.assert_allclose(out["auprc"], 50 / 1000, atol=1e-12)

    def test_sign_flipped_ranks_negatives_first(self):
        out = metrics.direction_auprc_matrix(self.true, -self.true, k=50)
        both = np.vstack([out["auprc_up"], out["auprc_down"]])
        self.assertTrue(np.all(both < 50 / 1000 + 1e-9))
        self.assertTrue(np.all(both >= 0.0))

    def test_random_prediction_below_perfect(self):
        rng = np.random.default_rng(5)
        noise = rng.normal(size=self.true.shape)
        out = metrics.direction_auprc_matrix(self.true, noise, k=50)
        self.assertLess(float(np.nanmean(out["auprc"])), 1.0)
        self.assertGreater(float(np.nanmean(out["auprc"])), 0.0)

    def test_replicate_aggregation_happens_before_ranking(self):
        true = self.true[:1]
        base = np.zeros((1, 1000))
        # single prediction vector; averaged replicates must reproduce it
        reps_pred = np.repeat(base, 3, axis=0)
        reps_true = np.repeat(true, 3, axis=0)
        cond_pred = metrics.aggregate_rows(reps_pred, np.zeros(3, dtype=int), 1)
        cond_true = metrics.aggregate_rows(reps_true, np.zeros(3, dtype=int), 1)
        np.testing.assert_allclose(cond_pred[0], base[0])
        np.testing.assert_allclose(cond_true[0], true[0])
        direct = metrics.direction_auprc_matrix(true, base, k=50)
        aggregated = metrics.direction_auprc_matrix(cond_true, cond_pred, k=50)
        np.testing.assert_allclose(direct["auprc"], aggregated["auprc"])

    def test_group_metrics_matches_matrix_version(self):
        for row in range(self.true.shape[0]):
            scalar = metrics.group_metrics(self.true[row], -self.true[row], k=50)
            matrix = metrics.direction_auprc_matrix(
                self.true[row : row + 1], -self.true[row : row + 1], k=50
            )
            self.assertAlmostEqual(scalar["auprc"], float(matrix["auprc"][0]), places=12)
            self.assertAlmostEqual(scalar["auprc_up"], float(matrix["auprc_up"][0]), places=12)
            self.assertAlmostEqual(
                scalar["auprc_down"], float(matrix["auprc_down"][0]), places=12
            )


class RegressionMetricTests(unittest.TestCase):
    def test_perfect_rmse_zero_and_correlation_one(self):
        rng = np.random.default_rng(6)
        true = rng.normal(size=(4, 200))
        out = metrics.regression_metrics_matrix(true, true.copy())
        np.testing.assert_allclose(out["rmse"], 0.0)
        np.testing.assert_allclose(out["pcc"], 1.0)
        np.testing.assert_allclose(out["spearman"], 1.0)

    def test_constant_prediction_nan_correlation_but_finite_rmse(self):
        rng = np.random.default_rng(7)
        true = rng.normal(size=(2, 100))
        pred = np.ones_like(true)
        out = metrics.regression_metrics_matrix(true, pred)
        self.assertTrue(np.all(np.isnan(out["pcc"])))
        self.assertTrue(np.all(np.isnan(out["spearman"])))
        self.assertTrue(np.all(np.isfinite(out["rmse"])))
        self.assertTrue(np.all(out["pred_constant"]))

    def test_constant_true(self):
        pred = np.zeros((1, 50))
        true = np.zeros((1, 50))
        out = metrics.regression_metrics_matrix(true, pred)
        self.assertTrue(np.isnan(out["pcc"][0]))
        self.assertTrue(out["true_constant"][0])
        self.assertAlmostEqual(float(out["rmse"][0]), 0.0)

    def test_macro_omits_nan(self):
        summary = metrics.macro([0.5, np.nan, 1.0])
        self.assertAlmostEqual(summary["mean"], 0.75)
        self.assertEqual(summary["nan"], 1)
        self.assertAlmostEqual(summary["min"], 0.5)

    def test_aggregate_rows_means_within_group(self):
        values = np.array([[1.0, 2.0], [3.0, 4.0], [10.0, 20.0]])
        groups = np.array([0, 0, 1])
        out = metrics.aggregate_rows(values, groups, 2)
        np.testing.assert_allclose(out[0], [2.0, 3.0])
        np.testing.assert_allclose(out[1], [10.0, 20.0])


if __name__ == "__main__":
    unittest.main()
