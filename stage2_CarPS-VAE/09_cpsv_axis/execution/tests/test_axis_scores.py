"""Behavior of the held-out, fixed-denominator disease comparison."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axis_scores import COVERAGE, score_deltas


class HeldOutScores(unittest.TestCase):
    def test_direction_coverage_counts_only_equal_nonzero_signs_on_axis(self):
        genes = np.array([9, 6, 1, 19])
        measured = np.array([[1., -2., 0., 7.], [0., 0., 0., 0.]])
        predicted = np.array([[2., 1., -4., 9.], [0., 0., 0., 0.]])
        conditions, scores = score_deltas(measured, predicted, genes, {1, 6, 9, 21}, 3481,
                                           {"axis": 0.3, "disease": 0.7})
        disease = conditions[conditions.scope.eq("disease")]
        self.assertEqual(disease[COVERAGE].tolist(), [100 / 3481, 0.0])
        self.assertTrue(math.isnan(disease["Δ PCC"].iloc[1]))
        coverage = scores[scores.metric.eq(COVERAGE)].iloc[0]
        self.assertEqual(coverage.value, 50 / 3481)
        self.assertEqual(coverage.overlap, 3)
        self.assertEqual(coverage.n_eligible, 2)
        full_mmd = scores[scores.scope.eq("axis") & scores.metric.eq("Δ MMD")].value.iloc[0]
        disease_mmd = scores[scores.scope.eq("disease") & scores.metric.eq("Δ MMD")].value.iloc[0]
        self.assertNotEqual(full_mmd, disease_mmd)

    def test_coverage_is_paired_by_fold_and_seed_and_averaged_equally(self):
        from axis_scores import paired_scores
        roster = [("hvg-2k", 0, 17), ("hvg-2k", 1, 17),
                  ("denoise-5k", 0, 17), ("denoise-5k", 1, 17)]
        rows = [
            dict(axis=axis, fold=fold, seed=seed, scope="disease", metric=COVERAGE,
                 value=value, n_conditions=count, n_eligible=count, n_undefined=0, overlap=3)
            for (axis, fold, seed), value, count in zip(roster, (25., 16., 17., 18.), (1, 5, 1, 5), strict=True)
        ]
        paired, summary = paired_scores(pd.DataFrame(rows).iloc[::-1], roster, metrics={"disease": (COVERAGE,)})
        self.assertEqual(paired.loc[paired.axis.eq("hvg-2k"), "paired_difference"].tolist(), [8., -2.])
        self.assertEqual(summary.iloc[0].paired_mean, 3.)
        self.assertAlmostEqual(summary.iloc[0].paired_sd, math.sqrt(50))
