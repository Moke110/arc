"""Behavior of the held-out full-axis comparison."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axis_scores import score_deltas


class HeldOutScores(unittest.TestCase):
    def test_full_axis_scores_signed_jaccard(self):
        genes = np.array([9, 6, 1, 19])
        measured = np.array([[1., -2., 0., 7.], [0., 0., 0., 0.]])
        predicted = np.array([[2., 1., -4., 9.], [0., 0., 0., 0.]])
        conditions, scores = score_deltas(measured, predicted, genes, 0.3)
        self.assertEqual(conditions["JAC_top200"].iloc[0], 2 / 5)
        self.assertTrue(math.isnan(conditions["JAC_top200"].iloc[1]))
        jaccard = scores.loc[scores.metric.eq("JAC_top200")].iloc[0]
        self.assertEqual(jaccard.n_eligible, 1)
        self.assertEqual(jaccard.n_undefined, 1)
        self.assertEqual(jaccard.width, 4)
        self.assertNotIn("scope", conditions.columns)

    def test_scores_are_paired_by_fold_and_seed_and_averaged_equally(self):
        from axis_scores import paired_scores
        roster = [("hvg-2k", 0, 17), ("hvg-2k", 1, 17),
                  ("denoise-5k", 0, 17), ("denoise-5k", 1, 17)]
        rows = [
            dict(axis=axis, fold=fold, seed=seed, metric="Δ PCC",
                 value=value, n_conditions=count, n_eligible=count, n_undefined=0, width=4)
            for (axis, fold, seed), value, count in zip(roster, (25., 16., 17., 18.), (1, 5, 1, 5), strict=True)
        ]
        paired, summary = paired_scores(pd.DataFrame(rows).iloc[::-1], roster, metrics=("Δ PCC",))
        self.assertEqual(paired.loc[paired.axis.eq("hvg-2k"), "paired_difference"].tolist(), [8., -2.])
        self.assertEqual(summary.iloc[0].paired_mean, 3.)
        self.assertAlmostEqual(summary.iloc[0].paired_sd, math.sqrt(50))
