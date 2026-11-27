"""Worked and boundary examples at the fixed-prediction scoring seam."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

STEP = (Path(__file__).resolve().parents[4] / "04_tsg_axis/execution/tests" / Path(__file__).name).parents[2]
sys.path[:0] = [str(STEP / "output/api"), str(STEP.parent / "util/output/api")]

from stage2_util import metrics
from tsg_axis.metrics import (distribution_scores, fold_bandwidth, signed_disease_f1, signed_jaccard,
                              summarize_conditions, summarize_folds)


class AxisMetricTest(unittest.TestCase):
    def test_top_ties_follow_frozen_index_and_zero_is_excluded(self):
        expected = signed_jaccard([1, 1, -1, 0, 0], [0, 1, 1, 0, 0],
                                   [13, 4, 8, 1, 17], limit=1)
        self.assertEqual(expected, (0.5, 1.0, 0.0))
        top, pos, neg = signed_jaccard([0, 0], [0, 1], [3, 2])
        self.assertEqual((top, pos), (0.0, 0.0))
        self.assertTrue(math.isnan(neg))
        self.assertTrue(all(math.isnan(x) for x in signed_jaccard([0], [0], [1])))

    def test_weighted_signed_jaccard_on_a_2k_axis(self):
        genes = np.arange(2000)
        truth, prediction = np.zeros(2000), np.zeros(2000)
        truth[:1000] = np.arange(1000, 0, -1)
        prediction[:100] = truth[:100]
        prediction[1000:1900] = np.arange(900, 0, -1)
        top200, weighted = metrics.signed_jaccard_scores(truth, prediction, genes)
        self.assertEqual(top200, 1)
        for cutoff, expected in ((200, 1), (500, 1 / 4), (1000, 1 / 9), (2000, 1 / 19)):
            with self.subTest(cutoff=cutoff):
                self.assertEqual(metrics.signed_jaccard(truth, prediction, genes, limit=cutoff // 2)[0],
                                 expected)
        self.assertAlmostEqual(weighted, (0.125 + 0.25 / 4 + 0.5 / 9 + 1 / 19) / 1.875)

    def test_weighted_jaccard_keeps_signs_separate_and_empty_undefined(self):
        genes = [13, 4, 8, 1, 17]
        observed = [1, 1, -1, 0, 0]
        predicted = [0, 1, 1, 0, 0]
        top200, weighted = metrics.signed_jaccard_scores(observed, predicted, genes)
        self.assertEqual(top200, 1 / 4)
        self.assertAlmostEqual(weighted, 1 / 4)
        self.assertTrue(all(math.isnan(value) for value in
                            metrics.signed_jaccard_scores([0], [0], [0])))

    def test_off_axis_and_reversed_disease_hits(self):
        score = signed_disease_f1(np.array([1, -1, 1]), [1, 2, 4],
                                  np.array([0, 1, 1, -1, 0]), [1, 2, 3])
        self.assertEqual((score.tp, score.fp, score.fn), (1, 2, 2))
        self.assertEqual(score.f1, 1 / 3)
        cases = ([0], [0], [0])
        self.assertTrue(math.isnan(signed_disease_f1(*cases, [0]).f1))
        self.assertEqual(signed_disease_f1([1], [0], [0], [0]).f1, 0.0)
        self.assertEqual(signed_disease_f1([0], [0], [1], [0]).f1, 0.0)

    def test_distribution_includes_diagonal_terms_and_scales_by_width(self):
        truth = np.array([[0, 0], [2, 2]])
        pred = np.array([[0, 0], [0, 0]])
        energy, mmd = distribution_scores(truth, pred, 2.0)
        self.assertAlmostEqual(energy, 1.0)
        self.assertAlmostEqual(mmd, math.sqrt((1 - math.exp(-0.5)) / 2))
        self.assertEqual(distribution_scores(truth, truth[::-1], 2.0), (0.0, 0.0))
        train = np.array([[0, 0], [0, 0], [1, 1], [4, 4]])
        self.assertAlmostEqual(fold_bandwidth(train), 2.0)

    def test_aggregation_weights_conditions_then_folds_and_reports_na(self):
        folds = [summarize_conditions([0.1, 0.3, float("nan")]), summarize_conditions([0.8]),
                 summarize_conditions([float("nan")])]
        rows = pd.DataFrame([dict(axis="3k", fold=fold, metric="metric", **row)
                             for fold, row in enumerate(folds)])
        score = summarize_folds(rows, ["axis", "metric"]).iloc[0]
        self.assertEqual(score.value, 0.5)
        self.assertEqual((score.n_folds, score.n_eligible_folds, score.n_undefined_folds), (3, 2, 1))
        self.assertEqual((score.n_conditions, score.n_eligible, score.n_undefined), (5, 3, 2))
