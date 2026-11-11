"""Step03 shares the canonical significance metric rather than ranking AP."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from metrics import METRICS, directional_curve, score_group, truth_labels


class MetricsTest(unittest.TestCase):
    def test_three_retrieval_metrics_and_continuous_metrics(self):
        self.assertEqual(METRICS, ('AUPRC', 'AUPRC_positive', 'AUPRC_negative', 'PCC', 'Spearman', 'RMSE', 'R2'))
        values, reasons, curves = score_group(np.array([2., -1., 1.]), np.array([2., -1., -1.]), np.array([1, -1, 1]))
        self.assertEqual(set(curves), set(METRICS[:3]))
        self.assertEqual(values['AUPRC_positive'], .5)
        self.assertEqual(values['AUPRC_negative'], .5)

    def test_rank_defined_scoring_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'retired'):
            directional_curve(np.ones(3), np.ones(3))
        with self.assertRaisesRegex(ValueError, 'DESeq2'):
            score_group(np.ones(3), np.ones(3), truth_labels(np.ones(3)))

    def test_frozen_input_rank_labels_still_reproduce(self):
        positive, negative = truth_labels(np.r_[np.arange(1, 110), -np.arange(1, 110)])
        self.assertEqual(positive.sum(), 100)
        self.assertEqual(negative.sum(), 100)
