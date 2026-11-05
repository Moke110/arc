"""Known metric outcomes at the accepted evaluator seam (ADR 0042)."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from metrics import score_group, directional_curve

class MetricSanity(unittest.TestCase):
    def test_constant_scores_give_prevalence_without_tie_break_gains(self):
        ap, recall, curve = directional_curve(np.array([1, 0, 0, 0], bool), np.zeros(4))
        self.assertEqual(ap, .25)
        self.assertEqual(recall, 0)
        np.testing.assert_array_equal(curve, np.full(101, .25))

    def test_perfect_reversed_tied_and_undefined(self):
        truth = np.array([2., 1., -1., -2.])
        values, reasons, curves = score_group(truth, truth)
        self.assertEqual(values['AUPRC_positive@100'], 1.)
        self.assertEqual(values['Recall_negative@P80'], 1.)
        self.assertAlmostEqual(values['R2'], 1.)
        reversed_values, _, _ = score_group(truth, -truth)
        self.assertAlmostEqual(reversed_values['R2'], -3.)
        self.assertEqual(reversed_values['Recall_positive@P80'], 0.)
        values, reasons, _ = score_group(np.zeros(4), np.ones(4))
        self.assertTrue(np.isnan(values['AUPRC_positive@100']))
        self.assertEqual(reasons['R2'], 'zero_variance_truth')
        ap, p80, _ = directional_curve([1,0,1,0], [3,3,1,0])
        self.assertAlmostEqual(ap, 7 / 12)
        self.assertEqual(p80, 0.)

    def test_ap_labels_are_not_top100_prediction_truncation(self):
        truth = np.zeros(200)
        truth[0] = 1
        prediction = np.arange(200.)
        values, _, _ = score_group(truth, prediction)
        self.assertAlmostEqual(values['AUPRC_positive@100'], 1 / 200)

if __name__ == '__main__': unittest.main()
