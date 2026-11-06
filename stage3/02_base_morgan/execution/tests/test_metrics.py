"""Known readout-direction examples at the standalone metric seam."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from metrics import evaluate_predictions


class MetricExamples(unittest.TestCase):
    def test_perfect_predictions_keep_absent_classes_in_macro_average(self):
        result = evaluate_predictions([-1, -1], [[1, 0, 0], [1, 0, 0]])
        self.assertAlmostEqual(result['macro_f1'], 1 / 3)
        self.assertEqual(result['precision'], [1, 0, 0])
        self.assertEqual(result['recall'], [1, 0, 0])
        self.assertEqual(result['confusion_matrix'], [[2, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertLess(result['log_loss'], 1e-14)

    def test_known_wrong_probabilities_have_fixed_confusion_and_loss(self):
        result = evaluate_predictions([-1, 0, 1], [[0.1, 0.8, 0.1], [0.1, 0.1, 0.8], [0.8, 0.1, 0.1]])
        self.assertEqual(result['macro_f1'], 0)
        self.assertEqual(result['confusion_matrix'], [[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        self.assertAlmostEqual(result['log_loss'], 2.302585092994046)

    def test_extreme_wrong_prediction_has_finite_frozen_float64_loss(self):
        self.assertAlmostEqual(evaluate_predictions([-1], [[0, 1, 0]])['log_loss'], 36.04365338911715)
        self.assertAlmostEqual(evaluate_predictions([-1, 0, 1], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])['macro_f1'], 1)

    def test_rejects_invalid_probability_and_truth_contracts(self):
        for probabilities in ([[float('nan'), 0, 1]], [[float('inf'), 0, 0]], [[-0.1, 0.1, 1]], [[0.2, 0.2, 0.2]], [[0.5, 0.5]], []):
            with self.subTest(probabilities=probabilities), self.assertRaises(ValueError):
                evaluate_predictions([-1], probabilities)
        for truth in ([2], [0.1], [], [[-1]]):
            with self.subTest(truth=truth), self.assertRaises(ValueError):
                evaluate_predictions(truth, [[1, 0, 0]])

    def test_argmax_consistency_includes_first_class_tie_break(self):
        self.assertEqual(evaluate_predictions([-1], [[0.5, 0.5, 0]], [-1])['precision'], [1, 0, 0])
        with self.assertRaisesRegex(ValueError, 'argmax'):
            evaluate_predictions([-1], [[0.5, 0.5, 0]], [0])
        with self.assertRaisesRegex(ValueError, 'argmax'):
            evaluate_predictions([-1], [[1, 0, 0]], [-1, -1])
