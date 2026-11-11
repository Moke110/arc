import unittest
import numpy as np
from stage2.evaluation.metrics import METRICS, PooledPR, average_precision, pooled_pr, retrieval


class RetrievalTests(unittest.TestCase):
    def test_three_metrics_have_full_axis_denominators_and_direction(self):
        truth = np.array([1, -1, 1, -1, 0, 1])
        result = retrieval(np.array([4., -3., -2., 1., 0.]), truth, np.arange(5))
        self.assertEqual(tuple(result), METRICS)
        self.assertAlmostEqual(result['AUPRC'].ap, 2 / 5)
        self.assertAlmostEqual(result['AUPRC_positive'].ap, 1 / 3)
        self.assertAlmostEqual(result['AUPRC_negative'].ap, 1 / 2)
        self.assertAlmostEqual(result['AUPRC'].recall[-1], 2 / 5)
        self.assertNotAlmostEqual(result['AUPRC'].ap,
                                 (result['AUPRC_positive'].ap + result['AUPRC_negative'].ap) / 2)

    def test_wrong_sign_and_missing_coverage_stay_false_negative(self):
        result = retrieval(np.array([-8., 7.]), np.array([1, -1, 1]), [0, 1])
        self.assertEqual(result['AUPRC'].ap, 0.)
        self.assertEqual(result['AUPRC_positive'].ap, 0.)
        self.assertEqual(result['AUPRC_negative'].ap, 0.)

    def test_ties_enter_together(self):
        curve = average_precision(np.array([2., -2., 1.]), np.array([1, 0, 1]), 4)
        self.assertAlmostEqual(curve.ap, .5 / 4 + (2 / 3) / 4)
        self.assertEqual(curve.recall[-1], .5)

    def test_undefined_and_abstention(self):
        result = retrieval(np.zeros(2), np.array([1, 0, 0]), [0, 1])
        self.assertEqual(result['AUPRC'].ap, 0.)
        self.assertEqual(result['AUPRC_positive'].ap, 0.)
        self.assertIsNone(result['AUPRC_negative'])
        self.assertTrue(all(x is None for x in retrieval(np.ones(2), np.zeros(3), [0, 1]).values()))

    def test_fixed_population_and_real_endpoint(self):
        a = average_precision(np.array([2.]), np.array([1]), 10)
        b = average_precision(np.array([1.]), np.array([1]), 1)
        frame = pooled_pr([a, b])
        self.assertAlmostEqual(frame.recall.iloc[-1], 2 / 11)
        self.assertEqual(frame.precision.iloc[-1], 1.)
        self.assertTrue(frame.contributing_groups.eq(2).all())
        self.assertTrue(frame.recall.is_monotonic_increasing)
        self.assertEqual(pooled_pr([]).shape[0], 0)

    def test_streaming_preserves_global_maximum_and_zero_predictions(self):
        pool = PooledPR()
        pool.add(average_precision(np.array([3., 1.]), np.array([1, 0]), 2))
        pool.add(average_precision(np.array([3., 0.]), np.array([0, 1]), 3))
        pool.add(average_precision(np.zeros(2), np.ones(2), 2))
        frame = pool.frame()
        self.assertEqual(frame.precision.iloc[0], .5)
        self.assertEqual(frame.precision.iloc[-1], 1 / 3)
        self.assertEqual(frame.recall.iloc[-1], 1 / 7)
        self.assertTrue(frame.contributing_groups.eq(3).all())

    def test_invalid_coverage_and_truth_rejected(self):
        for columns in ([0, 0], [-1, 0], [0, 3]):
            with self.assertRaises(ValueError):
                retrieval(np.ones(2), np.array([1, 0, -1]), columns)
        with self.assertRaises(ValueError):
            average_precision(np.ones(2), np.array([1, 1]), 1)


if __name__ == '__main__':
    unittest.main()
