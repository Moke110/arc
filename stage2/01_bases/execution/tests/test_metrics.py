import unittest
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.metrics import METRICS, PooledPR, average_precision, pooled_pr, retrieval


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

class GroupedRetrievalTests(unittest.TestCase):
    def test_mean_well_delta_then_rank_and_equal_group_run_weight(self):
        import pandas as pd
        from bases.metrics import score_groups, summarize_runs
        class Truth:
            def truth(self, parent, dose, sample_ids):
                assert (parent, dose, set(sample_ids)) in [('A', 1., {'a','b'}), ('B', 1., {'c'})]
                return np.array([1,-1,1,-1,0,1]) if parent == 'A' else np.zeros(6)
        samples = pd.DataFrame({'sample_id':['a','b','c'], 'drug_parent':['A','A','B'], 'dose_um':[1.,1.,1.]})
        prediction = np.array([[8.,-6.,-4.,2.,0.], [0.,0.,0.,0.,0.], [1.,1.,1.,1.,1.]])
        result = score_groups(prediction, samples, Truth(), np.arange(5))
        self.assertAlmostEqual(result['summary']['AUPRC'], .4)
        self.assertEqual(result['summary']['eligible_groups']['AUPRC'], 1)
        self.assertEqual(result['summary']['groups'], 2)
        summary = summarize_runs([result['summary'], {'AUPRC':.2, 'AUPRC_positive':None, 'AUPRC_negative':.5}])
        self.assertAlmostEqual(summary['AUPRC'], .3)
        self.assertAlmostEqual(summary['AUPRC_positive'], 1/3)
        self.assertTrue(result['curves']['AUPRC'].eligible_groups.eq(1).all())

    def test_continuous_metrics_keep_native_numeric_definition(self):
        from bases.metrics import continuous
        result = continuous([1,2,3], [2,4,6])
        self.assertAlmostEqual(result['PCC'], 1.)
        self.assertAlmostEqual(result['Spearman'], 1.)
        self.assertAlmostEqual(result['RMSE'], (14/3)**.5)
        self.assertAlmostEqual(result['R2'], -6.)
        self.assertIsNone(continuous([1,1], [2,2])['PCC'])

    def test_fractional_indices_and_tahoe_retrieval_are_rejected(self):
        from bases.metrics import score_groups
        with self.assertRaises(ValueError):
            retrieval(np.ones(2),np.array([1,0,-1]),[0.5,1])
        with self.assertRaisesRegex(ValueError,'tenv'):
            score_groups(None,None,None,None,source='tahoe')


if __name__ == '__main__':
    unittest.main()
