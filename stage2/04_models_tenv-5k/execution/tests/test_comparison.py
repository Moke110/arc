"""Approved OOF/parent/seed comparison seam; literal independent fixtures."""
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tenv5k.comparison import assemble_oof, seed_average, paired_comparisons, select_winner

class ComparisonTests(unittest.TestCase):
    def test_oof_reconstructs_population_order_and_rejects_duplicate(self):
        result = assemble_oof(['a','b','c'], [(['c'], np.array([[3.,4.]])), (['a','b'], np.array([[1.,2.],[2.,3.]]))])
        np.testing.assert_array_equal(result, [[1,2],[2,3],[3,4]])
        with self.assertRaises(ValueError):
            assemble_oof(['a','b'], [(['a'],np.ones((1,2))),(['a'],np.ones((1,2)))])
    def test_average_seed_metrics_before_parent_macro(self):
        rows = pd.DataFrame([dict(scenario='cold-drug',method='PRnet',seed=s,drug_parent=p,pcc=v,spearman=v,rmse=2.,r2=0.,constant_truth=False,constant_prediction=False,n_wells=1) for s,p,v in [(17,'a',0.),(42,'a',1.),(2026,'a',.5),(17,'b',1.),(42,'b',1.),(2026,'b',1.)]])
        parents = seed_average(rows)
        self.assertEqual(parents.pcc.mean(), .75)
        with self.assertRaises(ValueError): seed_average(rows.iloc[:-1])
    def test_pairing_constant_difference_exact_and_no_forced_winner(self):
        rows = pd.DataFrame([dict(scenario='cold-drug',method=m,drug_parent=p,pcc=v,spearman=v,rmse=2.,r2=0.,constant_truth=False,n_seeds=3) for m,v in [('PRnet',.4),('TranSiGen',.4),('Mean',.1),('Mean_cell',.1),('MLP_Morgan',.2)] for p in ['a','b','c']])
        pairs = paired_comparisons(rows, repetitions=10000, seed=2262026)
        pair=pairs[(pairs.candidate=='PRnet')&(pairs.reference=='Mean')&(pairs.metric=='pcc')].iloc[0]
        self.assertAlmostEqual(pair.low,.3); self.assertAlmostEqual(pair.high,.3)
        self.assertEqual(select_winner(pairs)['winner'], None)
    def test_nonfinite_is_not_a_pairwise_drop(self):
        rows=pd.DataFrame([dict(scenario='cold-drug',method=m,drug_parent='a',pcc=v,spearman=.2,rmse=1.,r2=.2,constant_truth=False,n_seeds=3) for m,v in [('PRnet',float('nan')),('TranSiGen',.1)]])
        with self.assertRaises(ValueError): paired_comparisons(rows)

if __name__=='__main__': unittest.main()
