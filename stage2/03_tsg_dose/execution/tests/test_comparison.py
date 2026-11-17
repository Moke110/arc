import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'output/api'))
from tsg_dose.comparison import choose_method,seed_average,paired_comparisons
from tsg_dose.protocol import METHODS,SEEDS

class Comparison(unittest.TestCase):
    def test_choice_uses_seed_metrics_and_frozen_ties(self):
        rows=[dict(method=m,seed=s,pcc=.2,parameter_count=20) for m in METHODS for s in SEEDS]
        rows[0]['parameter_count']=10
        for r in rows:
            if r['method']=='no-dose':r['parameter_count']=10
        self.assertEqual(choose_method(pd.DataFrame(rows)),'no-dose')
        for r in rows:
            if r['method']=='film':r['pcc']=.3
        self.assertEqual(choose_method(pd.DataFrame(rows)),'film')
        with self.assertRaises(ValueError):choose_method(pd.DataFrame(rows[:-1]))

    def test_seed_metrics_not_prediction_averages_and_parent_blocks(self):
        rows=[]
        for m in ['selected','no-dose']:
            for parent in ['a','b']:
                for s,value in zip(SEEDS,[0.,.3,.6]):
                    rows.append(dict(method=m,drug_parent=parent,seed=s,pcc=value+(0.25 if m=='selected' else 0.),spearman=value,rmse=1.,r2=value,n_doses=3,n_wells=9,truth_constant_count=0,prediction_constant_count=0,undefined_pcc_count=0,undefined_r2_count=0))
        averaged=seed_average(pd.DataFrame(rows))
        self.assertAlmostEqual(averaged[averaged.method=='no-dose'].pcc.iloc[0],.3)
        pair=paired_comparisons(averaged,repetitions=10000)
        main=pair[(pair.candidate=='selected')&(pair.metric=='pcc')].iloc[0]
        self.assertAlmostEqual(main.estimate,.25)
        self.assertAlmostEqual(main.low,.25)
        self.assertAlmostEqual(main.high,.25)
        self.assertEqual(main.n_parents,2)
        self.assertIn('fixed',main.inference)
        with self.assertRaises(ValueError):seed_average(pd.DataFrame(rows[:-1]))

if __name__=='__main__':unittest.main()
