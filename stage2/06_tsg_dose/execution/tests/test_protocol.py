"""Approved freeze → split → ordered scoring boundaries from #232."""
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'output/api'))
from tsg_dose.protocol import make_splits, score_prediction, parent_metrics, summary


def population():
    return pd.DataFrame([dict(profile_id=f'{p}-{d}-{r}',drug_parent_inchikey=p,dose_um=d)
        for p in ['a','b','c','d','e','f'] for d in [0.2,1.,10.] for r in range(2)]).set_index('profile_id',drop=False)


class ProtocolTests(unittest.TestCase):
    def test_parent_and_held_dose_partitions_preserve_groups(self):
        obs=population(); splits=make_splits(obs)
        for (scenario,fold), frame in splits.groupby(['scenario','fold']):
            sets={p:set(frame.loc[frame.partition==p,'row']) for p in ['train','validation','test']}
            self.assertFalse(sets['train']&sets['validation']);self.assertFalse(sets['train']&sets['test'])
            self.assertEqual(set.union(*sets.values()),set(range(len(obs))))
            if scenario in ['cold-drug','final']:
                parents={k:set(obs.iloc[sorted(v)].drug_parent_inchikey) for k,v in sets.items()}
                self.assertFalse(parents['train']&parents['validation']);self.assertFalse(parents['train']&parents['test'])
            if scenario=='held-dose':
                self.assertEqual(set(obs.iloc[sorted(sets['test'])].dose_um),{[0.2,1.,10.][fold]})
                self.assertTrue(set(obs.iloc[sorted(sets['test'])].drug_parent_inchikey)<=set(obs.iloc[sorted(sets['train'])].drug_parent_inchikey))
        relabelled=obs.copy();relabelled['treated_expression']=np.arange(len(obs))*1000
        pd.testing.assert_frame_equal(splits,make_splits(relabelled))
        cold=splits[(splits.scenario=='cold-drug')&(splits.partition=='test')]
        self.assertEqual(sorted(cold.row),list(range(len(obs))))

    def test_repeats_before_metrics_doses_then_parents_equal(self):
        obs=pd.DataFrame({'drug_parent_inchikey':['a','a','a','b'],'dose_um':[.2,.2,1.,.2]},index=['a1','a2','a3','b1'])
        truth=np.array([[0,2],[2,0],[0,2],[0,2]],dtype=float)
        pred=np.array([[0,2],[0,2],[2,0],[0,2]],dtype=float)
        wells,groups=score_prediction(obs,truth,pred,np.zeros_like(truth),prediction_ids=obs.index)
        # a/.2 truth is constant after repeat aggregation; a/1 is anticorrelated.
        self.assertEqual(len(groups),3)
        self.assertTrue(groups.iloc[0].constant_truth)
        self.assertAlmostEqual(groups[(groups.drug_parent=='a')&(groups.dose_um==1)].iloc[0].pcc,-1.)
        self.assertEqual(wells.dose_um.tolist(),[.2,.2,1.,.2])
        parents=parent_metrics(groups).set_index('drug_parent')
        self.assertEqual(parents.loc['a','undefined_pcc_count'],1)
        self.assertAlmostEqual(parents.loc['a','rmse'],1.5)
        self.assertAlmostEqual(summary(groups)['rmse'],.75)
        self.assertAlmostEqual(summary(groups)['pcc'],0.)
        with self.assertRaisesRegex(ValueError,'identity order'):
            score_prediction(obs,truth,pred,np.zeros_like(truth),prediction_ids=list(reversed(obs.index)))

if __name__=='__main__':unittest.main()
