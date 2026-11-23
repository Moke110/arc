"""Approved #226 seams: shared splits, train-only Means, drug-level scoring."""
import sys
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tenv5k.protocol import make_splits, mean_prediction, score_prediction

class ProtocolTests(unittest.TestCase):
    def test_mean_counts_shared_original_control_once_and_ignores_test_labels(self):
        x=np.array([[2,4],[4,8],[100,200],[0,0],[6,6]],dtype=np.float32)
        obs=pd.DataFrame({'drug_parent_inchikey':['A','A','B','',''], 'cell_line':['hiPSC-CM']*5})
        controls={0:[3,4],1:[3,4],2:[3]}
        pred=mean_prediction(x,obs,controls,[0,1],[2],'Mean','cold-drug')
        np.testing.assert_array_equal(pred,[[3,4.5]])
        x[2]=9999
        np.testing.assert_array_equal(mean_prediction(x,obs,controls,[0,1],[2],'Mean','cold-drug'),pred)
        np.testing.assert_array_equal(mean_prediction(x,obs,controls,[0,1],[2],'Mean_cell','cold-drug'),pred)
        self.assertIsNone(mean_prediction(x,obs,controls,[0,1],[2],'Mean_drug','cold-drug'))

    def test_splits_cover_all_outer_wells_once_and_obey_scenario(self):
        obs=pd.DataFrame({'drug_parent_inchikey':np.repeat(list('ABCDEFGHIJ'),3)},index=['w'+str(i) for i in range(30)])
        splits=make_splits(obs)
        pd.testing.assert_frame_equal(splits,make_splits(obs))
        self.assertEqual(sorted(splits.fold.unique()), [0, 1, 2])
        historical = make_splits(obs, folds=5)
        self.assertNotEqual(set(splits.loc[(splits.scenario=='cold-drug') & (splits.fold==0) & (splits.partition=='test'), 'profile_id']), set(historical.loc[(historical.scenario=='cold-drug') & (historical.fold==0) & (historical.partition=='test'), 'profile_id']))
        for scenario in ['cold-drug','warm-start']:
            scene=splits[splits.scenario==scenario]
            self.assertEqual(sorted(scene[scene.partition=='test'].row),list(range(30)))
            for fold in range(3):
                parts=scene[scene.fold==fold]
                sets={r:set(parts[parts.partition==r].drug_parent) for r in ['train','validation','test']}
                if scenario=='cold-drug':
                    self.assertFalse(sets['train']&sets['test']);self.assertFalse(sets['validation']&sets['test']);self.assertFalse(sets['train']&sets['validation'])
                else:self.assertTrue((sets['test']|sets['validation'])<=sets['train'])
        obs.iloc[0,0]='singleton'
        with self.assertRaisesRegex(ValueError,'at least two'):make_splits(obs)

    def test_drug_aggregation_precedes_correlation_and_constants_stay_visible(self):
        obs=pd.DataFrame({'drug_parent_inchikey':['A','A','B','C']})
        truth=np.array([[1,0,0],[0,1,0],[2,2,2],[0,1,2]],dtype=float)
        pred=np.array([[0,1,0],[1,0,0],[1,1,1],[1,1,1]],dtype=float)
        wells,groups=score_prediction(obs,truth,pred,np.zeros_like(truth))
        self.assertAlmostEqual(groups.iloc[0].pcc,1)
        self.assertAlmostEqual(wells.iloc[0].pcc,-0.5)
        self.assertTrue(np.isnan(groups.iloc[1].pcc));self.assertEqual(groups.iloc[1].rmse,1)
        self.assertEqual(groups.iloc[2].pcc,0)
        pred[0,0]=np.nan
        with self.assertRaisesRegex(ValueError,'non-finite'):score_prediction(obs,truth,pred,np.zeros_like(truth))

    def test_prediction_identity_order_is_checked_before_scoring(self):
        obs = pd.DataFrame({'drug_parent_inchikey': ['A', 'B']}, index=['well-a', 'well-b'])
        x = np.array([[0., 1.], [1., 0.]])
        with self.assertRaisesRegex(ValueError, 'identity'):
            score_prediction(obs, x, x, np.zeros_like(x), prediction_ids=['well-b', 'well-a'])

class FrozenProtocolTests(unittest.TestCase):
    def test_inputs_reject_historical_five_fold_protocol_before_payload_access(self):
        from tenv5k.inputs import Inputs, sha256
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'config.json').write_text(json.dumps({'schema': 'tenv5k-protocol-v1', 'folds': 5}))
            (root / 'manifest.json').write_text(json.dumps({'schema': 'tenv5k-inputs-v1', 'outputs': {'config.json': sha256(root / 'config.json')}}))
            with self.assertRaisesRegex(ValueError, 'protocol'):
                Inputs(root)

    def test_new_three_fold_freeze_runs_fifteen_means_and_three_na(self):
        from tenv5k.inputs import Inputs, sha256
        from tenv5k.protocol import PROTOCOL_SCHEMA, INPUT_SCHEMA
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        from means import run
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'inputs'; root.mkdir()
            obs = pd.DataFrame({'profile_id':['w'+str(i) for i in range(30)],
                'drug_parent_inchikey':np.repeat(list('ABCDEFGHIJ'),3), 'cell_line':['hiPSC-CM']*30})
            obs.to_csv(root/'population.csv',index=False)
            obs.to_csv(root/'all_obs.csv',index=False)
            pd.DataFrame({'drug_parent_inchikey':list('ABCDEFGHIJ')}).to_csv(root/'parents.csv',index=False)
            pd.DataFrame({'gene_index':[0,1,2]}).to_csv(root/'genes.csv',index=False)
            splits=make_splits(obs.set_index('profile_id',drop=False))
            splits.to_csv(root/'splits.csv',index=False)
            (root/'matched_controls.json').write_text(json.dumps({str(i):[] for i in range(30)}))
            x=np.arange(90,dtype=np.float32).reshape(30,3)
            for name,value in [('expression',x),('control',np.zeros_like(x)),('delta',x),('morgan',np.zeros((10,2048),dtype=np.uint8))]:
                np.save(root/(name+'.npy'),value)
            (root/'config.json').write_text(json.dumps({'schema':PROTOCOL_SCHEMA,'folds':3,'split_seed':226}))
            def seal():
                (root/'manifest.json').write_text(json.dumps({'schema':INPUT_SCHEMA,'outputs':{p.name:sha256(p) for p in root.iterdir() if p.name!='manifest.json'}}))
            seal()
            statuses=run(root,Path(directory)/'means')
            self.assertEqual(len(statuses),18)
            self.assertEqual(sum(s['status']=='complete' for s in statuses),15)
            self.assertEqual(sum(s['status']=='not-applicable' for s in statuses),3)
            splits=splits.drop(splits[(splits.scenario=='cold-drug') & (splits.partition=='test')].index[0])
            splits.to_csv(root/'splits.csv',index=False);seal()
            with self.assertRaisesRegex(ValueError,'OOF'):
                Inputs(root)
