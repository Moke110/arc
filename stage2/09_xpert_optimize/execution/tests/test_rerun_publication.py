"""Frozen rerun publication boundary and canonical native-axis scoring."""
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from rerun_publication import score_run

class PublicationTest(unittest.TestCase):
    def test_native_axis_and_full_truth_denominator(self):
        cases = pd.DataFrame({'case_id':['a','b'], 'group_id':['g','g']})
        result = score_run(cases, np.array([[2.,-1.],[2.,-1.]]),
                           np.array([[1.,-1.],[1.,-1.]]), [0,1],
                           lambda group, ids: np.array([1,-1,1,0]), 3.)
        self.assertAlmostEqual(result['metrics']['AUPRC'], 2/3)
        self.assertAlmostEqual(result['metrics']['AUPRC_positive'], .5)
        self.assertEqual(result['metrics']['AUPRC_negative'], 1.)
        self.assertAlmostEqual(result['metrics']['rmse'], np.sqrt(.5))
        self.assertEqual(result['metrics']['test_loss'], 3.)
        self.assertEqual(result['coverage']['covered_genes'], 2)

    def test_four_run_publication_replays_and_never_extrapolates(self):
        import json
        from tempfile import TemporaryDirectory
        from rerun_publication import MODELS, METRICS, publish, replay, ENTRIES
        with TemporaryDirectory() as directory:
            root = Path(directory)
            bundles = []
            for model in MODELS:
                for fold in (0,1):
                    for seed in (17,42):
                        checkpoint = root / f'{model}-{fold}-{seed}.pt'
                        checkpoint.write_bytes(f'{model}/{fold}/{seed}'.encode())
                        end = .25 if fold == 0 and seed == 17 else .75
                        points = [{'recall':0.,'precision':1.}, {'recall':end,'precision':.5}]
                        history = [{'update':10,'exposure':160,'train_terms':[1,1,1,1], 'validation_loss':.8},
                                   {'update':20,'exposure':320,'train_terms':[1,1,1,1], 'validation_loss':.9}]
                        bundles.append({'model':model,'fold':fold,'seed':seed,'checkpoint':str(checkpoint),
                                        'status':'complete','metrics':dict.fromkeys(METRICS,.5),
                                        'curves':{'AUPRC_positive':points,'AUPRC_negative':points},
                                        'phases':{p:{'history':history,'best_update':10} for p in ('tahoe','tenv')}})
            output = root/'publication'
            publish(bundles, output)
            self.assertEqual(set(p.name for p in output.iterdir()), ENTRIES)
            self.assertEqual(len(list((output/'checkpoints').iterdir())),24)
            data = json.loads((output/'figure_data/source.json').read_text())
            self.assertEqual(data['pr_summary'][MODELS[0]]['AUPRC_positive'][-1]['recall'],.25)
            self.assertTrue(replay(output))
            with self.assertRaises(ValueError):
                publish(bundles[:-1], root/'incomplete')
            self.assertFalse((root/'incomplete').exists())

    def test_missing_truth_is_pending_and_zero_significance_is_undefined(self):
        cases=pd.DataFrame({'case_id':['a'],'group_id':['g']})
        def missing(group,ids):
            raise KeyError('missing significance contrast')
        with self.assertRaises(KeyError):
            score_run(cases,np.array([[0.,1.]]),np.array([[0.,1.]]),[0,1],missing,1.)
        result=score_run(cases,np.array([[0.,1.]]),np.array([[0.,1.]]),[0,1],lambda g,ids:np.zeros(4),1.)
        self.assertIsNone(result['metrics']['AUPRC'])
        self.assertEqual(result['counts']['AUPRC'],{'defined':0,'undefined':1})
        self.assertEqual(result['groups'][0]['AUPRC_status'],'zero_significant')

    def test_ties_zero_predictions_and_constant_native_axis(self):
        cases=pd.DataFrame({'case_id':['a'],'group_id':['g']})
        result=score_run(cases,np.array([[1.,1.,0.]]),np.array([[2.,2.,2.]]),[0,1,2],
                         lambda g,ids:np.array([1,-1,1,0]),1.)
        self.assertAlmostEqual(result['metrics']['AUPRC'],1/6)
        self.assertIsNone(result['metrics']['pcc'])
        self.assertIsNone(result['metrics']['r2'])
        self.assertEqual(result['curves']['AUPRC'][-1]['recall'],1/3)

    def test_pr_summary_uses_only_the_shared_measured_recall_domain(self):
        from rerun_publication import summarize_pr, pr_summary_status
        curves = [[{'recall': start, 'precision': .8}, {'recall': .75, 'precision': .5}]
                  for start in (.1, .2, .15, .25)]
        summary = summarize_pr(curves)
        self.assertEqual(summary[0]['recall'], .25)
        self.assertEqual(summary[-1]['recall'], .75)
        self.assertEqual(pr_summary_status(curves), 'defined')
        curves[0] = [{'recall': .1, 'precision': .8}, {'recall': .2, 'precision': .5}]
        self.assertEqual(summarize_pr(curves), [])
        self.assertEqual(pr_summary_status(curves), 'no_common_recall_domain')
        curves[0] = []
        self.assertEqual(pr_summary_status(curves), 'undefined_run')
