import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from stage2.evaluation.evaluate import bound, score
from stage2.evaluation.inventory import REPO, binding, sha
from stage2.evaluation.metrics import METRICS


def fixture(root):
    pd.DataFrame({'case_id': ['a', 'b'], 'group_id': ['g', 'g'], 'dose_um': [1., 1.]}).to_csv(root / 'cases.csv', index=False)
    pd.DataFrame({'hgnc_id': [1, 2], 'symbol': ['A', 'B']}).to_csv(root / 'genes.csv', index=False)
    np.save(root / 'prediction.npy', np.array([[2., -1.], [2., -1.]], dtype=np.float32))
    (root / 'manifest.json').write_text('{}')
    class FakeLabels:
        def __init__(self, dataset):
            self.root = root
            self.manifest_digest = sha(root / 'manifest.json')
            self.genes = pd.DataFrame({'hgnc_id': [1, 2, 3], 'symbol': ['A', 'B', 'C']})
        def verify(self):
            self.assert_digest = sha(root / 'manifest.json')
            if self.assert_digest != self.manifest_digest:
                raise ValueError('label manifest changed')
        def truth(self, group, members):
            assert group == 'g' and members == ['a', 'b']
            return np.array([1, -1, 1])
    item = {'step': 9, 'dataset': 'tenv', 'regime': 'cold-drug', 'role': 'test',
            'model': 'covered', 'fold': 0, 'seed': 17, 'cases': binding(root / 'cases.csv'),
            'genes': binding(root / 'genes.csv'), 'prediction': binding(root / 'prediction.npy')}
    return item, FakeLabels


class ScoringTest(unittest.TestCase):
    def test_signed_coverage_artifacts_and_deterministic_replay(self):
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary)
            item, labels = fixture(root)
            with patch('stage2.evaluation.evaluate.Labels', labels):
                score([item], root / 'first')
                score([item], root / 'second')
            first = json.loads((root / 'first/manifest.json').read_text())
            second = json.loads((root / 'second/manifest.json').read_text())
            self.assertEqual(first, second)
            rows = pd.read_csv(root / 'first/per_run.csv')
            self.assertAlmostEqual(rows.AUPRC.iloc[0], 2 / 3)
            self.assertAlmostEqual(rows.AUPRC_positive.iloc[0], 1 / 2)
            self.assertEqual(rows.AUPRC_negative.iloc[0], 1.)
            self.assertEqual(set(pd.read_csv(root / 'first/pr.csv').metric), set(METRICS))
            self.assertFalse(any('@' in column or 'ap100' in column for column in rows.columns))
            for name, digest in first['files'].items():
                self.assertEqual(sha(root / 'second' / name), digest)

    def test_changed_prediction_is_rejected_even_after_cached_verification(self):
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary)
            item, labels = fixture(root)
            bound(item['prediction'])
            np.save(root / 'prediction.npy', np.array([[3., -1.], [3., -1.]], dtype=np.float32))
            with patch('stage2.evaluation.evaluate.Labels', labels), self.assertRaisesRegex(ValueError, 'changed'):
                score([item], root / 'result')
            self.assertFalse((root / 'result/manifest.json').exists())

    def test_rehashed_prediction_with_wrong_case_order_or_shape_is_rejected(self):
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary)
            item, labels = fixture(root)
            (root / 'order.json').write_text(json.dumps([['b', 'g'], ['a', 'g']]))
            item['prediction_order'] = binding(root / 'order.json')
            with patch('stage2.evaluation.evaluate.Labels', labels), self.assertRaisesRegex(ValueError, 'order'):
                score([item], root / 'wrong_order')
            del item['prediction_order']
            np.save(root / 'prediction.npy', np.ones((3, 2), dtype='float32'))
            item['prediction'] = binding(root / 'prediction.npy')
            with patch('stage2.evaluation.evaluate.Labels', labels), self.assertRaisesRegex(ValueError, 'shape'):
                score([item], root / 'wrong_shape')

    def test_reference_seed_copies_keep_outcomes_but_not_extra_variance_samples(self):
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary); template, labels = fixture(root)
            entries = []
            for fold, vector in enumerate(([2., -1.], [2., 1.], [-2., 1.])):
                path = root / f'prediction-{fold}.npy'
                np.save(path, np.array([vector, vector], dtype='float32'))
                for seed in (17, 42, 2026):
                    entries.append({**template, 'model': 'delta_mean', 'fold': fold, 'seed': seed,
                                    'prediction': binding(path)})
            with patch('stage2.evaluation.evaluate.Labels', labels):
                score(entries, root / 'result')
            self.assertEqual(len(pd.read_csv(root / 'result/per_run.csv')), 9)
            summary = pd.read_csv(root / 'result/summary.csv').iloc[0]
            self.assertEqual(summary.AUPRC_count, 3)
            self.assertAlmostEqual(summary.AUPRC_std, 1 / 3)

    def test_rehashed_prediction_cannot_change_its_selected_checkpoint(self):
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary); item, labels = fixture(root)
            selected = {'checkpoint_sha256': 'a' * 64}
            (root / 'run.json').write_text(json.dumps(selected))
            item['selected_run'] = binding(root / 'run.json')
            manifest = {**selected, 'prediction_sha256': item['prediction']['sha256'],
                        'run_sha256': item['selected_run']['sha256']}
            manifest['checkpoint_sha256'] = 'b' * 64
            (root / 'prediction.json').write_text(json.dumps(manifest))
            item['prediction_manifest'] = binding(root / 'prediction.json')
            with patch('stage2.evaluation.evaluate.Labels', labels), self.assertRaisesRegex(ValueError, 'selected checkpoint'):
                score([item], root / 'result')
            self.assertFalse((root / 'result/manifest.json').exists())

    def test_completed_release_rejects_changed_physical_pr_points(self):
        from stage2.evaluation.publish import verify_existing
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary); item, labels = fixture(root)
            with patch('stage2.evaluation.evaluate.Labels', labels):
                score([item], root / 'result')
            verify_existing(root / 'result')
            path = root / 'result/pr.csv'
            path.write_text(path.read_text() + '\n')
            with self.assertRaisesRegex(ValueError, 'published artifact changed'):
                verify_existing(root / 'result')

    def test_per_group_and_dose_continuous_metrics_keep_undefined_support(self):
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary); item, labels = fixture(root)
            cases = pd.read_csv(root / 'cases.csv')
            cases['dose_um'] = .2
            cases = pd.concat([cases, pd.DataFrame({'case_id':['c'], 'group_id':['h'], 'dose_um':[1.]})], ignore_index=True)
            cases.to_csv(root / 'cases.csv', index=False)
            np.save(root / 'prediction.npy', np.array([[2.,-1.],[2.,-1.],[0.,0.]], dtype='float32'))
            axis = pd.DataFrame({'hgnc_id': range(1,18165), 'symbol': [f'G{i}' for i in range(18164)]})
            axis.to_csv(root / 'axis.csv', index=False)
            pd.DataFrame({'group_id':['g','h']}).to_csv(root / 'groups.csv', index=False)
            truth = np.zeros((2,18164), dtype='float32'); truth[0,:2]=[2.,-1.];truth[1,:2]=[1.,3.]
            np.save(root / 'truth.npy', truth)
            item.update(cases=binding(root/'cases.csv'), prediction=binding(root/'prediction.npy'),
                        continuous_axis=binding(root/'axis.csv'),
                        continuous_truth=[{'groups':binding(root/'groups.csv'),'delta':binding(root/'truth.npy')}])
            with patch('stage2.evaluation.evaluate.Labels', labels), patch.object(labels, 'truth', lambda self,g,ids: np.array([1,-1,1]) if g=='g' else np.zeros(3)):
                score([item], root / 'result')
            groups = pd.read_csv(root / 'result/per_group.csv').set_index('group_id')
            self.assertEqual(groups.loc['g','PCC'], 1.)
            self.assertTrue(np.isnan(groups.loc['h','PCC']))
            doses = pd.read_csv(root / 'result/per_dose.csv').set_index('dose_um')
            self.assertEqual(doses.loc[.2,'PCC_defined_groups'], 1)
            self.assertEqual(doses.loc[1.,'PCC_undefined_groups'], 1)
