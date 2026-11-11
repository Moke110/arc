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
            self.genes = pd.DataFrame({'hgnc_id': [1, 2, 3], 'symbol': ['A', 'B', 'C']})
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
