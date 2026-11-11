"""CLI contract for fold-local rerun inputs and the 24-slot study."""
import json
from tempfile import TemporaryDirectory
from pathlib import Path
import unittest
import hashlib
import anndata as ad
import numpy as np
import pandas as pd
from test_landmark import make_fixture, run


class RerunInputsTest(unittest.TestCase):
    def test_two_folds_have_nested_axes_and_exact_frozen_membership(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            native, tenv = make_fixture(root)
            landmark = root / 'landmark'
            run('01_corpus/build.py', '--native', native, '--tenv-release', tenv,
                '--landmarks', root / 'landmarks.csv', '--tahoe-drugs', root / 'tahoe_drugs.parquet',
                '--output', landmark, '--fixture')
            for fold in (0, 1):
                output = root / f'fold_{fold}'
                run('01_corpus/gene_axes.py', '--native', native, '--landmark', landmark,
                    '--output', output, '--fixture-sizes', 1, 2, 3, '--rerun', '--fold', fold)
                previous = set()
                for name, count in [('G512', 1), ('G1024', 2), ('G2048', 3), ('G18164', 3)]:
                    predictor = output / 'views' / name / 'predictor'
                    manifest = json.loads((predictor / 'manifest.json').read_text())
                    self.assertEqual(manifest['ready_folds'], [fold])
                    genes = set(pd.read_csv(predictor / 'genes.csv').hgnc_id)
                    self.assertEqual(len(genes), count)
                    self.assertTrue(previous <= genes)
                    previous = genes
                    for role in ('train', 'validation', 'test'):
                        relative = f'tenv/fold_{fold}/{role}.csv'
                        self.assertEqual((predictor / relative).read_bytes(),
                                         (landmark / 'predictor' / relative).read_bytes())
                    self.assertFalse(list(predictor.rglob('test_expression.npy')))
                    actual = np.load(predictor / f'tenv/fold_{fold}/train_expression.npy', mmap_mode='r')
                    self.assertEqual(actual.shape[1], count)
            plan = root / 'plan.json'
            run('script/rerun_inputs.py', '--inputs', root, '--landmark', landmark,
                '--tenv-release', tenv, '--native', native, '--output', plan, '--fixture')
            payload = json.loads(plan.read_text())
            self.assertEqual(len(payload['runs']), 24)
            self.assertEqual(len({r['id'] for r in payload['runs']}), 24)
            self.assertEqual({r['fold'] for r in payload['runs']}, {0, 1})
            self.assertEqual({r['seed'] for r in payload['runs']}, {17, 42})

    def test_fold_one_panels_ignore_all_nontraining_expression(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            original = root / 'original'
            native, tenv = make_fixture(original)
            landmark = original / 'landmark'
            run('01_corpus/build.py', '--native', native, '--tenv-release', tenv,
                '--landmarks', original / 'landmarks.csv', '--tahoe-drugs', original / 'tahoe_drugs.parquet',
                '--output', landmark, '--fixture')
            first = root / 'first'
            run('01_corpus/gene_axes.py', '--native', native, '--landmark', landmark,
                '--output', first, '--fixture-sizes', 1, 2, 3, '--rerun', '--fold', 1)
            roles = pd.read_csv(tenv / 'folds.csv').query('fold == 1').set_index('drug_key').role
            for source in ('tahoe', 'tenv'):
                path = native / f'{source}.h5ad'
                data = ad.read_h5ad(path)
                values = data.X.toarray() if hasattr(data.X, 'toarray') else data.X.copy()
                if source == 'tahoe':
                    validation = min(data.obs.plate.unique(), key=lambda p: hashlib.sha256(f'20261110|{p}'.encode()).digest())
                    outside = data.obs.plate.eq(validation) & ~data.obs.is_control
                else:
                    outside = data.obs.drug_key.map(roles).isin(('validation', 'test'))
                rows = np.flatnonzero(outside)
                self.assertGreater(len(rows), 0)
                values[rows, 0] = np.arange(len(rows)) * 500 + 50
                data.X = values
                data.write_h5ad(path)
            changed_landmark = original / 'changed_landmark'
            run('01_corpus/build.py', '--native', native, '--tenv-release', tenv,
                '--landmarks', original / 'landmarks.csv', '--tahoe-drugs', original / 'tahoe_drugs.parquet',
                '--output', changed_landmark, '--fixture')
            second = root / 'second'
            run('01_corpus/gene_axes.py', '--native', native, '--landmark', changed_landmark,
                '--output', second, '--fixture-sizes', 1, 2, 3, '--rerun', '--fold', 1)
            for name in ('ranks', 'G512', 'G1024', 'G2048'):
                self.assertEqual((first / f'fold_1/{name}.csv').read_bytes(),
                                 (second / f'fold_1/{name}.csv').read_bytes())
