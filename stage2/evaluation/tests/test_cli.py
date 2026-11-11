"""Every step's public reporting CLI must use the same three-metric scorer."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
from stage2.evaluation.evaluate import membership
from stage2.evaluation.inventory import REPO, binding, sha
from stage2.evaluation.metrics import METRICS


class PublicCLITest(unittest.TestCase):
    def test_all_step_reporting_commands_use_significance_not_rank_labels(self):
        wrappers = {3: '03_modes_tenv/execution/script/benchmark.py',
                    4: '04_pdiff_hyperparameter/execution/script/evaluate.py',
                    5: '05_pdiff_transfer/execution/script/publish.py',
                    6: '06_pdiff_transfer_optimize/execution/script/paired_publication.py',
                    9: '09_xpert_optimize/execution/script/evaluate.py'}
        with tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests') as temporary:
            root = Path(temporary)
            labels = root / 'labels/tenv'
            labels.mkdir(parents=True)
            genes = pd.DataFrame({'hgnc_id': [f'HGNC:{i}' for i in range(18165)],
                                  'symbol': [f'G{i}' for i in range(18164)] + ['HDGFL3']})
            genes.to_csv(labels / 'genes.csv', index=False)
            pd.DataFrame({'group_id': ['g', 'z'], 'status': ['eligible', 'eligible'],
                          'case_ids_sha256': [membership(['a', 'b']), membership(['c'])]}).to_csv(labels / 'groups.csv', index=False)
            padj = np.ones((2, 18165)); padj[0, :3] = .001
            direction = np.zeros_like(padj); direction[0, :3] = [1, -1, 1]
            np.save(labels / 'padj.npy', padj); np.save(labels / 'log2_fold_change.npy', direction)
            manifest = {'schema': 'stage2-significance-v1', 'fixture': True, 'fdr': .05,
                        'files': {p.name: sha(p) for p in labels.iterdir()}}
            (labels / 'manifest.json').write_text(json.dumps(manifest))
            genes.iloc[:2].to_csv(root / 'genes.csv', index=False)
            pd.DataFrame({'case_id': ['a', 'b', 'c'], 'group_id': ['g', 'g', 'z'],
                          'dose_um': [1., 1., 1.]}).to_csv(root / 'cases.csv', index=False)
            np.save(root / 'prediction.npy', np.array([[2., -1.], [2., -1.], [4., -3.]], dtype='float32'))
            runs = [{'step': step, 'dataset': 'tenv', 'regime': 'cold-drug', 'role': 'test', 'model': 'covered',
                     'fold': 0, 'seed': 17, 'cases': binding(root / 'cases.csv'),
                     'genes': binding(root / 'genes.csv'), 'prediction': binding(root / 'prediction.npy')}
                    for step in wrappers]
            inventory = root / 'inventory.json'; inventory.write_text(json.dumps({'runs': runs}))
            for step, script in wrappers.items():
                with self.subTest(step=step):
                    destination = root / f'step{step}'
                    args = ['--inventory', str(inventory), '--labels-root', str(root / 'labels'),
                            '--destination', str(destination)]
                    for action in ('score', 'replay'):
                        run = subprocess.run([sys.executable, str(REPO / 'stage2' / script), action, *args],
                                             capture_output=True, text=True)
                        self.assertEqual(run.returncode, 0, run.stderr)
                    result = pd.read_csv(destination / 'per_run.csv').iloc[0]
                    self.assertAlmostEqual(result.AUPRC, 2 / 3)
                    self.assertAlmostEqual(result.AUPRC_positive, .5)
                    self.assertEqual(result.AUPRC_negative, 1.)
                    self.assertEqual(result.groups, 2)
                    self.assertEqual(result.AUPRC_defined_groups, 1)
                    self.assertEqual(result.AUPRC_undefined_groups, 1)
                    points = pd.read_csv(destination / 'pr.csv')
                    self.assertEqual(set(points.metric), set(METRICS))
                    self.assertAlmostEqual(points[points.metric.eq('AUPRC')].recall.iloc[-1], 2 / 3)
