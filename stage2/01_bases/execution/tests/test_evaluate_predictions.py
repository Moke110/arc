"""Public evaluator CLI preserves the declared role's complete population."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))
from bases.significance import group_key, publish_statistics


class EvaluatorPopulationContract(unittest.TestCase):
    def test_cli_rejects_omitted_group_duplicate_and_foreign_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            groups = pd.DataFrame([
                dict(group_id=group_key(parent, 1), drug_parent=parent, dose_um=1,
                     status='eligible', reason='', sample_ids_sha256=hashlib.sha256(
                         json.dumps([sample], separators=(',', ':')).encode()).hexdigest())
                for parent, sample in [('A', 'a'), ('B', 'b'), ('C', 'c')]])
            samples = pd.DataFrame([
                dict(group_id=group_key(parent, 1), sample_id=sample,
                     source_profile_id='native-' + sample, condition='treated')
                for parent, sample in [('A', 'a'), ('B', 'b'), ('C', 'c')]])
            # Tiny row count on the real axis; hand-worked significant first gene.
            pvalue = np.ones((3, 18164)); pvalue[:, 0] = 1e-10
            padj = np.ones((3, 18164)); padj[:, 0] = 18164e-10
            publish_statistics(root/'truth', pd.DataFrame({
                'gene_index': range(18164), 'hgnc_id': [f'HGNC:{i}' for i in range(18164)]}),
                groups, samples, dict(pvalue=pvalue, padj=padj,
                                      log2_fold_change=np.ones((3, 18164))), {})
            pd.DataFrame({'sample_id':['a','b','c'], 'role':['test','test','train'],
                          'fold':[0,0,0]}).to_csv(root/'assignments.csv', index=False)
            pd.DataFrame({'gene_index':[0]}).to_csv(root/'genes.csv', index=False)
            for ids in (['a'], ['a','a','b'], ['a','b','c'], ['a','b','x'], ['b','a']):
                with self.subTest(ids=ids):
                    pd.DataFrame({'sample_id': ids, 'drug_parent': [x.upper() for x in ids],
                                  'dose_um': [1]*len(ids)}).to_csv(root/'samples.csv', index=False)
                    np.save(root/'predictions.npy', np.ones((len(ids), 1)))
                    output = root/('output-' + ''.join(ids))
                    result = subprocess.run([sys.executable, str(EXECUTION/'evaluate_predictions.py'),
                        '--predictions', str(root/'predictions.npy'), '--samples', str(root/'samples.csv'),
                        '--genes', str(root/'genes.csv'), '--truth', str(root/'truth'),
                        '--assignments', str(root/'assignments.csv'), '--fold', '0', '--role', 'test',
                        '--label', 'Contract example', '--output', str(output)], capture_output=True, text=True)
                    if ids == ['b','a']:
                        self.assertEqual(result.returncode, 0, result.stderr)
                        report = json.loads((output/'report.json').read_text())
                        self.assertEqual(report['summary']['groups'], 2)
                        self.assertEqual(report['summary']['AUPRC'], 1.)
                    else:
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn('exactly cover', result.stderr)
                        self.assertFalse(output.exists())
