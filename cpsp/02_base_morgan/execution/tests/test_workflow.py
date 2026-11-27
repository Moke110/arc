"""Acceptance tests at the public workflow subprocess seam."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = Path(__file__).resolve().parents[1] / 'workflow.py'


class PrepareAcceptance(unittest.TestCase):
    def test_cold_drug_manifest_is_grouped_complete_and_order_independent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = []
            for condition in ['baseline', 'H-R']:
                for parent, smiles in enumerate(['C', 'CC', 'CCC', 'CCCC', 'CCO', 'CCN', 'CO', 'CN', 'CCCl']):
                    for dose in [0.04, 0.2, 1, 5]:
                        for alias in range(2):
                            rows.append(dict(condition=condition, phenotype='cell_viability', canonical_smiles=smiles,
                                             dose_um=dose, direction=parent % 3 - 1, source_label=f'{parent}-{alias}',
                                             source_asset='fixture', source_cell=f'{parent}-{dose}-{alias}'))
            for reverse, output in [(False, 'a'), (True, 'b')]:
                source = root / f'input-{output}.csv'
                with source.open('w') as handle:
                    writer = csv.DictWriter(handle, fieldnames=rows[0])
                    writer.writeheader()
                    writer.writerows(reversed(rows) if reverse else rows)
                result = subprocess.run([sys.executable, str(CLI), 'prepare', '--input', str(source), '--output', str(root / output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
            a = json.loads((root / 'a' / 'prepared.json').read_text())
            b = json.loads((root / 'b' / 'prepared.json').read_text())
            self.assertEqual(a, b)
            self.assertEqual(len(a['rows']), len(rows))
            self.assertEqual(a['class_order'], [-1, 0, 1])
            tested = []
            for fold in a['splits']:
                roles = fold['roles']
                parents = {role: {a['rows'][i]['parent_id'] for i, r in enumerate(roles) if r == role} for role in ['train', 'validation', 'test']}
                self.assertFalse(parents['train'] & parents['test'])
                self.assertFalse(parents['train'] & parents['validation'])
                self.assertFalse(parents['validation'] & parents['test'])
                tested.extend(parents['test'])
                for row, role in zip(a['rows'], roles):
                    self.assertEqual(role, fold['parent_roles'][row['parent_id']])
            self.assertEqual(len(tested), 9)
            self.assertEqual(len(set(tested)), 9)
            self.assertEqual(len(a['expected_runs']), 18)
            # Existing releases cannot be silently replaced, even when their file is damaged.
            (root / 'a' / 'prepared.json').write_text('{}')
            result = subprocess.run([sys.executable, str(CLI), 'prepare', '--input', str(root / 'input-a.csv'), '--output', str(root / 'a')], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('incompatible', result.stderr)
            rows[0]['dose_um'] = 0
            with (root / 'invalid.csv').open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0])
                writer.writeheader()
                writer.writerows(rows)
            result = subprocess.run([sys.executable, str(CLI), 'prepare', '--input', str(root / 'invalid.csv'), '--output', str(root / 'invalid')], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('positive and finite', result.stderr)
            self.assertFalse((root / 'invalid').exists())


    def test_public_smoke_roundtrip_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [dict(condition=condition, phenotype='cell_viability', canonical_smiles=smiles,
                         dose_um=dose, direction=parent % 3 - 1, source_label=f'{parent}-{alias}',
                         source_asset='fixture', source_cell=f'{parent}-{dose}-{alias}')
                    for condition in ['baseline', 'H-R']
                    for parent, smiles in enumerate(['C', 'CC', 'CCC', 'CCCC', 'CCO', 'CCN', 'CO', 'CN', 'CCCl'])
                    for dose in [0.04, 0.2, 1, 5] for alias in range(2)]
            source = root / 'fixture.csv'
            with source.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0])
                writer.writeheader()
                writer.writerows(rows)
            def invoke(*args, ok=True):
                result = subprocess.run([sys.executable, str(CLI), *map(str, args)], capture_output=True, text=True)
                if ok:
                    self.assertEqual(result.returncode, 0, result.stderr)
                else:
                    self.assertNotEqual(result.returncode, 0)
                return result
            prepared = root / 'prepared'
            run = root / 'run'
            invoke('prepare', '--input', source, '--output', prepared)
            invoke('train', '--prepared', prepared / 'prepared.json', '--output', run,
                   '--condition', 'baseline', '--phenotype', 'cell_viability', '--setting', 'cold_drug',
                   '--fold', 0, '--seed', 17, '--smoke-epochs', 2)
            payload = json.loads((prepared / 'prepared.json').read_text())
            doses = [math.log10(row['dose_um']) for row, role in zip(payload['rows'], payload['splits'][0]['roles'])
                     if role == 'train' and row['condition'] == 'baseline']
            statistics = json.loads((run / 'run.json').read_text())['dose_statistics']
            self.assertAlmostEqual(statistics['mean'], sum(doses) / len(doses))
            self.assertAlmostEqual(statistics['std'], (sum((d-statistics['mean'])**2 for d in doses) / len(doses))**0.5)
            invoke('train', '--prepared', prepared / 'prepared.json', '--output', root / 'hr-run',
                   '--condition', 'H-R', '--phenotype', 'cell_viability', '--setting', 'cold_drug',
                   '--fold', 0, '--seed', 17, '--smoke-epochs', 2)
            invoke('verify', '--run', root / 'hr-run', '--prepared', prepared / 'prepared.json')
            first = json.loads(invoke('evaluate', '--run', run).stdout)
            self.assertGreater(first['n'], 0)
            self.assertEqual(first, json.loads(invoke('evaluate', '--run', run).stdout))
            invoke('verify', '--run', run, '--prepared', prepared / 'prepared.json')
            payload = json.loads((prepared / 'prepared.json').read_text())
            payload['rows'][0]['direction'] = 1 - payload['rows'][0]['direction']
            (prepared / 'prepared.json').write_text(json.dumps(payload))
            invoke('verify', '--run', run, '--prepared', prepared / 'prepared.json', ok=False)
            predictions = run / 'predictions.json'
            records = json.loads(predictions.read_text())
            records[0]['probabilities'] = [1, 1, 1]
            predictions.write_text(json.dumps(records))
            invoke('verify', '--run', run, ok=False)
            invoke('evaluate', '--run', run, ok=False)

    def test_directory_requires_verified_corpus_and_impossible_split_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake = root / 'fake'
            fake.mkdir()
            result = subprocess.run([sys.executable, str(CLI), 'prepare', '--input', str(fake),
                                     '--output', str(root / 'out')], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Frozen corpus verification failed', result.stderr)
            source = root / 'small.csv'
            with source.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=['condition', 'phenotype', 'canonical_smiles',
                                        'dose_um', 'direction', 'source_asset', 'source_cell'])
                writer.writeheader()
                writer.writerow(dict(condition='baseline', phenotype='cell_viability', canonical_smiles='C',
                                     dose_um=1, direction=0, source_asset='fixture', source_cell='1'))
            result = subprocess.run([sys.executable, str(CLI), 'prepare', '--input', str(source),
                                     '--output', str(root / 'out')], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('at least six parents', result.stderr)



if __name__ == '__main__':
    unittest.main()
