"""Warm-start split acceptance through prepared public artifacts."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from warm_start import warm_start

CLI = Path(__file__).resolve().parents[1] / 'workflow.py'


class WarmStartAcceptance(unittest.TestCase):
    def prepare(self, root, rows, name='prepared', ok=True):
        source = root / (name + '.csv')
        with source.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
        result = subprocess.run([sys.executable, str(CLI), 'prepare', '--input', str(source),
                                 '--output', str(root / name), '--settings', 'warm_start'],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, ok, result.stderr)
        return json.loads((root / name / 'prepared.json').read_text()) if ok else result.stderr

    def fixture(self):
        return [dict(condition=condition, phenotype='cell_viability', canonical_smiles=smiles,
                     dose_um=dose, direction=parent % 3 - 1, source_label=f'{parent}-{alias}',
                     source_asset='fixture', source_cell=f'{parent}-{dose}-{alias}')
                for condition in ['baseline', 'H-R']
                for parent, smiles in enumerate(['C', 'CC', 'CCC', 'CCCC', 'CCO', 'CCN'])
                for dose in [0.04, 0.2, 1, 5] for alias in range(2)]

    def test_all_combinations_tested_once_with_seen_training_marginals(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = self.fixture()
            payload = self.prepare(root, rows)
            self.assertEqual(payload, self.prepare(root, list(reversed(rows)), 'reversed'))
            self.assertEqual(len(payload['rows']), len(rows))
            tested = []
            for split in payload['splits']:
                groups = {role: {(r['parent_id'], r['dose_um']) for r, assigned in
                                 zip(payload['rows'], split['roles']) if assigned == role}
                          for role in ['train', 'validation', 'test']}
                self.assertEqual(len(groups['validation']), math.ceil(.2 * (24 - len(groups['test']))))
                self.assertFalse(groups['test'] & groups['validation'])
                self.assertFalse(groups['train'] & (groups['test'] | groups['validation']))
                for role in ['validation', 'test']:
                    for axis in [0, 1]:
                        self.assertLessEqual({g[axis] for g in groups[role]}, {g[axis] for g in groups['train']})
                mapping = {(r['parent_id'], r['dose_um']): r['role'] for r in split['combination_roles']}
                for row, role in zip(payload['rows'], split['roles']):
                    self.assertEqual(role, mapping[row['parent_id'], row['dose_um']])
                self.assertEqual(sum(c['measurements'] for c in split['counts'].values()), len(rows))
                tested.extend(groups['test'])
            self.assertEqual(len(tested), 24)
            self.assertEqual(len(set(tested)), 24)

    def test_impossible_parent_coverage_fails_without_dropping_measurements(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.fixture()
            rows = [r for r in rows if r['canonical_smiles'] != 'C' or r['dose_um'] == .04]
            message = self.prepare(Path(tmp), rows, ok=False)
            self.assertIn('warm_start', message)
            self.assertIn('parent', message)
            self.assertIn('coverage', message)
            self.assertFalse((Path(tmp) / 'prepared').exists())

    def test_target_specific_singleton_cannot_borrow_training_from_other_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = [r for r in self.fixture()
                    if r['condition'] != 'H-R' or r['canonical_smiles'] != 'C'
                    or r['dose_um'] == .04]
            message = self.prepare(Path(tmp), rows, ok=False)
            self.assertIn('warm_start impossible parent coverage', message)
            self.assertIn('H-R', message)
            self.assertFalse((Path(tmp) / 'prepared').exists())

    def test_joint_coverage_can_be_impossible_despite_non_singleton_marginals(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = [r for r in self.fixture() if r['canonical_smiles'] in ['C', 'CC']
                    and r['dose_um'] in [.04, .2]]
            message = self.prepare(Path(tmp), rows, ok=False)
            self.assertIn('warm_start impossible', message)
            self.assertIn('20% validation', message)

    def test_dense_release_sized_grid_needs_no_feasibility_solver(self):
        rows = [dict(parent_id=f'p{parent:04}', dose_um=dose, condition=condition,
                     phenotype='cell_viability', direction=0)
                for parent in range(663) for dose in [.04, .2, 1, 5]
                for condition in ['baseline', 'H-R']]
        with patch('warm_start._solve', side_effect=AssertionError('dense grid invoked solver')):
            splits = warm_start(rows)
            reversed_splits = warm_start(list(reversed(rows)))
        tested = []
        for split, reverse in zip(splits, reversed_splits):
            self.assertEqual(split['combination_roles'], reverse['combination_roles'])
            combinations = split['combination_roles']
            tests = [r for r in combinations if r['role'] == 'test']
            validation = [r for r in combinations if r['role'] == 'validation']
            self.assertEqual(len(validation), math.ceil(.2 * (2652 - len(tests))))
            for axis in ['parent_id', 'dose_um']:
                train = {r[axis] for r in combinations if r['role'] == 'train'}
                self.assertEqual(train, {r[axis] for r in combinations})
            tested.extend((r['parent_id'], r['dose_um']) for r in tests)
        self.assertEqual(len(tested), 2652)
        self.assertEqual(len(set(tested)), 2652)

    def test_public_train_evaluate_verify_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = self.prepare(root, self.fixture())
            run = root / 'run'
            def invoke(*args, ok=True):
                result = subprocess.run([sys.executable, str(CLI), *map(str, args)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode == 0, ok, result.stderr)
                return result
            invoke('train', '--prepared', root / 'prepared/prepared.json', '--output', run,
                   '--condition', 'baseline', '--phenotype', 'cell_viability', '--setting', 'warm_start',
                   '--fold', 0, '--seed', 17, '--smoke-epochs', 2)
            invoke('evaluate', '--run', run)
            invoke('verify', '--run', run, '--prepared', root / 'prepared/prepared.json')
            statistics = json.loads((run / 'run.json').read_text())['dose_statistics']
            doses = [math.log10(r['dose_um']) for r, role in zip(payload['rows'], payload['splits'][0]['roles'])
                     if role == 'train' and r['condition'] == 'baseline']
            self.assertAlmostEqual(statistics['mean'], sum(doses) / len(doses))
            self.assertAlmostEqual(statistics['std'], (sum((d-statistics['mean'])**2 for d in doses) / len(doses))**.5)
            prepared = root / 'prepared/prepared.json'
            payload['splits'][0]['roles'][0] = 'corrupted'
            prepared.write_text(json.dumps(payload))
            invoke('verify', '--run', run, '--prepared', prepared, ok=False)
