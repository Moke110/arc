"""Bounded fixture for the Step 06 source-audit table interface."""
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from source_audit import PATHS, Release, read_json, rows, sha256, write_audit, write_csv


class SourceAuditTest(unittest.TestCase):
    def fixture_release(self, root, members, expressions):
        root.mkdir()
        drugs = {}
        for role, entries in members.items():
            cases = []
            for number, (parent, dose, line, smiles) in enumerate(entries):
                drugs[parent] = {'drug_key': parent, 'parent_connectivity': parent,
                                 'parent_smiles': smiles}
                cases.append(dict(case_id=f'{role}-{number}', drug_key=parent,
                                  dose_um=dose, cell_line_id=line,
                                  group_id=f'{parent}-{dose}-{line}'))
            write_csv(root / PATHS[f'{role}_cases'], cases)
        write_csv(root / PATHS['drugs'], drugs.values())
        write_csv(root / PATHS['genes'], [{'symbol': 'A'}, {'symbol': 'B'}, {'symbol': 'C'}])
        for role, values in expressions.items():
            np.save(root / PATHS[f'{role}_expression'], np.asarray(values, dtype=np.float32))
        return Release(root, {key: sha256(root / path) for key, path in PATHS.items()
                              if (root / path).exists()})

    def test_separates_fitting_roles_for_historical_and_cold_exposure(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            a, b, v, t = [(parent, 0.2, 'CM', parent) for parent in 'ABVT']
            tenv = self.fixture_release(
                root / 'tenv', {'train': [a, b], 'validation': [v], 'test': [t]},
                {'train': [[1, -1, 0], [1, 1, 0]]})
            a, b, v, t = [(parent, 0.05, 'L1', parent) for parent in 'ABVT']
            historical = self.fixture_release(root / 'historical',
                                              {'train': [a, v], 'validation': [b, t]}, {})
            cold = self.fixture_release(root / 'cold', {'train': [a], 'validation': [b]},
                                        {'train': [[2, -2, 0]], 'validation': [[-1, -1, 0]]})
            output = root / 'audit'
            write_audit(output, tenv, historical, cold, 3)
            exposures = {(row['source_release'], row['tahoe_role'], row['tenv_role']): row
                         for row in rows(output / 'exposure.csv')}
            self.assertEqual(len(exposures), 12)
            self.assertEqual(exposures['historical', 'train', 'validation']['shared_parents'], 'V')
            self.assertEqual(exposures['historical', 'validation', 'validation']['shared_parent_count'], '0')
            self.assertEqual(exposures['historical', 'validation', 'test']['shared_parents'], 'T')
            for role in ('train', 'validation'):
                for target in ('validation', 'test'):
                    self.assertEqual(exposures['globally_cold', role, target]['shared_parent_count'], '0')
            self.assertEqual(exposures['globally_cold', 'validation', 'train']['shared_parents'], 'B')

    def test_keeps_observed_brackets_line_role_and_exact_smiles_sensitivity(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            tenv = self.fixture_release(
                root / 'tenv',
                {'train': [(parent, dose, 'CM', parent) for parent in 'AB'
                           for dose in (0.2, 1, 10)],
                 'validation': [('V', 0.2, 'CM', 'V')], 'test': [('T', 0.2, 'CM', 'T')]},
                {'train': [[3, -3, 0]] * 6})
            membership = {
                'train': [('A', 0.05, 'L1', 'A'), ('A', 0.5, 'L1', 'A'), ('A', 5, 'L1', 'A')],
                'validation': [('B', 0.5, 'L2', 'B-other')],
            }
            historical = self.fixture_release(root / 'historical', membership, {})
            cold = self.fixture_release(
                root / 'cold', membership,
                {'train': [[2, 2, 3], [-2, 2, 0], [2, -2, 0]], 'validation': [[0, 0, 0]]})
            output = root / 'audit'
            write_audit(output, tenv, historical, cold, 3)
            panel = rows(output / 'delta_panel.csv')
            sensitivity = rows(output / 'exact_smiles_panel.csv')
            cross = [row for row in panel if row['parent_connectivity'] == 'A' and row['tahoe_role']]
            self.assertEqual([(row['tenv_dose_um'], row['tahoe_dose_um'], row['bracket']) for row in cross],
                             [('0.2', '0.05', 'lower'), ('0.2', '0.5', 'upper'),
                              ('1.0', '0.5', 'lower'), ('1.0', '5.0', 'upper')])
            self.assertTrue(all(row['tahoe_cell_line_id'] == 'L1' and row['tahoe_role'] == 'train'
                                for row in cross))
            self.assertAlmostEqual(float(cross[0]['tenv_delta_rms']), 6 ** 0.5)
            self.assertEqual(cross[0]['sign_concordance'], '0.5')
            self.assertEqual(cross[0]['joint_nonzero_genes'], '2')
            self.assertEqual(cross[1]['sign_concordance'], '0.0')
            self.assertEqual(cross[3]['sign_concordance'], '1.0')
            self.assertEqual({row['parent_connectivity'] for row in sensitivity}, {'A'})
            b_upper = next(row for row in panel if row['parent_connectivity'] == 'B'
                           and row['tenv_dose_um'] == '0.2' and row['bracket'] == 'upper')
            self.assertEqual(b_upper['status'], 'no_joint_nonzero_genes')
            self.assertEqual(b_upper['tahoe_role'], 'validation')
            self.assertEqual(b_upper['tahoe_cell_line_id'], 'L2')
            b_lower = next(row for row in panel if row['parent_connectivity'] == 'B'
                           and row['tenv_dose_um'] == '0.2' and row['bracket'] == 'lower')
            self.assertEqual(b_lower['status'], 'missing_tahoe_dose')
            ten_only = [row for row in panel if row['tenv_dose_um'] == '10.0']
            self.assertEqual(len(ten_only), 2)
            self.assertTrue(all(row['status'] == 'tenv_only' and row['tahoe_dose_um'] == ''
                                and row['sign_concordance'] == '' for row in ten_only))
            evidence = read_json(output / 'lineage.json')
            self.assertEqual(evidence['outputs']['delta_panel.csv'], sha256(output / 'delta_panel.csv'))
