"""Step-owned bounded predictor-input construction (#163).

These tests cover the two pure seams of the input release build: the binary
Morgan-2048 representation derived from the accepted corpus parent structures,
and the deterministic bounded well selection inside one declared fold.  The
corpus-reading call itself runs in the shared environment, because the public
``bases`` readers need its newer ``anndata``; the module stays importable here.
"""

import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'execution/script'))

from inputs import (  # noqa: E402
    CASE_FIELDS, MORGAN_SPEC, InputError, drug_fingerprints, role_selection)


class MorganRepresentation(unittest.TestCase):
    def setUp(self):
        self.obs = pd.DataFrame({
            'drug_key': ['aspirin', 'ethanol', 'benzene', 'ethanol'],
            'parent_smiles': ['CC(=O)Oc1ccccc1C(=O)O', 'CCO',
                              'c1ccccc1', 'CCO'],
        })

    def test_representation_is_the_declared_binary_morgan_2048(self):
        result = drug_fingerprints(self.obs)
        rows = result['drugs']
        morgan = result['morgan']
        self.assertEqual(MORGAN_SPEC['radius'], 2)
        self.assertEqual(MORGAN_SPEC['bits'], 2048)
        self.assertTrue(MORGAN_SPEC['chirality'])
        self.assertTrue(MORGAN_SPEC['binary'])
        self.assertEqual([row['drug_key'] for row in rows],
                         ['aspirin', 'benzene', 'ethanol'])
        self.assertEqual([row['morgan_index'] for row in rows], [0, 1, 2])
        self.assertEqual(morgan.shape, (3, 2048))
        self.assertEqual(morgan.dtype, np.uint8)
        self.assertTrue(np.isin(morgan, (0, 1)).all())
        self.assertGreater(int(morgan.sum()), 0)
        again = drug_fingerprints(self.obs)['morgan']
        self.assertTrue(np.array_equal(morgan, again))

    def test_a_drug_with_two_parent_structures_is_rejected(self):
        obs = pd.DataFrame({
            'drug_key': ['x', 'x'],
            'parent_smiles': ['CCO', 'c1ccccc1'],
        })
        with self.assertRaises(InputError):
            drug_fingerprints(obs)

    def test_an_invalid_structure_is_rejected(self):
        obs = pd.DataFrame({'drug_key': ['x'], 'parent_smiles': ['not-a-smiles']})
        with self.assertRaises(InputError):
            drug_fingerprints(obs)


class BoundedSelection(unittest.TestCase):
    def setUp(self):
        self.cohorts = pd.DataFrame({
            'sample_id': [f'plate-{index:02d}' for index in range(9)],
            'fold': [0, 0, 0, 0, 0, 0, 1, 1, 0],
            'role': ['train', 'test', 'train', 'validation', 'train',
                     'test', 'train', 'test', 'validation'],
        })

    def test_prefix_is_sorted_capped_and_fold_scoped(self):
        selected = role_selection(self.cohorts, 0, {'train': 2, 'test': 1,
                                                    'validation': 1})
        self.assertEqual(selected['train'], ['plate-00', 'plate-02'])
        self.assertEqual(selected['test'], ['plate-01'])
        self.assertEqual(selected['validation'], ['plate-03'])

    def test_limits_must_fit_the_available_wells(self):
        with self.assertRaises(InputError):
            role_selection(self.cohorts, 0, {'train': 99})
        with self.assertRaises(InputError):
            role_selection(self.cohorts, 7, {'train': 1})

    def test_unknown_roles_are_rejected(self):
        cohorts = self.cohorts.copy()
        cohorts.loc[0, 'role'] = 'holdout'
        with self.assertRaises(InputError):
            role_selection(cohorts, 0, {'train': 1})

    def test_case_fields_are_the_predictor_allowlist(self):
        self.assertEqual(CASE_FIELDS, ['case_id', 'drug_key', 'dose_um',
                                       'group_id', 'control_group',
                                       'source_row', 'morgan_index'])


if __name__ == '__main__':
    unittest.main()
