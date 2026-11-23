"""KPGT input identity and exact parent alignment at the preparation boundary."""

import unittest
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tsg_kpgt.inputs import align_vectors


class ParentAlignment(unittest.TestCase):
    def test_aligns_inchikey_and_smiles_without_changing_vector_values(self):
        columns = ['drug_parent_inchikey', 'drug_parent_smiles']
        expected = pd.DataFrame([('A', 'CO'), ('B', 'CC')], columns=columns)
        supplied = expected.iloc[::-1].reset_index(drop=True)
        values = np.stack([np.full(2304, 5), np.full(2304, 7)]).astype(np.float32)
        np.testing.assert_array_equal(
            align_vectors(expected, supplied, values), values[::-1]
        )
        mismatch = supplied.assign(drug_parent_smiles=['CO', 'CC'])
        for bad in (mismatch, supplied.iloc[[0, 0]], supplied.iloc[:1]):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    align_vectors(expected, bad, values)
        values[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'finite'):
            align_vectors(expected, supplied, values)
