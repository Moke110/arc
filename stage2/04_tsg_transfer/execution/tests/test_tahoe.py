"""Public source fold construction and retained training coverage."""
from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))

from tsg_transfer.tahoe import make_tahoe_splits


class TahoeSplitsTest(unittest.TestCase):
    def test_replicates_are_atomic_and_rare_conditions_remain_in_training(self):
        regular = [(parent, line, dose)
                   for parent in ('P1', 'P2', 'P3')
                   for line in ('C1', 'C2')
                   for dose in (0.05, .5)
                   for _ in range(2)]
        obs = pd.DataFrame([*regular, ('P4', 'C1', 0.05)],
                           columns=['drug_parent_inchikey', 'cell_line', 'dose_um'],
                           index=[f'profile-{i}' for i in range(len(regular) + 1)])
        actual = make_tahoe_splits(obs)
        key = ['fold', 'profile_id']
        shuffled = make_tahoe_splits(obs.sample(frac=1, random_state=99))
        pd.testing.assert_frame_equal(actual.set_index(key).drop(columns='row').sort_index(),
                                      shuffled.set_index(key).drop(columns='row').sort_index())
        held = actual[actual.partition == 'validation']
        self.assertEqual(set(held.fold), {0, 1, 2})
        self.assertFalse(held.profile_id.duplicated().any())
        self.assertEqual(set(held.profile_id), set(obs.index[:-1]))
        for _, fold in actual.groupby('fold'):
            self.assertEqual(set(fold.profile_id), set(obs.index))
            self.assertEqual(len(fold), len(obs))
            self.assertTrue((fold.groupby(['drug_parent', 'cell_line', 'dose_um'])
                             .partition.nunique() == 1).all())
            for column in ['drug_parent', 'cell_line', 'dose_um']:
                self.assertLessEqual(set(fold.loc[fold.partition == 'validation', column]),
                                     set(fold.loc[fold.partition == 'train', column]))
