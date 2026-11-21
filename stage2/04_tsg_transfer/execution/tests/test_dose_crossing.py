"""Public design and parent-paired report checks for the dose crossing."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_crossing import crossing_splits, pair_scores


class DoseCrossingTest(unittest.TestCase):
    def test_same_held_parents_and_balanced_training_without_leaking_dose_or_parent(self):
        records, roles = [], {"train": [], "validation": []}
        for role, parents in (("train", ("a", "b")), ("validation", ("c", "d"))):
            for parent in parents:
                for dose, count in ((0.2, 2), (1.0, 3), (10.0, 4)):
                    for _ in range(count):
                        roles[role].append(len(records))
                        records.append((parent, dose))
        records.extend([("outer", .2), ("outer", 1.), ("outer", 10.)])
        data = SimpleNamespace(obs=pd.DataFrame(records, columns=["drug_parent_inchikey", "dose_um"]))
        splits = crossing_splits(data, **roles)
        held_sets = []
        for dose, split in splits.items():
            assessment = data.obs.iloc[split["assessment"]]
            held_sets.append(set(assessment.drug_parent_inchikey))
            self.assertTrue(assessment.dose_um.eq(dose).all())
            for role in ("train", "selection"):
                part = data.obs.iloc[split[role]]
                self.assertNotIn(dose, set(part.dose_um))
                self.assertNotIn("outer", set(part.drug_parent_inchikey))
                self.assertFalse(set(part.drug_parent_inchikey) & held_sets[-1])
            self.assertEqual(len(split["train"]), 8)
            self.assertEqual(len(split["selection"]), 4)
        self.assertTrue(all(held == held_sets[0] for held in held_sets))

    def test_pairing_takes_mean_profile_difference_on_common_finite_rows(self):
        scores = pd.DataFrame(dict(drug_parent_inchikey=["a", "a", "b"], dose_um=[1., 1., 1.],
                                   pcc=[.1, .3, .6], rmse=[1., 2., 3.],
                                   direction_error=[.2, .4, .5], amplitude_error=[1., 3., 2.]))
        transfer = scores.assign(pcc=[.3, np.nan, .8], rmse=[.5, 1., 2.],
                                 direction_error=[.1, .3, .1], amplitude_error=[.5, 2., 1.])
        paired = pair_scores(scores, transfer)
        self.assertEqual(int(paired.iloc[0].pcc_n), 1)
        self.assertAlmostEqual(paired.iloc[0].pcc_change, .2)
        self.assertAlmostEqual(paired.iloc[0].rmse_change, -.75)

    def test_mismatched_profile_id_is_rejected(self):
        scores = pd.DataFrame(dict(drug_parent_inchikey=["a"], dose_um=[1.],
                                   pcc=[.1], rmse=[1.], direction_error=[.2], amplitude_error=[1.]))
        with self.assertRaises(ValueError):
            pair_scores(scores, scores.set_axis([2]))


if __name__ == "__main__":
    unittest.main()
