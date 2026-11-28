"""Verify the released cold-dose-drug row boundary and its two strips."""

from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_inputs import cold_dose_partitions


class PartitionTests(unittest.TestCase):
    def test_holdout_intersection_and_strips(self):
        rows = pd.DataFrame({
            "drug_parent": ["held", "held", "train", "train", "validate"],
            "dose_um": [0.2, 1, 0.2, 1, 10],
            "cold_drug_0": ["test", "test", "train", "train", "valid"],
            "cold_dose_drug_0": ["test", "", "", "train", "valid"],
            "cold_dose_drug_0_reason": ["", "held_parent", "held_dose", "", ""],
        })
        partitions = cold_dose_partitions(rows, 0)
        self.assertEqual({key: value.tolist() for key, value in partitions.items()},
                         {"train": [3], "valid": [4], "test": [0]})
        rows.loc[2, "cold_dose_drug_0_reason"] = "held_parent"
        with self.assertRaisesRegex(ValueError, "exclusion"):
            cold_dose_partitions(rows, 0)


if __name__ == "__main__":
    unittest.main()
