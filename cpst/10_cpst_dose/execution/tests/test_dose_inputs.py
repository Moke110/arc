"""Verify the crossed Step 10 row boundary through the frozen ledger."""

from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_inputs import cold_dose_partitions
from dose_protocol import protocol


class PartitionTests(unittest.TestCase):
    def test_crossed_holdouts_and_excluded_strips(self):
        rows = pd.DataFrame(
            [(row.plate, row.drug_parent, row.dose_um) for row, _ in protocol().folds[0].partitions],
            columns=["plate", "drug_parent", "dose_um"],
        )
        tests = []
        for index, fold in enumerate(protocol().folds):
            with self.subTest(index=index):
                splits = cold_dose_partitions(rows, index)
                self.assertEqual(sum(map(len, splits.values())), 1383)
                held_parents = rows.drug_parent.isin(fold.held_parents)
                held_dose = rows.dose_um.eq((0.2, 1.0, 10.0)[index % 3])
                for label, mask in (("test", held_parents & held_dose),
                                    ("excluded-held-parent", held_parents & ~held_dose),
                                    ("excluded-held-dose", ~held_parents & held_dose)):
                    self.assertEqual(splits[label].tolist(), rows.index[mask].tolist())
                development = sorted([*splits["train"], *splits["valid"]])
                self.assertEqual(development, rows.index[~held_parents & ~held_dose].tolist())
                tests.extend(splits["test"])
        self.assertEqual(sorted(tests), list(rows.index))
