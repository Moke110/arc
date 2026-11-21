import sys
from pathlib import Path
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from publish_report import IDENTITY, parent_pcc, paired_pcc


class ReportEstimands(unittest.TestCase):
    def test_primary_weights_parents_equally(self):
        scores = pd.DataFrame({"fold": [0, 0, 0], "drug_parent": ["A", "A", "B"],
                               "pcc": [.5, .5, -.5]})
        parents = parent_pcc(scores)
        self.assertEqual(parents.pcc.tolist(), [.5, -.5])
        self.assertEqual(parents.pcc.mean(), 0.)

    def test_paired_effect_uses_only_defined_matches(self):
        baseline = pd.DataFrame([(0, "1", "A", .2, "L", .5),
                                 (0, "2", "A", .2, "L", None),
                                 (0, "3", "B", .2, "L", -.5)],
                                columns=[*IDENTITY, "pcc"])
        arm = baseline.assign(pcc=[.75, -.75, 0.])
        pairs = paired_pcc(arm, baseline)
        self.assertEqual(pairs.delta_pcc.tolist(), [.25, .5])
        self.assertEqual(pairs.defined_pairs.tolist(), [1, 1])

    def test_pairing_rejects_changed_held_out_identity(self):
        baseline = pd.DataFrame([(0, "1", "A", .2, "L", .5)],
                                columns=[*IDENTITY, "pcc"])
        with self.assertRaisesRegex(ValueError, "held-out"):
            paired_pcc(baseline.assign(drug_parent="B"), baseline)
