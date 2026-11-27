from pathlib import Path
import sys
import unittest

import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/05_tsg_loss/similarity"))

from report_similarity import METRICS, verify_reference_scores


class FormalReferenceTests(unittest.TestCase):
    def setUp(self):
        self.reference = pd.DataFrame([
            {"alpha": "1", "fold": fold, "seed": 42, "metric": metric, "value": 0.1,
             "n_conditions": 10, "n_eligible": 9, "n_undefined": 1}
            for fold in range(3) for metric in METRICS
        ])
        self.fresh = self.reference.assign(alpha="0")

    def test_exact_reference_scores_accept_all_three_folds(self):
        self.assertEqual(verify_reference_scores(self.fresh, self.reference), 27)

    def test_changed_reference_eligibility_fails(self):
        self.fresh.loc[0, ["n_eligible", "n_undefined"]] = [8, 2]
        with self.assertRaisesRegex(ValueError, "eligibility"):
            verify_reference_scores(self.fresh, self.reference)


if __name__ == "__main__":
    unittest.main()
