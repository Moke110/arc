"""Exercise the exact KPGT score gate used before reference fit reuse."""

from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chem_reference import verify_reference_scores
from chem_report import METRICS


class ReferenceTest(unittest.TestCase):
    def test_changed_score_rejects_reuse(self):
        recorded = pd.DataFrame({
            "metric": METRICS, "value": [.3] * len(METRICS),
            "n_conditions": 7, "n_eligible": 7, "n_undefined": 0,
            "run_sha256": "original-registration",
        })
        independently_scored = recorded.drop(columns="run_sha256").copy()
        independently_scored.loc[independently_scored.metric.eq("Δ PCC"), "value"] = .30000000000000004
        with self.assertRaisesRegex(ValueError, "exactly"):
            verify_reference_scores(independently_scored, recorded, "original-registration")


if __name__ == "__main__":
    unittest.main()
