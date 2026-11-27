from pathlib import Path
import sys
import unittest

import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
EXECUTION = REPO / "stage2_CarPS-VAE/06_tsg_shape/06.2_round2/execution"
sys.path.insert(0, str(EXECUTION))
from prior_architecture import verify_original_scores

PUBLISHED = REPO / "stage2_CarPS-VAE/06_tsg_shape/06.1_round1/figure/data/fold-scores.csv"


class FormalArchitectureTest(unittest.TestCase):
    def test_formal_report_rejects_changed_original_eligibility(self):
        scores = pd.read_csv(PUBLISHED, float_precision="round_trip")
        original = scores[scores["shape"] == "original"].copy()
        verify_original_scores(original, scores)
        original.loc[original.index[0], ["n_eligible", "n_undefined"]] = [461, 1]
        with self.assertRaisesRegex(ValueError, "original scores or eligibility"):
            verify_original_scores(original, scores)
