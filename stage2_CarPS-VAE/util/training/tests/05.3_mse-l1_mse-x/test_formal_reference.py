"""Exercise the formal acceptance boundary for matching held-out fold scores."""

from itertools import product
from pathlib import Path
import sys
from unittest import TestCase

import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2_CarPS-VAE/05_tsg_loss/05.3_mse-l1_mse-x/execution"))

from report_reconstruction import METRICS, verify_reference_scores


class FormalReferenceTests(TestCase):
    def setUp(self):
        self.accepted = pd.DataFrame([
            {"alpha": "1.75", "fold": fold, "seed": 42, "metric": metric, "value": value,
             "n_conditions": 16, "n_eligible": 14, "n_undefined": 2}
            for fold, (metric, value) in product(range(3), zip(METRICS, range(9), strict=True))
        ])
        self.fresh = self.accepted.assign(alpha="1")

    def test_reference_requires_all_exact_fold_scores(self):
        self.assertEqual(verify_reference_scores(self.fresh, self.accepted), 27)
        with self.assertRaisesRegex(ValueError, "fold scores"):
            verify_reference_scores(self.fresh.iloc[:-1], self.accepted.iloc[:-1])

    def test_reference_rejects_different_eligibility(self):
        self.fresh.loc[0, ["n_eligible", "n_undefined"]] = [15, 1]
        with self.assertRaisesRegex(ValueError, "eligibility"):
            verify_reference_scores(self.fresh, self.accepted)
