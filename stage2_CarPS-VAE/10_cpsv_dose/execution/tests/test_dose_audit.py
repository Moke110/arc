"""Check independently audited dose contrasts at the formal report seam."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from audit_dose_report import audit_aggregates


class AuditTests(unittest.TestCase):
    def test_rejects_cross_seed_pairing_even_when_absolute_scores_are_unchanged(self):
        metric = "Δ PCC"
        scores = {
            (arm, 0, seed, metric): {
                "value": observed, "reference_value": reference, "paired_difference": difference,
                "n_eligible": "2", "n_undefined": "0",
            }
            for arm, seed, observed, reference, difference in (
                ("no-dose", 17, "0.25", "0.25", "0"),
                ("no-dose", 42, "0.5", "0.5", "0"),
                ("Xmol-FiLM", 17, "0.5", "0.5", "0"),
                ("Xmol-FiLM", 42, "0.75", "0.25", "0.5"),
            )
        }
        with self.assertRaisesRegex(ValueError, "fold and seed"):
            audit_aggregates(scores, {}, {})


if __name__ == "__main__":
    unittest.main()
