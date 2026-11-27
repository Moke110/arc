"""Public formal matrix contracts around coordinate selection and release reads."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axis_formal import fresh_coordinates


class FormalRoster(unittest.TestCase):
    def test_failed_reference_is_fitted_while_successful_replays_are_reused(self):
        decisions = [
            {"coordinate": ["denoise-5k", fold, seed], "decision": "reuse",
             "validation_replay": "exact", "test_replay": "exact"}
            for fold in range(3) for seed in (17, 42, 2026)
        ]
        decisions[3] = {"coordinate": ["denoise-5k", 1, 17], "decision": "refit",
                        "reason": "checkpoint replay differs"}
        remaining = fresh_coordinates({"decisions": decisions})
        self.assertEqual(len(remaining), 64)
        self.assertIn(("denoise-5k", 1, 17), remaining)
        self.assertNotIn(("denoise-5k", 1, 42), remaining)
        with self.assertRaises(ValueError):
            fresh_coordinates({"decisions": decisions[:-1] + decisions[:1]})


class ScientificAudit(unittest.TestCase):
    def test_audit_rejects_weighted_coverage_and_interchanged_pairs(self):
        from audit_axis_report import audit_aggregates

        scores = {
            ("hvg-2k", 0, 17, "disease", "direction-correct coverage"):
                dict(value="25", reference_value="17", paired_difference="8", overlap="1258",
                     n_conditions="1", n_eligible="1", n_undefined="0"),
            ("hvg-2k", 1, 17, "disease", "direction-correct coverage"):
                dict(value="16", reference_value="18", paired_difference="-2", overlap="1258",
                     n_conditions="5", n_eligible="5", n_undefined="0"),
            ("denoise-5k", 0, 17, "disease", "direction-correct coverage"):
                dict(value="17", reference_value="17", paired_difference="0", overlap="1513",
                     n_conditions="1", n_eligible="1", n_undefined="0"),
            ("denoise-5k", 1, 17, "disease", "direction-correct coverage"):
                dict(value="18", reference_value="18", paired_difference="0", overlap="1513",
                     n_conditions="5", n_eligible="5", n_undefined="0"),
        }
        summary = {
            ("hvg-2k", "disease", "direction-correct coverage"):
                dict(mean="20.5", sample_sd="6.363961030678928", n_defined="2", paired_mean="3",
                     paired_sd="7.0710678118654755", n_paired="2", n_coordinates="2",
                     n_eligible="6", n_undefined="0", overlap="1258"),
            ("denoise-5k", "disease", "direction-correct coverage"):
                dict(mean="17.5", sample_sd="0.7071067811865476", n_defined="2", paired_mean="0",
                     paired_sd="0", n_paired="2", n_coordinates="2",
                     n_eligible="6", n_undefined="0", overlap="1513"),
        }
        audit_aggregates(scores, summary)
        group = summary["hvg-2k", "disease", "direction-correct coverage"]
        group["paired_mean"] = "-0.3333333333333333"
        with self.assertRaises(ValueError):
            audit_aggregates(scores, summary)
        group["paired_mean"] = "3"
        scores["hvg-2k", 0, 17, "disease", "direction-correct coverage"]["reference_value"] = "18"
        with self.assertRaises(ValueError):
            audit_aggregates(scores, summary)
