"""Test the matched published-score aggregation used by the coverage audit."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "axis_local"))
from coverage_audit import nearest_similarity, summarize_scores


class CoverageAuditTest(unittest.TestCase):
    def test_training_neighbor_excludes_self_but_accepts_identical_distinct_parent(self):
        fingerprints = pd.DataFrame([[1, 1, 0], [1, 1, 0], [1, 0, 0], [1, 0, 0]],
                                    index=list("ABCD"))
        nearest = nearest_similarity(fingerprints, {"A", "B", "C"})
        pd.testing.assert_series_equal(nearest, pd.Series([1., 1., .5, 1.], index=list("ABCD")),
                                       check_dtype=False)

    def test_profile_and_parent_weights_apply_to_paired_baseline_contrasts(self):
        # Parent A has two wells and no gain. Parent B contributes a gain of 0.6.
        scores = pd.DataFrame({
            "drug_parent": ["A", "A", "B"],
            "Mean_cell": [0.1, 0.2, 0.3],
            "no-dose": [0.1, 0.2, 0.3],
            "latent-embedding": [0.1, 0.2, 0.9],
        })
        report, _ = summarize_scores(scores, 19, 25)
        points = report.set_index(["weight", "metric"])["value"]
        self.assertAlmostEqual(points["profile", "latent-embedding_minus_Mean_cell"], 0.2)
        self.assertAlmostEqual(points["parent", "latent-embedding_minus_Mean_cell"], 0.3)
        self.assertEqual(points["parent", "no-dose_minus_Mean_cell"], 0.0)

    def test_common_finite_rows_are_used_for_every_method(self):
        scores = pd.DataFrame({
            "drug_parent": ["A", "B"],
            "Mean_cell": [0.3, -0.4],
            "no-dose": [0.4, float("nan")],
            "latent-embedding": [0.5, 0.9],
        })
        report, _ = summarize_scores(scores, 19, 25)
        self.assertTrue((report["n_profiles"] == 1).all())
        self.assertTrue((report["n_parents"] == 1).all())
        self.assertTrue(report["ci_low"].isna().all())


if __name__ == "__main__":
    unittest.main()
