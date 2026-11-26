"""Recheck the published report independently of the pandas aggregation."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_chem_report import audit
from chem_formal import select_configuration
from chem_report import METRICS, paired_scores
from chem_runs import coordinates


class AcceptanceTest(unittest.TestCase):
    def test_a_different_winner_on_each_primary_metric_keeps_kpgt(self):
        groups = [
            ("KPGT", "Δ PCC", 0., 0.), ("KPGT", "Δ RMSE", 0., 0.),
            ("E3FP", "Δ PCC", .020, .001), ("E3FP", "Δ RMSE", -.010, .001),
            ("UniMol", "Δ PCC", .010, .001), ("UniMol", "Δ RMSE", -.020, .001),
        ]
        summary = pd.DataFrame(groups, columns=["arm", "metric", "score_difference", "score_difference_sd"])
        selected = select_configuration(summary.assign(n_paired=9))
        self.assertEqual(selected["working_default"], "KPGT")
        self.assertEqual(selected["metric_winners"], {"Δ PCC": "E3FP", "Δ RMSE": "UniMol"})

    def test_changed_paired_summary_fails_independent_arithmetic_check(self):
        roster = coordinates(True)
        observed = pd.DataFrame([
            {"arm": arm, "fold": fold, "seed": seed, "metric": metric, "value": .5,
             "n_conditions": 4, "n_eligible": 4, "n_undefined": 0}
            for arm, fold, seed in roster for metric in METRICS
        ])
        paired, by_fold, summary = paired_scores(observed, roster)
        with TemporaryDirectory() as temporary:
            report = Path(temporary)
            for name, table in (("fold-scores", paired), ("paired-by-fold", by_fold),
                                ("paired-summary", summary)):
                table.to_csv(report / f"{name}.csv", index=False)
            self.assertEqual(audit(report, True)["score_rows"], 63)
            summary.loc[0, "score_difference"] = .25
            summary.to_csv(report / "paired-summary.csv", index=False)
            with self.assertRaisesRegex(ValueError, "aggregate"):
                audit(report, True)


if __name__ == "__main__":
    unittest.main()
