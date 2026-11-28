"""Recheck the published report independently of the pandas aggregation."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_chem_report import audit
from chem_formal import FIGURE_RECEIPT, select_configuration, verify_figure_receipt
from chem_report import METRICS, paired_scores
from chem_runs import STAGE, STEP, axis_run, coordinates


class AcceptanceTest(unittest.TestCase):
    def test_sealed_figures_accept_only_the_verified_historical_renderer_transition(self):
        sealed = axis_run.read_json(STEP / FIGURE_RECEIPT)
        source = "util/output/api/cpst_util/comparison.py"
        renderers = axis_run.file_hashes(STAGE, sealed["renderer_files_sha256"])
        replay = {**{key: value for key, value in sealed.items() if key != "evidence"},
                  "renderer_files_sha256": renderers}

        verify_figure_receipt(sealed, replay)
        verify_figure_receipt(replay, replay)
        other = next(path for path in renderers if path != source)
        for changed in (
            {**replay, "png_files_sha256": {}},
            {**replay, "renderer_files_sha256": {**renderers, source: "0" * 64}},
            {**replay, "renderer_files_sha256": {**renderers, other: "0" * 64}},
        ):
            with self.assertRaisesRegex(ValueError, "figure receipt differs"):
                verify_figure_receipt(sealed, changed)
        wrong_sealed_sources = {**sealed["renderer_files_sha256"], source: "0" * 64}
        with self.assertRaisesRegex(ValueError, "figure receipt differs"):
            verify_figure_receipt({**sealed, "renderer_files_sha256": wrong_sealed_sources}, replay)

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
