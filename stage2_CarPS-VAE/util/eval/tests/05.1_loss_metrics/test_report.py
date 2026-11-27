"""A published curve requires its observed epochs and selected checkpoint."""

from pathlib import Path
import sys
from unittest import TestCase
import pandas as pd

sys.path.insert(0, str((Path(__file__).resolve().parents[4] / "05_tsg_loss/05.1_loss_metrics/execution/tests" / Path(__file__).name).parents[1]))
from evaluate_loss import check_trace


class TraceTests(TestCase):
    def test_validation_minimum_binds_to_the_selected_checkpoint(self):
        trace = pd.DataFrame({"epoch": [1, 2, 3], "validation_loss": [3., 1., 2.],
                              "validation_mse": [1., 2., 3.]})
        trained = {"best_epoch": 2, "epochs": 3, "validation_loss": 1., "validation_mse": 2.}
        config = {"max_epochs": 3, "patience": 30}
        self.assertEqual(check_trace(trace, trained, config), 2)
        with self.assertRaisesRegex(ValueError, "selected checkpoint"):
            check_trace(trace, {**trained, "best_epoch": 1}, config)


class ScratchTests(TestCase):
    def test_verify_uses_requested_scratch_and_score_preserves_placement(self):
        from evaluate_loss import axis_run, report_scratch

        checkout = axis_run.REPO
        output = checkout / "stage2_CarPS-VAE/05_tsg_loss/05.1_loss_metrics/output/formal/report"
        scratch = checkout / "tmp/session"
        self.assertEqual(report_scratch("verify", output, scratch), scratch)
        self.assertEqual(report_scratch("verify", output, None), checkout / "tmp")
        self.assertEqual(report_scratch("score", output, None), output.parent)


class ScoreSummaryTests(TestCase):
    def test_score_differences_are_signed_and_folds_have_equal_weight(self):
        from evaluate_loss import summarize_loss_scores

        scores = pd.DataFrame([
            {"loss": loss, "alpha": alpha, "fold": fold, "metric": "Δ RMSE", "value": score,
             "n_conditions": count, "n_eligible": count, "n_undefined": 0}
            for loss, alpha, fold, score, count in (
                ("mse", "1", 0, 0.5, 10), ("mse", "1", 1, 0.75, 200),
                ("l1", "0.75", 0, 0.25, 10), ("l1", "0.75", 1, 0.5, 200),
            )
        ])
        rows = summarize_loss_scores(scores).set_index("loss")
        self.assertEqual(rows.loc["mse", "value"], 0.625)
        self.assertEqual(rows.loc["l1", "value"], 0.375)
        self.assertEqual(rows.loc["l1", "native_mse_score"], 0.625)
        self.assertEqual(rows.loc["mse", "score_difference"], 0.0)
        self.assertEqual(rows.loc["l1", "score_difference"], -0.25)
