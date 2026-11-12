"""The #165 matrix runner: declared plan, resource parsing, fail-closed publish."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

STEP_ROOT = Path(__file__).resolve().parents[1].parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_matrix  # noqa: E402

BASES = Path("/home/chang/arc/stage2/01_bases/output")


class PlanTests(unittest.TestCase):
    def test_plan_declares_the_two_by_two_matrix_in_order(self):
        cells = run_matrix.declared_identities(BASES / "fold_n_seed")
        self.assertEqual([(cell["fold"], cell["seed"]) for cell in cells],
                         [(0, 17), (0, 42), (1, 17), (1, 42)])
        self.assertEqual([cell["fold_row"] for cell in cells], [0, 0, 1, 1])
        self.assertEqual([cell["seed_row"] for cell in cells], [0, 1, 0, 1])
        self.assertEqual(run_matrix.ENTRIES, ("prnet", "mean"))

    def test_role_limits_are_the_full_fold_cohorts(self):
        import pandas as pd
        cohorts = pd.read_csv(BASES / "fold_n_seed/cohorts.csv", keep_default_na=False)
        limits = run_matrix.role_limits(cohorts, 0)
        self.assertEqual(limits, {"train": 2508, "validation": 297, "test": 1401})
        self.assertEqual(sum(limits.values()), len(cohorts[cohorts.fold.eq(0)]))


class TimeParsingTests(unittest.TestCase):
    SAMPLE = """
\tCommand being timed: "python"
\tUser time (seconds): 12.34
\tSystem time (seconds): 1.23
\tElapsed (wall clock) time (m:ss): 1:02.50
\tMaximum resident set size (kbytes): 1048576
"""

    def test_parses_wall_time_and_peak_rss(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.time"
            path.write_text(self.SAMPLE)
            value = run_matrix.parse_time(path)
        self.assertAlmostEqual(value["elapsed_seconds"], 62.5, places=3)
        self.assertEqual(value["peak_rss_kib"], 1048576)
        self.assertAlmostEqual(value["user_seconds"], 12.34)

    def test_missing_resource_fields_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.time"
            path.write_text("no measurement here")
            with self.assertRaisesRegex(run_matrix.MatrixError, "resident"):
                run_matrix.parse_time(path)


class FailClosedTests(unittest.TestCase):
    def test_only_a_complete_terminal_status_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "run"
            run.mkdir()
            for status in ("failed", "timeout", "over_budget", "interrupted"):
                (run / "terminal.json").write_text(json.dumps({"status": status}))
                with self.assertRaisesRegex(run_matrix.MatrixError, "not complete"):
                    run_matrix.check_terminal(run)
            (run / "terminal.json").write_text(json.dumps({"status": "complete"}))
            self.assertEqual(run_matrix.check_terminal(run)["status"], "complete")

    def test_publish_tree_replaces_a_stable_path_only_from_a_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            candidate = Path(tmp) / "candidate/inputs"
            candidate.mkdir(parents=True)
            (candidate / "manifest.json").write_text('{"ok": true}')
            stable = Path(tmp) / "output/inputs"
            stable.mkdir(parents=True)
            (stable / "old.txt").write_text("superseded")
            run_matrix.publish_tree(candidate, stable)
            self.assertEqual((stable / "manifest.json").read_text(), '{"ok": true}')
            self.assertFalse((stable / "old.txt").exists())
            self.assertFalse((Path(tmp) / "output/inputs.publishing").exists())
            self.assertFalse((Path(tmp) / "output/inputs.superseded").exists())

    def test_missing_candidate_product_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(run_matrix.MatrixError, "missing"):
                run_matrix.publish_tree(Path(tmp) / "absent", Path(tmp) / "output/inputs")


class FigureValueContractTests(unittest.TestCase):
    def values(self):
        sys.path.insert(0, str(run_matrix.repo_root() / "stage2/01_bases/output/api"))
        from bases import plots
        data = STEP_ROOT / "figure/data"
        return plots, plots.read_plot_values(data / "metrics_values.csv"), \
            plots.read_plot_values(data / "learning_curve_values.csv")

    def test_recorded_tables_follow_the_adr_layout(self):
        plots, metrics, curve = self.values()
        self.assertEqual(metrics["schema"], plots.METRICS_SCHEMA)
        self.assertEqual([panel["panel"] for panel in metrics["panels"]],
                         list(plots.PANEL_ORDER))
        self.assertEqual([panel["kind"] for panel in metrics["panels"]],
                         ["bar"] * 4 + ["curve"] * 2 + ["bar"] * 2)
        self.assertEqual(len(metrics["panels"][0]["entries"]), 2)
        self.assertEqual({entry["series"] for entry in metrics["panels"][0]["entries"]},
                         {"PRnet", "Mean"})
        self.assertEqual(curve["schema"], plots.LEARNING_CURVE_SCHEMA)
        self.assertEqual(len(curve["panels"]), 1)
        self.assertEqual(curve["panels"][0]["panel"], "PRnet")
        self.assertEqual(len(curve["panels"][0]["runs"]), 4)


if __name__ == "__main__":
    unittest.main()
