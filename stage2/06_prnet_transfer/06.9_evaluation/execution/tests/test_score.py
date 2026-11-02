"""Score-driver tests: the test-phase gate, the plan, and the resume predicate."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import _support  # noqa: E402
import score  # noqa: E402


class ScoreGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "score")
        cls.pop = _support.write_synthetic_population(cls.root, drugs_per_fold=2, folds=(0,))

    def test_test_phase_requires_explicit_confirmation(self):
        code = score.main(
            [
                "--units", str(self.pop["units_path"]),
                "--phase", "test",
                "--out", str(self.root / "out"),
                "--plan",
            ]
        )
        self.assertEqual(code, 4)

    def test_plan_lists_units_and_scores_nothing(self):
        import contextlib
        import io

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = score.main(
                [
                    "--units", str(self.pop["units_path"]),
                    "--phase", "test",
                    "--out", str(self.root / "out"),
                    "--plan",
                    "--confirm-frozen-selection",
                ]
            )
        self.assertEqual(code, 0)
        self.assertEqual(buffer.getvalue().count("e9_f"), 6)

    def test_subset_limit_is_honoured_in_plan(self):
        import contextlib
        import io

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = score.main(
                [
                    "--units", str(self.pop["units_path"]),
                    "--phase", "test",
                    "--out", str(self.root / "out"),
                    "--plan",
                    "--confirm-frozen-selection",
                    "--limit", "4",
                ]
            )
        self.assertEqual(code, 0)
        self.assertEqual(buffer.getvalue().count("e9_f"), 4)

    def test_unknown_unit_id_is_rejected(self):
        code = score.main(
            [
                "--units", str(self.pop["units_path"]),
                "--phase", "test",
                "--out", str(self.root / "out"),
                "--plan",
                "--confirm-frozen-selection",
                "--only-unit", "does-not-exist",
            ]
        )
        self.assertEqual(code, 2)


class CumulativeRunRecordTests(unittest.TestCase):
    """A subset re-run must not drop the earlier units from the run record."""

    def setUp(self) -> None:
        self.root = _support.clean_dir(_support.tmp_root() / "cumulative")

    def _record(self, unit_id, status="completed", wall=1.0, peak=100.0):
        return {
            "unit_id": unit_id, "phase": "test", "status": status,
            "fold": 0, "seed": 20261021, "init": "random", "lr": 0.001,
            "created_utc": "2026-01-02T00:00:00Z",
            "steps": {
                "predict": {"exit_code": 0, "wall_s": wall, "peak_rss_mb": peak,
                            "started_utc": "2026-01-02T00:00:00Z"},
                "evaluate": {"exit_code": 0, "wall_s": wall, "peak_rss_mb": peak},
            },
            "gpu": {"peak_mb": 2500},
            "verification": {"eval_manifest_sha256": "a" * 64},
        }

    def test_subset_rerun_merges_with_earlier_rows(self):
        out = self.root / "merge_out"
        rows = score._write_runs_csv(
            out, [self._record("u1"), self._record("u2", wall=2.0)]
        )
        self.assertEqual(len(rows), 2)
        rows = score._write_runs_csv(out, [self._record("u3", wall=3.0)])
        self.assertEqual([r["unit_id"] for r in rows], ["u1", "u2", "u3"])
        rows = score._write_runs_csv(out, [self._record("u2", wall=9.0)])
        self.assertEqual(len(rows), 3)
        self.assertAlmostEqual(float(rows[1]["predict_wall_s"]), 9.0)

    def test_cumulative_summary_counts_all_rows(self):
        rows = [
            {"unit_id": "u1", "status": "completed", "predict_wall_s": "1.0",
             "evaluate_wall_s": "1.0", "predict_peak_rss_mb": "100",
             "evaluate_peak_rss_mb": "120", "gpu_peak_mb": "2000",
             "started_utc": "2026-01-02T00:00:00Z"},
            {"unit_id": "u2", "status": "failed", "predict_wall_s": "0.5",
             "evaluate_wall_s": "", "predict_peak_rss_mb": "90",
             "evaluate_peak_rss_mb": "", "gpu_peak_mb": "2100",
             "started_utc": "2026-01-02T00:01:00Z", "error": "boom"},
        ]
        summary = score._summary_from_rows(
            rows, phase="test", preflight={"mem_available_gib": 40.0},
            units_manifest=Path("units.json"),
            units_manifest_sha256="f" * 64, selection_complete=True,
        )
        self.assertEqual(summary["units_planned"], 2)
        self.assertEqual(summary["units_completed"], 1)
        self.assertEqual(summary["units_failed"], 1)
        self.assertEqual(summary["failures"][0]["unit_id"], "u2")
        self.assertAlmostEqual(summary["wall_s"], 2.5)
        self.assertEqual(summary["host_peak_rss_mb"], 120.0)
        self.assertEqual(summary["gpu_peak_mb"], 2100)
        self.assertEqual(summary["started_utc"], "2026-01-02T00:00:00Z")


if __name__ == "__main__":
    unittest.main()
