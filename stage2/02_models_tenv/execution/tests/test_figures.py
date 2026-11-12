"""CSV-only figure reproduction for the #169 report (ADR 0062 testing seam 4).

Each renderer rebuilds its PNG from the tracked ``figure/data`` plot values with
the step's ``output/`` unreadable, and the drawn values equal the recorded ones.
The three figures are the four direct-metric panels, the eight-panel ``@d``
reliability scan and the learning curve.
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

STEP_ROOT = Path(__file__).resolve().parents[1].parent
PYTHON = sys.executable
RENDERER = STEP_ROOT / "figure" / "script" / "render_figures.py"
DATA = STEP_ROOT / "figure" / "data"
PNGS = {"metrics": STEP_ROOT / "figure" / "metrics.png",
        "metrics_indirect": STEP_ROOT / "figure" / "metrics_indirect.png",
        "learning_curve": STEP_ROOT / "figure" / "learning_curve.png"}
PRODUCTS = {"metrics_values.csv": ("metrics", (1, 4)),
            "metrics_indirect_values.csv": ("metrics_indirect", (2, 4)),
            "learning_curve_values.csv": ("learning_curve", (1, 1))}


def render(data, out, cwd):
    completed = subprocess.run(
        [PYTHON, str(RENDERER), "--figure-data", str(data), "--out-dir", str(out)],
        cwd=str(cwd), capture_output=True, text=True, timeout=300)
    if completed.returncode != 0:
        raise AssertionError(f"renderer failed: {completed.stderr.strip()}")
    return completed


class ReplayTests(unittest.TestCase):
    def test_all_three_pngs_rebuild_byte_identically_from_recorded_values(self):
        self.assertTrue(RENDERER.is_file(), "figure/script/render_figures.py is missing")
        with tempfile.TemporaryDirectory() as tmp:
            render(DATA, tmp, tmp)
            for name, tracked in PNGS.items():
                rebuilt = Path(tmp) / tracked.name
                self.assertTrue(rebuilt.is_file(), name)
                self.assertEqual(rebuilt.read_bytes(), tracked.read_bytes(), name)

    def test_replay_succeeds_with_the_step_output_unreadable(self):
        output = STEP_ROOT / "output"
        hidden = STEP_ROOT / "output.test-hidden"
        moved = False
        if output.exists():
            os.replace(output, hidden)
            moved = True
        try:
            with tempfile.TemporaryDirectory() as tmp:
                render(DATA, tmp, tmp)
                for name, tracked in PNGS.items():
                    self.assertEqual((Path(tmp) / tracked.name).read_bytes(),
                                     tracked.read_bytes(), name)
        finally:
            if moved:
                os.replace(hidden, output)

    def test_replay_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "a", Path(tmp) / "b"
            render(DATA, first, tmp)
            render(DATA, second, tmp)
            for name in PNGS:
                self.assertEqual((first / PNGS[name].name).read_bytes(),
                                 (second / PNGS[name].name).read_bytes(), name)

    def test_recorded_values_are_the_shared_plot_schema(self):
        sys.path.insert(0, str(STEP_ROOT.parents[1] / "stage2/01_bases/output/api"))
        from bases import plots
        for name, (kind, shape) in PRODUCTS.items():
            recorded = plots.read_plot_values(DATA / name)
            if kind == "metrics":
                self.assertEqual(recorded["schema"], plots.METRICS_SCHEMA)
                self.assertEqual(recorded["kind"], "metrics")
                self.assertEqual(len(recorded["panels"]), 4)
                self.assertEqual([panel["kind"] for panel in recorded["panels"]],
                                 ["bar"] * 4)
            elif kind == "metrics_indirect":
                self.assertEqual(recorded["schema"], plots.INDIRECT_SCHEMA)
                self.assertEqual(recorded["kind"], "metrics_indirect")
                self.assertEqual(len(recorded["panels"]), 8)
                self.assertEqual([panel["kind"] for panel in recorded["panels"]],
                                 ["curve"] * 4 + ["bar", "bar", "heatmap", "bar"])
            else:
                self.assertEqual(recorded["schema"], plots.LEARNING_CURVE_SCHEMA)
                self.assertEqual(len(recorded["panels"]), 1)
            self.assertEqual((recorded["layout"]["rows"],
                              recorded["layout"]["columns"]), shape, name)

    def test_the_indirect_table_retires_the_at_k_family(self):
        text = (DATA / "metrics_indirect_values.csv").read_text()
        self.assertNotIn("@k", text)
        for retired in ("Precision@k", "Spearman@k", "max_k@P80S40", "Recall@max_k"):
            self.assertNotIn(retired, text)
        for required in ("Precision@d", "Recall@d", "Spearman@d", "Reliability@d",
                         "Reliability@d*", "d*", "Confusion@d*", "NullSplit@d*"):
            self.assertIn(required, text)


if __name__ == "__main__":
    unittest.main()
