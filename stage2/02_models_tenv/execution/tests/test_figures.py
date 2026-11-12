"""CSV-only figure reproduction for the #165 report (ADR 0062 testing seam 4).

Each renderer rebuilds its PNG from the tracked ``figure/data`` plot values with
the step's ``output/`` unreadable, and the drawn values equal the recorded ones.
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
        "learning_curve": STEP_ROOT / "figure" / "learning_curve.png"}


def render(data, out, cwd):
    completed = subprocess.run(
        [PYTHON, str(RENDERER), "--figure-data", str(data), "--out-dir", str(out)],
        cwd=str(cwd), capture_output=True, text=True, timeout=300)
    if completed.returncode != 0:
        raise AssertionError(f"renderer failed: {completed.stderr.strip()}")
    return completed


class ReplayTests(unittest.TestCase):
    def test_both_pngs_rebuild_byte_identically_from_recorded_values(self):
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
        for name, expected in (("metrics_values.csv", plots.METRICS_SCHEMA),
                               ("learning_curve_values.csv", plots.LEARNING_CURVE_SCHEMA)):
            frame = plots.read_plot_values(DATA / name)
            self.assertEqual(frame["schema"], expected)
            self.assertTrue(frame["panels"])


if __name__ == "__main__":
    unittest.main()
