"""CSV-only figure contract: byte-identical PNG re-render in fresh isolated processes.

The fixture is the synthetic report tree from ``test_report``; the renderer is exercised
as a subprocess under ``python -I`` with an isolated ``MPLCONFIGDIR`` inside a temporary
tree that contains only the renderer and ``data/``.  No run directory, matrix, checkpoint
or truth product is reachable from that tree.
"""
import ast
import hashlib
import io
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
STEP = EXECUTION.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import render  # noqa: E402
import test_report as fixture  # noqa: E402

PYTHON = sys.executable
TENV_FIGURES = {"learning_curves_tenv.png", "metrics_tenv.png",
                "metrics_per_dose_tenv.png"}
TAHOE_FIGURES = {"learning_curves_tahoe_bulk.png", "metrics_tahoe_bulk.png",
                 "metrics_per_dose_tahoe_bulk.png"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_isolated(base, data, suite, name="run", output=None):
    """Render in a fresh process in a tree holding only the renderer, data/ and output/.

    ``common.py`` is the renderer's whole import closure beyond the standard library
    (the schedule constants); it imports no science stack and holds no payload.
    """
    isolated = Path(base) / name
    isolated.mkdir(parents=True, exist_ok=True)
    for module in ("render.py", "common.py"):
        shutil.copy2(EXECUTION / module, isolated / module)
    shutil.copytree(data, isolated / "data")
    environment = dict(os.environ, MPLCONFIGDIR=str(isolated / "mpl"))
    environment.pop("PYTHONPATH", None)
    destination = Path(output) if output is not None else isolated / "out"
    result = subprocess.run(
        [PYTHON, "-I", "render.py", "--data", "data", "--output", str(destination),
         "--suite", suite],
        cwd=isolated, env=environment, capture_output=True, text=True)
    return result, destination


def build_report(fixture_data, out):
    with redirect_stdout(io.StringIO()):
        report_result = fixture.report.build_report(
            runs_root=fixture_data["runs"], truth_root=fixture_data["truth"],
            repo_root=fixture_data["root"],
            mean_cell_source=fixture_data["mean_cell"],
            mean_cell_sha256=fixture_data["mean_cell_sha"],
            mean_cell_receipt=fixture_data["receipt"])
    fixture.report.write_report(report_result, Path(out) / "report", Path(out) / "figure")
    return report_result


class TestCsvOnlyRerender(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temporary = tempfile.TemporaryDirectory(prefix="tsg-view-render-")
        cls.tmp = Path(cls._temporary.name)
        cls.fixture = fixture.build_fixture(cls.tmp / "fixture")
        cls.out = cls.tmp / "report-out"
        build_report(cls.fixture, cls.out)
        cls.data = cls.out / "figure"

    @classmethod
    def tearDownClass(cls):
        cls._temporary.cleanup()

    def test_byte_identical_png_rerender_in_fresh_processes(self):
        first, first_out = run_isolated(self.tmp / "a", self.data, "all")
        self.assertEqual(first.returncode, 0, first.stderr)
        second, second_out = run_isolated(self.tmp / "b", self.data, "all")
        self.assertEqual(second.returncode, 0, second.stderr)
        expected = TENV_FIGURES | TAHOE_FIGURES | {"views_effect.png"}
        self.assertEqual({path.name for path in first_out.glob("*.png")}, expected)
        self.assertEqual({path.name for path in second_out.glob("*.png")}, expected)
        for name in sorted(expected):
            self.assertEqual((first_out / name).read_bytes(),
                             (second_out / name).read_bytes(), name)
        def digests(directory):
            return {path.name: sha256(path) for path in sorted(Path(directory).glob("*.png"))}
        self.assertEqual(digests(first_out), digests(second_out))
        status_a = json.loads((first_out / "render_status.json").read_text())
        status_b = json.loads((second_out / "render_status.json").read_text())
        self.assertEqual(status_a, status_b)
        self.assertEqual(status_a["suite"], "all")
        self.assertEqual(status_a["pending"], {})
        self.assertEqual(status_a["figures"], digests(first_out))

    def test_isolated_tree_holds_no_payload_and_renderer_imports_no_science_stack(self):
        first, output = run_isolated(self.tmp / "isolation", self.data, "tenv")
        self.assertEqual(first.returncode, 0, first.stderr)
        isolated = self.tmp / "isolation" / "run"
        payloads = [path for path in isolated.rglob("*")
                    if path.suffix in (".npy", ".pt", ".h5ad", ".pkl")]
        self.assertEqual(payloads, [])
        source = (EXECUTION / "render.py").read_text()
        imported = set()
        for module in ("render.py", "common.py"):
            tree = ast.parse((EXECUTION / module).read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0])
        self.assertLessEqual(imported, {"argparse", "csv", "hashlib", "json", "sys", "pathlib",
                                        "matplotlib", "__future__", "common"})
        for forbidden in (".npy", "checkpoint.pt", "torch", "pandas", "numpy", "output/truth"):
            self.assertNotIn(forbidden, source)
            self.assertNotIn(forbidden, (EXECUTION / "common.py").read_text())
        self.assertIn('metadata={"Software": None}', source)
        self.assertIn("matplotlib.use(\"Agg\")", source)

    def test_tenv_suite_without_tahoe_is_explicitly_pending(self):
        with tempfile.TemporaryDirectory(prefix="tsg-view-render-tenv-") as temporary:
            root = Path(temporary)
            tenv_only = fixture.build_fixture(root / "fixture", run_sources=("tenv",))
            out = root / "report-out"
            build_report(tenv_only, out)
            data = out / "figure"
            self.assertNotIn("tahoe_bulk", set(pd.read_csv(data / "metrics_bars.csv")["source"]))
            result, output = run_isolated(root / "iso", data, "tenv")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual({path.name for path in output.glob("*.png")}, TENV_FIGURES)
            status = json.loads((output / "render_status.json").read_text())
            self.assertIn("tahoe_bulk", status["pending"])
            self.assertIn("pending", status["pending"]["tahoe_bulk"])
            cross, _ = run_isolated(root / "cross", data, "cross", name="cross-run")
            self.assertNotEqual(cross.returncode, 0)
            self.assertIn("requires both sources", cross.stderr)
            every, _ = run_isolated(root / "all", data, "all", name="all-run")
            self.assertNotEqual(every.returncode, 0)
            self.assertIn("requires both sources", every.stderr)

    def test_cross_suite_leaves_accepted_tenv_panels_untouched(self):
        root = self.tmp / "two-pass"
        result, tenv_out = run_isolated(root, self.data, "tenv", name="tenv-pass")
        self.assertEqual(result.returncode, 0, result.stderr)
        before = {path.name: sha256(path) for path in tenv_out.glob("*.png")}
        self.assertEqual(set(before), TENV_FIGURES)
        cross, cross_out = run_isolated(root, self.data, "cross", name="cross-pass",
                                        output=tenv_out)
        self.assertEqual(cross.returncode, 0, cross.stderr)
        after = {path.name: sha256(path) for path in cross_out.glob("*.png")}
        self.assertEqual(set(after), TENV_FIGURES | {"views_effect.png"})
        for name, digest in before.items():
            self.assertEqual(after[name], digest, name)
        status = json.loads((cross_out / "render_status.json").read_text())
        self.assertEqual(status["suite"], "cross")
        self.assertEqual(set(status["figures"]), {"views_effect.png"})

    def test_missing_view_is_a_hard_error(self):
        with tempfile.TemporaryDirectory(prefix="tsg-view-render-missing-") as temporary:
            root = Path(temporary)
            data = root / "data"
            shutil.copytree(self.data, data)
            bars = pd.read_csv(data / "metrics_bars.csv")
            bars = bars[~((bars["source"] == "tenv") & (bars["view"] == "solo"))]
            bars.to_csv(data / "metrics_bars.csv", index=False, lineterminator="\n")
            result, _ = run_isolated(root, data, "tenv")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("grid mismatch", result.stderr)

    def test_documented_wrapper_runs_the_same_renderer(self):
        with tempfile.TemporaryDirectory(prefix="tsg-view-render-wrapper-") as temporary:
            output = Path(temporary) / "out"
            result = subprocess.run(
                [PYTHON, str(STEP / "figure" / "script" / "render.py"),
                 "--data", str(self.data), "--output", str(output), "--suite", "tenv"],
                capture_output=True, text=True,
                env=dict(os.environ, MPLCONFIGDIR=str(Path(temporary) / "mpl")))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual({path.name for path in output.glob("*.png")}, TENV_FIGURES)

    def test_renderer_validates_curve_grid_and_best_epochs(self):
        rows = render.read_rows(self.data, "learning_curves")
        curves = [row for row in rows if row["source"] == "tenv" and row["view"] == "full"
                  and row["fold"] == "0" and row["seed"] == "42"]
        self.assertEqual([int(row["epoch"]) for row in curves],
                         list(range(1, len(curves) + 1)))
        self.assertEqual(sum(int(row["best_validation"]) for row in curves), 1)
        with tempfile.TemporaryDirectory() as temporary:
            (Path(temporary) / "learning_curves.csv").write_text("bad,header\n1,2\n")
            with self.assertRaisesRegex(render.RenderError, "header does not match"):
                render.read_rows(temporary, "learning_curves")
        source = (EXECUTION / "render.py").read_text()
        for knob in ("DPI = 150", "facecolor=\"white\"", "metadata={\"Software\": None}",
                     "subplots_adjust", "figure.savefig"):
            self.assertIn(knob, source)

    def test_status_pins_inputs_and_figures(self):
        result, output = run_isolated(self.tmp / "status", self.data, "all")
        self.assertEqual(result.returncode, 0, result.stderr)
        status = json.loads((output / "render_status.json").read_text())
        self.assertEqual(status["schema"], "tsg-view-render-v1")
        self.assertEqual(status["renderer_version"], "tsg-view-render-v1")
        self.assertEqual(status["matplotlib"], render.matplotlib.__version__)
        self.assertEqual(status["dpi"], 150)
        self.assertIsNone(status["metadata_software"])
        for name, digest in status["figures"].items():
            self.assertEqual(sha256(output / name), digest)
        for name in ("metrics_bars.csv", "learning_curves.csv", "views_effect.csv"):
            self.assertEqual(status["csv_inputs"][name], sha256(self.data / name))
        self.assertIn("confounded by design", status["estimands"]["per_dose"])


if __name__ == "__main__":
    unittest.main()
