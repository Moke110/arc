import csv
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
STEP = REPO / "cpst/06_tsg_module/06.1_round1"
RENDERER = STEP / "execution/render_shape.py"
SHAPES = ("original", "widen_z-110", "widen_z-120", "widen_z-130",
          "widen_enc-2k", "widen_enc-4k", "widen_enc-6k", "deepen_enc", "deepen_dec", "deepen_both")
METRICS = ("JAC_top200", "JAC_pos100", "JAC_neg100", "Δ PCC", "Δ Spearman",
           "Δ RMSE", "Δ E-distance", "Δ MMD", "Δ R²")


class ShapeRendererTest(unittest.TestCase):
    def test_keeps_wide_encoder_spikes_visible_without_expanding_the_shared_scale(self):
        spec = spec_from_file_location("shape_renderer", RENDERER)
        renderer = module_from_spec(spec)
        spec.loader.exec_module(renderer)
        spikes = {"widen_enc-2k": 1e10, "widen_enc-4k": 1e12, "widen_enc-6k": 1e20}
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.write_table(directory / "epoch-objectives.csv",
                             ("shape", "fold", "seed", "epoch", "train_loss", "validation_loss"),
                             ((shape, fold, 42, epoch, spikes.get(shape, 2) if epoch == 1 else 1, 3 / epoch)
                              for shape, fold, epoch in product(SHAPES, range(3), (1, 2))))
            with patch.object(renderer.plotting, "save_figure") as save:
                renderer.render_curves(directory, directory)
            panels = {axis.get_title(): axis for axis in save.call_args.args[0].axes}
            common = panels["original"].get_ylim()
            self.assertLess(common[1], 10)
            for shape in (shape for shape in SHAPES if shape not in spikes):
                self.assertEqual(panels[shape].get_ylim(), common)
            for shape, peak in spikes.items():
                axis = panels[f"{shape} (separate y-scale)"]
                self.assertEqual(len(axis.lines), 6)
                self.assertGreaterEqual(axis.get_ylim()[1], peak)
                self.assertEqual(axis.get_yscale(), "log")

    def test_renders_and_replays_the_ten_shape_comparison_from_csvs(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.write_table(directory / "delta-metrics.csv", ("shape", "metric", "value", "score_difference"),
                             ((*key, .25, 0) for key in product(SHAPES, METRICS)))
            self.write_table(directory / "fold-scores.csv", ("shape", "metric", "fold", "value", "score_difference"),
                             ((*key, .25, 0) for key in product(SHAPES, METRICS, range(3))))
            curves = directory / "epoch-objectives.csv"
            self.write_table(curves, ("shape", "fold", "seed", "epoch", "train_loss", "validation_loss"),
                             ((shape, fold, 42, epoch, 2 / epoch, 3 / epoch)
                              for shape, fold, epoch in product(SHAPES, range(3), (1, 2))))
            for output in ("first", "repeat"):
                subprocess.run([sys.executable, RENDERER, "--data", directory,
                                "--output", directory / output], check=True, capture_output=True)
            expected = {"delta-metrics.png", "delta-metrics-paired-difference.png", "learning-curves.png"}
            self.assertEqual({path.name for path in (directory / "first").iterdir()}, expected)
            for name in expected:
                self.assertEqual((directory / "first" / name).read_bytes(),
                                 (directory / "repeat" / name).read_bytes())
            rows = curves.read_text().splitlines()
            curves.write_text("\n".join(row for row in rows if not row.startswith("deepen_both,")) + "\n")
            attempt = subprocess.run([sys.executable, RENDERER, "--data", directory,
                                      "--output", directory / "rejected"], capture_output=True, text=True)
            self.assertNotEqual(attempt.returncode, 0)
            self.assertIn("common three-fold run roster", attempt.stderr)

    @staticmethod
    def write_table(path, columns, rows):
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            writer.writerows(rows)
