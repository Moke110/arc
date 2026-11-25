"""Exercise CSV-only metrics rendering at the command boundary."""

import csv
from itertools import product
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

RENDER = Path(__file__).resolve().parents[2] / "figure/script/render_metrics.py"
METRICS = ("JAC_top200", "JAC_pos100", "JAC_neg100", "Δ PCC", "Δ Spearman",
           "Δ RMSE", "Δ E-distance", "Δ MMD", "Δ R²")


class MetricsFigureTests(TestCase):
    def test_requires_all_scores_and_renders_the_complete_grid_from_only_csv(self):
        with TemporaryDirectory() as scratch:
            root = Path(scratch)
            data = root / "data"
            data.mkdir()
            rows = [(loss, alpha, metric, score / 108)
                    for score, (loss, alpha, metric) in enumerate(product(
                        ("mse", "cos", "pcc", "l1"), ("1/4", "1/3", "1/2"), METRICS))]
            command = [sys.executable, str(RENDER), "--data", str(data), "--output", str(root / "figures")]
            with (data / "delta-metrics.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("loss", "alpha", "metric", "value"))
                writer.writerows(rows[:-1])
            rejected = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("complete loss grid", rejected.stderr)

            with (data / "delta-metrics.csv").open("a", newline="") as stream:
                csv.writer(stream).writerow(rows[-1])
            subprocess.run(command, capture_output=True, text=True, check=True)
            self.assertEqual((root / "figures/delta-metrics.png").read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
