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
            rows = [(loss, alpha, metric, score - 20, 9 + METRICS.index(metric) - 20, score - 9 - METRICS.index(metric))
                    for score, (loss, alpha, metric) in enumerate(product(
                        ("mse", "cos", "pcc", "l1"), ("0.75", "1", "1.25"), METRICS))]
            command = [sys.executable, str(RENDER), "--data", str(data), "--output", str(root / "figures")]
            with (data / "delta-metrics.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("loss", "alpha", "metric", "value", "native_mse_score", "score_difference"))
                writer.writerows(rows[:-1])
            rejected = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("complete loss grid", rejected.stderr)

            with (data / "delta-metrics.csv").open("a", newline="") as stream:
                csv.writer(stream).writerow(rows[-1])
            with (data / "fold-scores.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("loss", "alpha", "metric", "fold", "seed", "value"))
                writer.writerows((*row[:3], fold, 42, row[3] + fold - 1)
                                 for row in rows for fold in range(3))
            subprocess.run(command, capture_output=True, text=True, check=True)
            scores = root / "figures/delta-metrics.png"
            differences = root / "figures/delta-metrics-improvement.png"
            for path in (scores, differences):
                self.assertEqual(path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertNotEqual(scores.read_bytes(), differences.read_bytes())
