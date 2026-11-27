import csv
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from stage2_util.comparison import METRIC_PANELS, render_comparison, write_paired_metrics


class ComparisonTests(unittest.TestCase):
    def test_csv_replay_pairs_only_matching_fold_and_seed(self):
        with TemporaryDirectory(dir=os.environ.get("TMPDIR")) as scratch:
            data = Path(scratch)
            arms = {"base": ("Baseline", 20, 2), "model": ("Model", 220, 2)}
            fits = [("base", 0, 42, None), ("base", 1, 42, 0.2),
                    ("model", 0, 42, 0.4), ("model", 1, 42, 0.5), ("model", 1, 17, 0.9)]
            with (data / "metrics.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("arm", "fold", "seed", "metric", "value"))
                writer.writerows((arm, fold, seed, metric, value)
                                 for arm, fold, seed, value in fits for metric, _ in METRIC_PANELS)
            with (data / "learning-curve.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("arm", "fold", "seed", "epoch", "train", "validation"))
                writer.writerows((arm, fold, seed, epoch, 0.25 / epoch, 0.5 / epoch)
                                 for arm, fold, seed, _ in fits for epoch in (1, 2))
            write_paired_metrics(data, "base")
            with (data / "metrics-paired.csv").open(newline="") as stream:
                rows = [row for row in csv.DictReader(stream) if row["metric"] == "Δ PCC"]
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["value"], "")
            self.assertAlmostEqual(float(rows[1]["value"]), 0.3)
            self.assertEqual(rows[2]["value"], "")
            for subdir in ("first", "second"):
                render_comparison(data, data / subdir, arms=arms,
                                  reference="base", title="Fixture")
            self.assertEqual({path.name for path in (data / "first").iterdir()},
                             {"learning-curve.png", "metrics.png", "metrics-paired.png"})
            for figure in (data / "first").iterdir():
                self.assertEqual(figure.read_bytes(), (data / "second" / figure.name).read_bytes())
