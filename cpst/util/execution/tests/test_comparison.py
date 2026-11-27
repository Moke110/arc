import csv
import json
import os
import runpy
import sys
import unittest
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev
from tempfile import TemporaryDirectory
from unittest.mock import patch

from matplotlib.collections import PathCollection
from matplotlib.container import BarContainer, ErrorbarContainer

from cpst_util.comparison import METRIC_PANELS, render_comparison, write_paired_metrics

STAGE = Path(__file__).resolve().parents[3]
PANELS = (
    ("Δ Spearman", "ΔSpearman"), ("JAC_top200", "JAC-top200"),
    ("Δ PCC", "ΔPCC"), ("Δ RMSE", "ΔRMSE"), ("Δ R²", "ΔR²"),
    ("Δ E-distance", "ΔE-distance"), ("Δ MMD", "ΔMMD"),
)


class ComparisonTests(unittest.TestCase):
    def test_metric_roster_matches_adr_0109(self):
        self.assertEqual(METRIC_PANELS, PANELS)

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

    def test_published_metric_figures_show_every_point_and_full_sample_sd(self):
        receipts = sorted(STAGE.glob("*/log/figure-render.json"))
        receipts += sorted(STAGE.glob("*/*/log/figure-render.json"))
        self.assertEqual(len(receipts), 13)
        with TemporaryDirectory(dir=os.environ.get("TMPDIR")) as scratch:
            for receipt in receipts:
                with self.subTest(study=receipt.parent.parent.relative_to(STAGE)):
                    self.check_published_metrics(receipt, scratch)

    def check_published_metrics(self, receipt, scratch):
        step = receipt.parent.parent
        data = step / "figure/data"
        scripts = json.loads(receipt.read_text())["renderer_files_sha256"]
        prefix = step.relative_to(STAGE) / "figure/script"
        renderer, = (STAGE / path for path in scripts if Path(path).parent == prefix)
        checked = set()

        def check_saved_figure(figure, path, **_):
            if path.name == "learning-curve.png":
                return
            self.assertNotIn(path.name, checked)
            checked.add(path.name)
            with (data / path.with_suffix(".csv").name).open(newline="") as stream:
                scores = list(csv.DictReader(stream))
            self.assertEqual(len(figure.axes), 7)
            for position, ((metric, title), axis) in enumerate(zip(PANELS, figure.axes, strict=True)):
                self.assertEqual(axis.get_title(), title)
                self.assertEqual(axis.get_subplotspec().get_geometry(), (2, 4, position, position))
                self.check_metric_marks(axis, scores, metric)

        arguments = [str(renderer), "--data", str(data), "--output", scratch]
        with patch.object(sys, "argv", arguments), patch.object(sys, "path", list(sys.path)), \
             patch("cpst_util.comparison.save_figure", side_effect=check_saved_figure), \
             patch("cpst_util.plotting.save_figure", side_effect=check_saved_figure):
            runpy.run_path(str(renderer), run_name="__main__")
        self.assertEqual(checked, {"metrics.png", "metrics-paired.png"})

    def check_metric_marks(self, axis, scores, metric):
        expected = defaultdict(list)
        for row in scores:
            if row["metric"] == metric and row["value"]:
                expected[row["arm"]].append(float(row["value"]))

        dots = defaultdict(list)
        for artist in axis.collections:
            if isinstance(artist, PathCollection):
                for x, y in artist.get_offsets():
                    dots[round(x)].append(y)
        self.assertCountEqual([sorted(values) for values in dots.values()],
                              [sorted(values) for values in expected.values()])
        bars = {round(bar.get_x() + bar.get_width() / 2): bar.get_height()
                for container in axis.containers if isinstance(container, BarContainer)
                for bar in container.patches}
        errors = {round(segment[0, 0]): sorted(segment[:, 1])
                  for container in axis.containers if isinstance(container, ErrorbarContainer)
                  for lines in container.lines[2] for segment in lines.get_segments()}
        self.assertEqual(set(bars), set(dots))
        self.assertEqual(set(errors), {group for group, values in dots.items() if len(values) > 1})
        lower, upper = axis.get_ylim()
        self.assertLessEqual(lower, 0)
        self.assertGreaterEqual(upper, 0)
        for group, values in dots.items():
            center = mean(values)
            self.assertAlmostEqual(bars[group], center, delta=1e-12)
            self.assertLessEqual(lower, min(values))
            self.assertGreaterEqual(upper, max(values))
            if group in errors:
                deviation = stdev(values)
                bottom, top = errors[group]
                self.assertAlmostEqual(bottom, center - deviation, delta=1e-12)
                self.assertAlmostEqual(top, center + deviation, delta=1e-12)
                self.assertLessEqual(lower, bottom)
                self.assertGreaterEqual(upper, top)
