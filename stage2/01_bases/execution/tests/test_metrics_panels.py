"""Eight-panel tenv metrics figure: cross-run values, roles and CSV-only replay.

Acceptance for ADR 0063 / issue #162: indirect panels draw a cross-run mean line
with a min-max band, the other six panels draw a cross-run mean bar with one dot
per run, every colour comes from the stage palette, undefined entries are drawn
explicitly and never as zero, and a drawing call returns exactly the values it
drew so a renderer can rebuild the PNG from the recorded CSV alone.
"""
import math
import numpy as np
import sys
import tempfile
import unittest
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import palette
from bases import plots

FIXTURES = Path(__file__).resolve().parent / "fixtures"
FIXTURE_PNG = FIXTURES / "metrics_fixture.png"
FIXTURE_CSV = FIXTURES / "metrics_fixture.csv"
K = (1.0, 2.0, 3.0, 4.0, 5.0)


def rgba(value):
    """Normalise a matplotlib colour to a rounded RGBA tuple."""
    return tuple(round(float(channel), 4) for channel in value)


def canonical(value):
    """Make NaN and tuples comparable across a CSV round-trip."""
    if isinstance(value, dict):
        return {key: canonical(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return tuple(canonical(item) for item in value)
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return value


def run(label, direct, precision, spearman, max_k, recall, reasons=None):
    """One synthetic run sample for an entry."""
    return {
        "run": label,
        "direct": direct,
        "indirect": {"k": K, "precision": precision, "spearman": spearman,
                     "max_k": max_k, "recall": recall},
        **({"reasons": reasons} if reasons else {}),
    }


def synthetic_entries():
    """Two synthetic entrants covering defined and undefined values."""
    nan = float("nan")
    prnet = {
        "label": "PRnet", "hue": "blue",
        "runs": [
            run("fold_0-seed_17", {"PCC": 0.4, "Spearman": 0.35, "RMSE": 2.0, "R2": 0.1},
                (0.2, 0.4, 0.6, 0.8, 1.0), (nan, 0.2, 0.4, 0.6, 0.8), 37.0, 0.5),
            run("fold_0-seed_42", {"PCC": 0.5, "Spearman": 0.45, "RMSE": 1.8, "R2": 0.2},
                (0.4, 0.6, 0.8, 1.0, 1.0), (nan, 0.4, 0.6, 0.8, 1.0), 42.0, 0.6),
            run("fold_1-seed_17", {"PCC": 0.6, "Spearman": 0.55, "RMSE": 1.6, "R2": 0.3},
                (0.6, 0.8, 1.0, 1.0, 1.0), (nan, 0.6, 0.8, 1.0, 1.0), 51.0, 0.7),
            run("fold_1-seed_42", {"PCC": 0.7, "Spearman": 0.65, "RMSE": 1.4, "R2": 0.4},
                (0.8, 1.0, 1.0, 1.0, 1.0), (nan, 0.8, 1.0, 1.0, 1.0), 63.0, 0.8),
        ],
    }
    reference = {
        "label": "Mean", "hue": "gray",
        "runs": [
            run(f"fold_{fold}-seed_{seed}", {"PCC": nan, "Spearman": nan, "RMSE": 3.0, "R2": nan},
                (0.0, 0.0, 0.0, 0.0, 0.0), (nan, nan, nan, nan, nan), nan, nan,
                {"PCC": "constant_prediction", "Spearman": "constant_prediction",
                 "R2": "constant_prediction", "max_k@P80S40": "no_qualifying_k",
                 "Recall@max_k": "no_qualifying_k"})
            for fold in (0, 1) for seed in (17, 42)
        ],
    }
    return [prnet, reference]


class MetricsValueTests(unittest.TestCase):
    def setUp(self):
        self.entries = synthetic_entries()
        self.values = plots.metrics_values(self.entries)

    def panel(self, name):
        return next(panel for panel in self.values["panels"] if panel["panel"] == name)

    def entry(self, name, panel):
        return next(entry for entry in self.panel(panel)["entries"] if entry["series"] == name)

    def test_layout_is_two_rows_by_four_panels_in_adr_order(self):
        self.assertEqual(self.values["schema"], "stage2-tenv-metrics-values-v1")
        self.assertEqual(self.values["kind"], "metrics")
        self.assertEqual(self.values["layout"]["rows"], 2)
        self.assertEqual(self.values["layout"]["columns"], 4)
        self.assertEqual(self.values["layout"]["panels"],
                         ("PCC", "Spearman", "RMSE", "R2",
                          "Precision@k", "Spearman@k", "max_k@P80S40", "Recall@max_k"))
        self.assertEqual([panel["kind"] for panel in self.values["panels"]],
                         ["bar"] * 4 + ["curve"] * 2 + ["bar"] * 2)

    def test_configured_and_threshold_panels_are_bars_of_the_cross_run_mean(self):
        self.assertAlmostEqual(self.entry("PRnet", "PCC")["mean"], 0.55)
        self.assertEqual(self.entry("PRnet", "PCC")["values"], (0.4, 0.5, 0.6, 0.7))
        self.assertAlmostEqual(self.entry("PRnet", "RMSE")["mean"], 1.7)
        self.assertAlmostEqual(self.entry("PRnet", "max_k@P80S40")["mean"], 48.25)
        self.assertAlmostEqual(self.entry("PRnet", "Recall@max_k")["mean"], 0.65)
        self.assertEqual(self.entry("PRnet", "PCC")["run_labels"],
                         ("fold_0-seed_17", "fold_0-seed_42", "fold_1-seed_17", "fold_1-seed_42"))
        self.assertTrue(all(self.entry("PRnet", "PCC")["defined"]))

    def test_undefined_bar_is_nan_with_a_reason_and_never_zero(self):
        entry = self.entry("Mean", "PCC")
        self.assertTrue(math.isnan(entry["mean"]))
        self.assertFalse(entry["mean_defined"])
        self.assertEqual(entry["reason"], "constant_prediction")
        self.assertEqual(self.values["kind"], "metrics")
        self.assertFalse(any(entry["defined"]))
        self.assertAlmostEqual(self.entry("Mean", "RMSE")["mean"], 3.0)
        nodispatch = self.entry("Mean", "max_k@P80S40")
        self.assertTrue(math.isnan(nodispatch["mean"]))
        self.assertEqual(nodispatch["reason"], "no_qualifying_k")

    def test_indirect_panels_are_cross_run_mean_lines_with_min_max_bands(self):
        precision = self.entry("PRnet", "Precision@k")
        self.assertEqual(self.panel("Precision@k")["k"], K)
        self.assertEqual(precision["mean"], (0.5, 0.7, 0.85, 0.95, 1.0))
        self.assertEqual(precision["minimum"], (0.2, 0.4, 0.6, 0.8, 1.0))
        self.assertEqual(precision["maximum"], (0.8, 1.0, 1.0, 1.0, 1.0))
        self.assertEqual(precision["defined_count"], (4, 4, 4, 4, 4))
        spearman = self.entry("PRnet", "Spearman@k")
        self.assertTrue(math.isnan(spearman["mean"][0]))
        self.assertTrue(math.isnan(spearman["minimum"][0]))
        self.assertEqual(spearman["defined_count"][0], 0)
        self.assertEqual(spearman["mean"][1:], (0.5, 0.7, 0.85, 0.95))
        self.assertEqual(spearman["reason"], "undefined")

    def test_every_entry_records_the_hue_it_draws_with(self):
        self.assertEqual([entry["series"] for entry in self.values["panels"][0]["entries"]],
                         ["PRnet", "Mean"])
        self.assertEqual(self.entry("PRnet", "PCC")["hue"], "blue")
        self.assertEqual(self.entry("Mean", "PCC")["hue"], "gray")

    def test_duplicate_or_unknown_entries_are_rejected(self):
        entries = synthetic_entries()
        entries[1]["label"] = "PRnet"
        with self.assertRaisesRegex(ValueError, "unique"):
            plots.metrics_values(entries)
        entries = synthetic_entries()
        entries[0]["hue"] = "red"
        with self.assertRaisesRegex(ValueError, "palette"):
            plots.metrics_values(entries)
        entries = synthetic_entries()
        entries[0]["runs"][0]["indirect"]["k"] = (1.0, 2.0, 3.0)
        with self.assertRaisesRegex(ValueError, "k axis"):
            plots.metrics_values(entries)


class MetricsDrawingTests(unittest.TestCase):
    def setUp(self):
        self.values = plots.metrics_values(synthetic_entries())

    def figure(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        return plt.subplots(2, 4, figsize=(16, 7.5), layout="constrained")

    def test_draw_returns_exactly_the_values_it_drew(self):
        fig, axes = self.figure()
        try:
            drawn = plots.draw_metrics(self.values, axes)
            self.assertEqual(canonical(drawn), canonical(self.values))
            # Audit the actual artists against the returned values, so the test
            # fails if the drawing call ignores or alters the values it returns.
            for panel, ax in zip(self.values["panels"], axes.flat):
                if panel["kind"] == "bar":
                    patches = {round(float(patch.get_x() + patch.get_width() / 2), 6): patch
                               for patch in ax.patches}
                    for position, entry in enumerate(panel["entries"]):
                        if entry["mean_defined"]:
                            self.assertAlmostEqual(patches[float(position)].get_height(),
                                                   entry["mean"], places=12)
                        else:
                            self.assertNotIn(float(position), patches)
                    dots = [line for line in ax.lines if line.get_marker() == "o"]
                    for line, entry in zip(dots, [item for item in panel["entries"]
                                                  if item["mean_defined"]]):
                        expected = [value for value, defined
                                    in zip(entry["values"], entry["defined"]) if defined]
                        self.assertEqual(tuple(float(value) for value in line.get_ydata()),
                                         tuple(expected))
                else:
                    k = np.asarray(panel["k"], dtype=float)
                    for band_index, entry in enumerate(panel["entries"]):
                        color = palette.color(entry["hue"], palette.LINE_SHADE)
                        mean_lines = [line for line in ax.lines
                                      if line.get_color() == color
                                      and line.get_marker() in (None, "None", "")
                                      and len(line.get_xdata()) == len(k)]
                        self.assertEqual(len(mean_lines), 1, entry["series"])
                        np.testing.assert_array_equal(
                            np.asarray(mean_lines[0].get_xdata(), dtype=float), k)
                        np.testing.assert_array_equal(
                            np.asarray(mean_lines[0].get_ydata(), dtype=float),
                            np.asarray(entry["mean"], dtype=float))
                        if not any(entry["mean_defined"]):
                            continue
                        vertices = ax.collections[band_index].get_paths()[0].vertices
                        for index in range(len(k)):
                            if not entry["mean_defined"][index]:
                                continue
                            for bound in (entry["minimum"][index], entry["maximum"][index]):
                                matched = any(abs(float(vx) - k[index]) < 1e-9
                                              and abs(float(vy) - bound) < 1e-9
                                              for vx, vy in vertices)
                                self.assertTrue(matched, (panel["panel"], entry["series"], index, bound))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_bar_panels_use_the_entry_hue_edge_and_lighter_fill_with_one_dot_per_run(self):
        from matplotlib import colors as mcolors
        fig, axes = self.figure()
        try:
            plots.draw_metrics(self.values, axes)
            panel = axes[0, 2]  # RMSE: both entries are defined
            self.assertEqual(len(panel.patches), 2)
            self.assertEqual(rgba(panel.patches[0].get_edgecolor()),
                             rgba(mcolors.to_rgba(palette.color("blue", palette.LINE_SHADE))))
            self.assertEqual(rgba(panel.patches[0].get_facecolor()),
                             rgba(mcolors.to_rgba(palette.color("blue", palette.FILL_SHADE))))
            self.assertAlmostEqual(panel.patches[0].get_height(), 1.7)
            dots = [line for line in panel.lines if line.get_marker() == "o"][0]
            self.assertEqual(tuple(dots.get_xdata()), (0.0, 0.0, 0.0, 0.0))
            self.assertEqual(tuple(float(value) for value in dots.get_ydata()), (2.0, 1.8, 1.6, 1.4))
            self.assertEqual(dots.get_color(), palette.color("blue", palette.LINE_SHADE))
            self.assertEqual(dots.get_markeredgecolor(), palette.color("gray", 8))
            reference = [line for line in panel.lines if line.get_marker() == "o"][1]
            self.assertEqual(reference.get_color(), palette.color("gray", palette.LINE_SHADE))
            self.assertEqual(tuple(float(value) for value in reference.get_ydata()), (3.0, 3.0, 3.0, 3.0))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_undefined_entry_is_annotated_and_drawn_without_a_zero_bar(self):
        fig, axes = self.figure()
        try:
            plots.draw_metrics(self.values, axes)
            panel = axes[0, 0]
            self.assertEqual([round(patch.get_height(), 6) for patch in panel.patches], [0.55])
            self.assertTrue(any("undefined" in text.get_text() for text in panel.texts))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_indirect_panels_draw_a_band_a_mean_line_and_explicit_undefined_marks(self):
        from matplotlib import colors as mcolors
        fig, axes = self.figure()
        try:
            plots.draw_metrics(self.values, axes)
            panel = axes[1, 1]  # Spearman@k: undefined for k < 3
            self.assertEqual(len(panel.collections), 2)
            self.assertEqual(rgba(panel.collections[0].get_facecolor()[0]),
                             rgba(mcolors.to_rgba(palette.color("blue", palette.FILL_SHADE))))
            self.assertEqual(rgba(panel.collections[1].get_facecolor()[0]),
                             rgba(mcolors.to_rgba(palette.color("gray", palette.FILL_SHADE))))
            self.assertEqual(panel.lines[0].get_color(), palette.color("blue", palette.LINE_SHADE))
            undefined_marks = [line for line in panel.lines if line.get_marker() == "x"]
            self.assertTrue(undefined_marks)
            self.assertTrue(any("undefined" in text.get_text() for text in panel.texts))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)


class MetricsReplayTests(unittest.TestCase):
    def test_recorded_values_round_trip_and_replay_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        from PIL import Image
        values = plots.metrics_values(synthetic_entries())
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            values_csv = tmp / "values.csv"
            plots.write_plot_values(values, values_csv)
            self.assertEqual(canonical(plots.read_plot_values(values_csv)), canonical(values))
            first = tmp / "first.png"
            second = tmp / "second.png"
            plots.render_metrics(synthetic_entries(), first)
            plots.replay_metrics(values_csv, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with Image.open(first) as image:
                self.assertGreater(image.width, image.height)
                self.assertGreater(image.width, 1400)

    def test_committed_fixture_is_the_synthetic_values_and_replays_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        from PIL import Image
        self.assertTrue(FIXTURE_PNG.exists(), "run the fixture build command to create metrics_fixture.png")
        self.assertTrue(FIXTURE_CSV.exists(), "run the fixture build command to create metrics_fixture.csv")
        self.assertEqual(canonical(plots.read_plot_values(FIXTURE_CSV)),
                         canonical(plots.metrics_values(synthetic_entries())))
        with tempfile.TemporaryDirectory() as tmp:
            replay = Path(tmp) / "replay.png"
            plots.replay_metrics(FIXTURE_CSV, replay)
            self.assertEqual(replay.read_bytes(), FIXTURE_PNG.read_bytes())
        with Image.open(FIXTURE_PNG) as image:
            self.assertGreater(image.width, image.height)


if __name__ == "__main__":
    unittest.main()
