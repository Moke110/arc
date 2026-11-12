"""ADR 0064 comparison figures: direct 1x4 and reliability 2x4 plotted values.

Acceptance for #169: ``metrics_values`` records the four direct-metric bars in
one row by four columns; ``indirect_metrics_values`` records the four ``@d``
curves over ``d'`` with the cross-run mean line and min-max band, each run's
``d*`` marker and the ``F1@d`` reference line, then the ``Reliability@d*`` and
``d*`` bars, the macro three-class confusion heat map and the null/changed
split.  Every colour comes from the stage palette, every recorded layout drives
the figure, and a drawing call returns exactly the values it drew so a renderer
can rebuild the PNG from the recorded CSV alone.
"""
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import palette
from bases import plots

FIXTURES = Path(__file__).resolve().parent / "fixtures"
FIXTURE_PNG = FIXTURES / "metrics_fixture.png"
FIXTURE_CSV = FIXTURES / "metrics_fixture.csv"
D_PRIME = (0.01, 0.25, 0.5, 1.0)
CLASSES = ("up", "not_significant", "down")
DIRECT = ("PCC", "Spearman", "RMSE", "R2")


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


def direct_run(label, direct, reasons=None):
    """One synthetic run sample for the direct-metric figure."""
    return {"run": label, "direct": direct,
            **({"reasons": reasons} if reasons else {})}


def synthetic_direct_entries():
    """Two synthetic entrants: PRnet defined everywhere, Mean undefined on PCC/R2."""
    nan = float("nan")
    prnet = {
        "label": "PRnet", "hue": "blue",
        "runs": [
            direct_run("fold_0-seed_17", {"PCC": 0.4, "Spearman": 0.35, "RMSE": 2.0, "R2": 0.1}),
            direct_run("fold_0-seed_42", {"PCC": 0.5, "Spearman": 0.45, "RMSE": 1.8, "R2": 0.2}),
            direct_run("fold_1-seed_17", {"PCC": 0.6, "Spearman": 0.55, "RMSE": 1.6, "R2": 0.3}),
            direct_run("fold_1-seed_42", {"PCC": 0.7, "Spearman": 0.65, "RMSE": 1.4, "R2": 0.4}),
        ],
    }
    reference = {
        "label": "Mean", "hue": "gray",
        "runs": [
            direct_run(f"fold_{fold}-seed_{seed}",
                       {"PCC": nan, "Spearman": nan, "RMSE": 3.0 + 0.1 * index, "R2": nan},
                       {"PCC": "constant_prediction", "Spearman": "constant_prediction",
                        "R2": "constant_prediction"})
            for index, (fold, seed) in enumerate(((0, 17), (0, 42), (1, 17), (1, 42)))
        ],
    }
    return [prnet, reference]


def confusion(scale=1):
    return {"classes": list(CLASSES),
            "matrix": [[1 * scale, 2 * scale, 3 * scale],
                       [4 * scale, 5 * scale, 6 * scale],
                       [7 * scale, 8 * scale, 9 * scale]]}


def indirect_run(label, *, precision, recall, spearman, f1, reliability,
                 d_star=0.5, reliability_at_d_star=0.5,
                 null=0.25, changed=0.75, scale=1):
    """One synthetic run sample for the reliability figure."""
    return {"run": label,
            "d_prime": list(D_PRIME),
            "precision": list(precision), "recall": list(recall),
            "spearman": list(spearman), "f1": list(f1),
            "reliability": list(reliability),
            "d_star_prime": d_star,
            "reliability_at_d_star": reliability_at_d_star,
            "null_reliability_at_d_star": null,
            "changed_reliability_at_d_star": changed,
            "confusion": confusion(scale)}


def synthetic_indirect_entries():
    prnet = {
        "label": "PRnet", "hue": "blue",
        "runs": [
            indirect_run("fold_0-seed_17",
                         precision=(0.1, 0.2, 0.3, 0.4), recall=(0.4, 0.5, 0.6, 0.7),
                         spearman=(0.0, 0.2, 0.4, 0.6), f1=(0.2, 0.3, 0.4, 0.5),
                         reliability=(0.2, 0.3, 0.4, 0.5), d_star=0.5,
                         reliability_at_d_star=0.3, null=0.2, changed=0.6, scale=1),
            indirect_run("fold_0-seed_42",
                         precision=(0.3, 0.4, 0.5, 0.6), recall=(0.6, 0.7, 0.8, 0.9),
                         spearman=(0.1, 0.3, 0.5, 0.7), f1=(0.4, 0.5, 0.6, 0.7),
                         reliability=(0.4, 0.5, 0.6, 0.7), d_star=0.75,
                         reliability_at_d_star=0.6, null=0.3, changed=0.7, scale=2),
        ],
    }
    reference = {
        "label": "Mean", "hue": "orange",
        "runs": [
            indirect_run("fold_0-seed_17",
                         precision=(0.01, 0.02, 0.03, 0.04), recall=(0.04, 0.05, 0.06, 0.07),
                         spearman=(0.0, 0.0, 0.0, 0.0), f1=(0.02, 0.03, 0.04, 0.05),
                         reliability=(0.02, 0.03, 0.04, 0.05), d_star=0.25,
                         reliability_at_d_star=0.03, null=0.0, changed=0.01, scale=3),
        ],
    }
    return [prnet, reference]


class DirectMetricsValueTests(unittest.TestCase):
    def setUp(self):
        self.entries = synthetic_direct_entries()
        self.values = plots.metrics_values(self.entries)

    def panel(self, name):
        return next(panel for panel in self.values["panels"] if panel["panel"] == name)

    def entry(self, name, panel):
        return next(entry for entry in self.panel(panel)["entries"] if entry["series"] == name)

    def test_layout_is_one_row_by_four_direct_panels(self):
        self.assertEqual(self.values["schema"], "stage2-tenv-metrics-values-v2")
        self.assertEqual(self.values["kind"], "metrics")
        self.assertEqual(self.values["layout"],
                         {"rows": 1, "columns": 4, "panels": DIRECT})
        self.assertEqual([panel["panel"] for panel in self.values["panels"]], list(DIRECT))
        self.assertEqual([panel["kind"] for panel in self.values["panels"]], ["bar"] * 4)
        self.assertEqual([panel["title"] for panel in self.values["panels"]], list(DIRECT))

    def test_panels_are_bars_of_the_cross_run_mean(self):
        self.assertAlmostEqual(self.entry("PRnet", "PCC")["mean"], 0.55)
        self.assertEqual(self.entry("PRnet", "PCC")["values"], (0.4, 0.5, 0.6, 0.7))
        self.assertAlmostEqual(self.entry("PRnet", "RMSE")["mean"], 1.7)
        self.assertEqual(self.entry("PRnet", "PCC")["run_labels"],
                         ("fold_0-seed_17", "fold_0-seed_42",
                          "fold_1-seed_17", "fold_1-seed_42"))
        self.assertTrue(all(self.entry("PRnet", "PCC")["defined"]))

    def test_undefined_bar_is_nan_with_a_reason_and_never_zero(self):
        entry = self.entry("Mean", "PCC")
        self.assertTrue(math.isnan(entry["mean"]))
        self.assertFalse(entry["mean_defined"])
        self.assertEqual(entry["reason"], "constant_prediction")
        self.assertFalse(any(entry["defined"]))
        self.assertAlmostEqual(self.entry("Mean", "RMSE")["mean"], 3.15)
        self.assertIsNone(self.entry("Mean", "RMSE")["reason"])

    def test_every_entry_records_the_hue_it_draws_with(self):
        self.assertEqual([entry["series"] for entry in self.values["panels"][0]["entries"]],
                         ["PRnet", "Mean"])
        self.assertEqual(self.entry("PRnet", "PCC")["hue"], "blue")
        self.assertEqual(self.entry("Mean", "PCC")["hue"], "gray")

    def test_duplicate_or_unknown_entries_are_rejected(self):
        entries = synthetic_direct_entries()
        entries[1]["label"] = "PRnet"
        with self.assertRaisesRegex(ValueError, "unique"):
            plots.metrics_values(entries)
        entries = synthetic_direct_entries()
        entries[0]["hue"] = "red"
        with self.assertRaisesRegex(ValueError, "palette"):
            plots.metrics_values(entries)
        with self.assertRaisesRegex(ValueError, "unknown metric panel"):
            plots.metrics_values(synthetic_direct_entries(),
                                 panel_titles={"Precision@k": "retired"})


class IndirectMetricsValueTests(unittest.TestCase):
    def setUp(self):
        self.entries = synthetic_indirect_entries()
        self.values = plots.indirect_metrics_values(self.entries)

    def panel(self, name):
        return next(panel for panel in self.values["panels"] if panel["panel"] == name)

    def entry(self, name, panel):
        return next(entry for entry in self.panel(panel)["entries"] if entry["series"] == name)

    def test_layout_is_two_rows_by_four_in_adr_order(self):
        self.assertEqual(self.values["schema"], "stage2-tenv-metrics-indirect-values-v1")
        self.assertEqual(self.values["kind"], "metrics_indirect")
        self.assertEqual(self.values["layout"], {
            "rows": 2, "columns": 4,
            "panels": ("Precision@d", "Recall@d", "Spearman@d", "Reliability@d",
                       "Reliability@d*", "d*", "Confusion@d*", "NullSplit@d*")})
        self.assertEqual([panel["kind"] for panel in self.values["panels"]],
                         ["curve"] * 4 + ["bar", "bar", "heatmap", "bar"])

    def test_curves_are_cross_run_mean_lines_with_min_max_bands(self):
        precision = self.entry("PRnet", "Precision@d")
        self.assertEqual(self.panel("Precision@d")["d_prime"], D_PRIME)
        for actual, expected in ((precision["mean"], (0.2, 0.3, 0.4, 0.5)),
                                 (precision["minimum"], (0.1, 0.2, 0.3, 0.4)),
                                 (precision["maximum"], (0.3, 0.4, 0.5, 0.6)),
                                 (self.entry("Mean", "Reliability@d")["mean"],
                                  (0.02, 0.03, 0.04, 0.05))):
            for got, want in zip(actual, expected):
                self.assertAlmostEqual(got, want)
        self.assertEqual(precision["defined_count"], (2, 2, 2, 2))
        self.assertEqual(self.panel("Precision@d")["x_label"], "d'")

    def test_reliability_panel_records_each_run_d_star_and_the_f1_reference(self):
        reliability = self.entry("PRnet", "Reliability@d")
        for got, want in zip(reliability["reference"], (0.3, 0.4, 0.5, 0.6)):
            self.assertAlmostEqual(got, want)
        self.assertEqual([(marker["run"], marker["x"], marker["y"])
                          for marker in reliability["d_star"]],
                         [("fold_0-seed_17", 0.5, 0.3),
                          ("fold_0-seed_42", 0.75, 0.6)])
        self.assertNotIn("reference", self.entry("PRnet", "Precision@d"))

    def test_summary_bars_carry_the_run_mean_and_each_run_value(self):
        reliability = self.entry("PRnet", "Reliability@d*")
        self.assertEqual(reliability["values"], (0.3, 0.6))
        self.assertAlmostEqual(reliability["mean"], 0.45)
        self.assertEqual(self.entry("Mean", "d*")["values"], (0.25,))
        self.assertEqual(self.panel("d*")["y_label"], "d*")

    def test_confusion_panel_sums_each_entrys_per_run_counts(self):
        panel = self.panel("Confusion@d*")
        self.assertEqual(panel["classes"], CLASSES)
        prnet = self.entry("PRnet", "Confusion@d*")
        self.assertEqual(prnet["matrix"],
                         ((3, 6, 9), (12, 15, 18), (21, 24, 27)))
        self.assertEqual(self.entry("Mean", "Confusion@d*")["matrix"],
                         ((3, 6, 9), (12, 15, 18), (21, 24, 27)))
        self.assertEqual(panel["y_label"], "predicted class")
        self.assertEqual(panel["x_label"], "measured class")

    def test_null_and_changed_split_has_one_series_per_class(self):
        panel = self.panel("NullSplit@d*")
        self.assertEqual([entry["series"] for entry in panel["entries"]],
                         ["PRnet | null", "PRnet | changed",
                          "Mean | null", "Mean | changed"])
        self.assertAlmostEqual(self.entry("PRnet | null", "NullSplit@d*")["mean"], 0.25)
        self.assertAlmostEqual(self.entry("PRnet | changed", "NullSplit@d*")["mean"], 0.65)

    def test_d_prime_axis_must_be_increasing_bounded_and_shared(self):
        for axis, message in (((0.25, 0.25, 0.5, 1.0), "strictly increase"),
                              ((0.0, 0.25, 0.5, 1.0), "0 < d' <= 1"),
                              ((0.25, 0.5, 0.75, 1.5), "0 < d' <= 1")):
            entries = synthetic_indirect_entries()
            for run in entries[0]["runs"]:
                run["d_prime"] = list(axis)
            with self.assertRaisesRegex(ValueError, message):
                plots.indirect_metrics_values(entries)
        entries = synthetic_indirect_entries()
        entries[1]["runs"][0]["d_prime"] = [0.1, 0.2, 0.3, 0.4]
        with self.assertRaisesRegex(ValueError, "one d_prime axis"):
            plots.indirect_metrics_values(entries)
        entries = synthetic_indirect_entries()
        entries[0]["runs"][0]["precision"] = (0.1, 0.2)
        with self.assertRaisesRegex(ValueError, "d' axis"):
            plots.indirect_metrics_values(entries)

    def test_confusion_classes_must_be_the_three_declared_classes(self):
        entries = synthetic_indirect_entries()
        entries[0]["runs"][0]["confusion"] = {
            "classes": ["up", "down", "flat"],
            "matrix": [[0, 1, 2], [3, 4, 5], [6, 7, 8]]}
        with self.assertRaisesRegex(ValueError, "confusion classes"):
            plots.indirect_metrics_values(entries)


class MetricsDrawingTests(unittest.TestCase):
    def figure(self, rows, columns):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        return plt.subplots(rows, columns, figsize=(4.2 * columns, 4.2 * rows),
                            layout="constrained")

    def close(self, fig):
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_direct_draw_returns_exactly_the_values_it_drew(self):
        values = plots.metrics_values(synthetic_direct_entries())
        fig, axes = self.figure(1, 4)
        try:
            drawn = plots.draw_metrics(values, axes)
            self.assertEqual(canonical(drawn), canonical(values))
            for panel, ax in zip(values["panels"], np.asarray(axes).flat):
                patches = {round(float(patch.get_x() + patch.get_width() / 2), 6): patch
                           for patch in ax.patches}
                for position, entry in enumerate(panel["entries"]):
                    if entry["mean_defined"]:
                        self.assertAlmostEqual(patches[float(position)].get_height(),
                                               entry["mean"], places=12)
                    else:
                        self.assertNotIn(float(position), patches)
        finally:
            self.close(fig)

    def test_direct_bar_panels_use_the_entry_hue_and_one_dot_per_run(self):
        from matplotlib import colors as mcolors
        values = plots.metrics_values(synthetic_direct_entries())
        fig, axes = self.figure(1, 4)
        try:
            plots.draw_metrics(values, axes)
            panel = np.asarray(axes).flat[2]  # RMSE: both entries defined
            self.assertEqual(len(panel.patches), 2)
            self.assertEqual(rgba(panel.patches[0].get_edgecolor()),
                             rgba(mcolors.to_rgba(palette.color("blue", palette.LINE_SHADE))))
            self.assertEqual(rgba(panel.patches[0].get_facecolor()),
                             rgba(mcolors.to_rgba(palette.color("blue", palette.FILL_SHADE))))
            dots = [line for line in panel.lines if line.get_marker() == "o"][0]
            self.assertEqual(tuple(float(value) for value in dots.get_ydata()),
                             (2.0, 1.8, 1.6, 1.4))
            self.assertEqual(dots.get_markeredgecolor(), palette.color("gray", 8))
        finally:
            self.close(fig)

    def test_direct_undefined_entry_is_annotated_and_drawn_without_a_zero_bar(self):
        values = plots.metrics_values(synthetic_direct_entries())
        fig, axes = self.figure(1, 4)
        try:
            plots.draw_metrics(values, axes)
            panel = np.asarray(axes).flat[0]
            self.assertEqual([round(patch.get_height(), 6) for patch in panel.patches], [0.55])
            self.assertTrue(any("undefined" in text.get_text() for text in panel.texts))
        finally:
            self.close(fig)

    def test_indirect_curves_draw_a_band_a_mean_line_a_reference_and_markers(self):
        from matplotlib import colors as mcolors
        values = plots.indirect_metrics_values(synthetic_indirect_entries())
        fig, axes = self.figure(2, 4)
        try:
            plots.draw_metrics(values, axes)
            panel = axes[0, 3]  # Reliability@d
            self.assertEqual(len(panel.collections), 2)
            self.assertEqual(rgba(panel.collections[0].get_facecolor()[0]),
                             rgba(mcolors.to_rgba(palette.color("blue", palette.FILL_SHADE))))
            blue_mean = [line for line in panel.lines
                         if line.get_color() == palette.color("blue", palette.LINE_SHADE)
                         and line.get_marker() in (None, "None", "")]
            self.assertTrue(blue_mean)
            np.testing.assert_allclose(
                np.asarray(blue_mean[0].get_ydata(), dtype=float),
                np.asarray((0.3, 0.4, 0.5, 0.6), dtype=float))
            reference = [line for line in panel.lines if line.get_linestyle() == "--"]
            self.assertTrue(reference)
            self.assertEqual(reference[0].get_color(),
                             palette.color("blue", palette.REFERENCE_SHADE))
            self.assertLessEqual(reference[0].get_linewidth(), 1.0)
            markers = [line for line in panel.lines
                       if line.get_marker() == "o"
                       and tuple(float(value) for value in line.get_xdata()) == (0.5,)]
            self.assertTrue(markers)
            self.assertEqual(markers[0].get_color(),
                             palette.color("blue", palette.MARKER_SHADE))
            self.assertEqual(tuple(float(value) for value in markers[0].get_ydata()), (0.3,))
        finally:
            self.close(fig)

    def test_heatmap_draws_each_entrys_recorded_count_matrix(self):
        values = plots.indirect_metrics_values(synthetic_indirect_entries())
        fig, axes = self.figure(2, 4)
        try:
            plots.draw_metrics(values, axes)
            panel = axes[1, 2]
            self.assertEqual(len(panel.images), 2)
            expected = [entry["matrix"] for entry in values["panels"][6]["entries"]]
            for image, matrix in zip(panel.images, expected):
                np.testing.assert_array_equal(np.asarray(image.get_array(), dtype=float),
                                              np.asarray(matrix, dtype=float))
            self.assertTrue(any(text.get_text() == "PRnet" for text in panel.texts))
        finally:
            self.close(fig)

    def test_summary_and_split_bars_draw_mean_and_run_dots(self):
        values = plots.indirect_metrics_values(synthetic_indirect_entries())
        fig, axes = self.figure(2, 4)
        try:
            plots.draw_metrics(values, axes)
            for axis, name in ((axes[1, 0], "Reliability@d*"), (axes[1, 1], "d*"),
                               (axes[1, 3], "NullSplit@d*")):
                panel = next(item for item in values["panels"] if item["panel"] == name)
                dots = [line for line in axis.lines if line.get_marker() == "o"]
                defined = [entry for entry in panel["entries"] if entry["mean_defined"]]
                self.assertEqual(len(dots), len(defined))
                for line, entry in zip(dots, defined):
                    np.testing.assert_allclose(
                        np.asarray(line.get_ydata(), dtype=float),
                        np.asarray([value for value, flag in zip(entry["values"],
                                                                entry["defined"]) if flag],
                                   dtype=float))
        finally:
            self.close(fig)


class MetricsReplayTests(unittest.TestCase):
    def test_direct_values_round_trip_and_replay_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        from PIL import Image
        entries = synthetic_direct_entries()
        values = plots.metrics_values(entries)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            values_csv = tmp / "values.csv"
            plots.write_plot_values(values, values_csv)
            self.assertEqual(canonical(plots.read_plot_values(values_csv)), canonical(values))
            first, second = tmp / "first.png", tmp / "second.png"
            plots.render_metrics(entries, first)
            plots.replay_metrics(values_csv, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with Image.open(first) as image:
                self.assertGreater(image.width, image.height * 2)

    def test_indirect_values_round_trip_and_replay_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        from PIL import Image
        entries = synthetic_indirect_entries()
        values = plots.indirect_metrics_values(entries)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            values_csv = tmp / "values.csv"
            plots.write_plot_values(values, values_csv)
            self.assertEqual(canonical(plots.read_plot_values(values_csv)), canonical(values))
            first, second = tmp / "first.png", tmp / "second.png"
            plots.render_indirect_metrics(entries, first)
            plots.replay_metrics(values_csv, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with Image.open(first) as image:
                self.assertGreater(image.width, image.height)
                self.assertGreater(image.width, 1400)

    def test_recorded_layout_drives_the_replay_grid(self):
        entries = synthetic_direct_entries()
        values = plots.metrics_values(entries)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "values.csv"
            plots.write_plot_values(values, path)
            rebuilt = plots.read_plot_values(path)
            self.assertEqual((rebuilt["layout"]["rows"], rebuilt["layout"]["columns"]), (1, 4))
            import pandas as pd
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            frame.loc[0, "x"] = ""
            frame.to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "layout is not recorded"):
                plots.read_plot_values(path)

    def test_committed_fixture_is_the_synthetic_values_and_replays_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        from PIL import Image
        self.assertTrue(FIXTURE_PNG.exists(),
                        "run the fixture build command to create metrics_fixture.png")
        self.assertTrue(FIXTURE_CSV.exists(),
                        "run the fixture build command to create metrics_fixture.csv")
        self.assertEqual(canonical(plots.read_plot_values(FIXTURE_CSV)),
                         canonical(plots.metrics_values(synthetic_direct_entries())))
        with tempfile.TemporaryDirectory() as tmp:
            replay = Path(tmp) / "replay.png"
            plots.replay_metrics(FIXTURE_CSV, replay)
            self.assertEqual(replay.read_bytes(), FIXTURE_PNG.read_bytes())
        with Image.open(FIXTURE_PNG) as image:
            self.assertGreater(image.width, image.height * 2)


if __name__ == "__main__":
    unittest.main()
