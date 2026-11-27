"""Exercise the shared plotting contract through its rendered artists."""

import unittest

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from cpst_util import palette, plotting


class PlottingTests(unittest.TestCase):
    def test_palette_coordinates_preserve_rgb_between_roles(self):
        self.assertEqual(palette.HUES, tuple(range(20, 341, 40)))
        self.assertEqual(palette.DEPTHS, tuple(range(4)))
        for hue in palette.HUES:
            colors = [palette.color_pair(hue, depth) for depth in palette.DEPTHS]
            self.assertEqual(len({pair.main for pair in colors}), 4)
            for pair in colors:
                self.assertEqual(pair.main[:3], pair.secondary[:3])
                self.assertEqual((pair.main[3], pair.secondary[3]), (.90, .40))
        for channel, expected in zip(palette.color_pair(20, 0).main[:3], (.5808, .2464, .0792)):
            self.assertAlmostEqual(channel, expected)
        with self.assertRaises(ValueError):
            palette.color_pair(0)
        with self.assertRaises(ValueError):
            palette.color_pair(20, 4)

    def setUp(self):
        self.figure, self.axis = plt.subplots()
        self.addCleanup(plt.close, self.figure)

    def test_bars_and_grouped_bars_use_role_colors(self):
        colors = palette.color_pair(220, 1)
        bar, = plotting.bar(self.axis, [1], [2.5], hue=220, depth=1, width=.4)
        self.assertEqual((bar.get_x(), bar.get_height()), (.8, 2.5))
        self.assertEqual(bar.get_facecolor(), colors.secondary)
        self.assertEqual(bar.get_edgecolor(), colors.main)

        bars = plotting.bar(self.axis, [0, 1], [2, 3], hue=[220, 140], depth=1)
        self.assertEqual(bars[1].get_facecolor(), palette.color_pair(140, 1).secondary)
        with self.assertRaises(ValueError):
            plotting.bar(self.axis, [0, 1], [2, 3], hue=[20, 60, 100])

        groups = plotting.grouped_bar(self.axis, [1, 2], [
            ([1, 3], 220, 1, "Group A"),
            ([2, 4], 220, 3, "Group B"),
        ], width=.8)
        self.assertEqual([group.get_label() for group in groups], ["Group A", "Group B"])
        self.assertAlmostEqual(groups[0][0].get_center()[0], .8)
        self.assertAlmostEqual(groups[1][0].get_center()[0], 1.2)
        self.assertEqual(groups[1][1].get_height(), 4)
        self.assertEqual(groups[1][0].get_edgecolor(), palette.color_pair(220, 3).main)

    def test_line_range_and_error_bars_share_a_coordinate(self):
        colors = palette.color_pair(140, 3)
        x, values = [1, 2], [2, 3]
        line, = plotting.line(self.axis, x, values, hue=140, depth=3, linestyle="--")
        self.assertEqual(list(line.get_ydata()), values)
        self.assertEqual(line.get_color(), colors.main)
        self.assertEqual(line.get_linestyle(), "--")

        band = plotting.range_band(self.axis, x, [1, 2], [3, 4], hue=140, depth=3)
        self.assertEqual(tuple(band.get_facecolors()[0]), colors.secondary)
        self.assertEqual((band.get_paths()[0].vertices[:, 1].min(),
                          band.get_paths()[0].vertices[:, 1].max()), (1, 4))

        errors = plotting.errorbar(self.axis, x, values, yerr=[.5, 1], capsize=3,
                                   hue=140, depth=3, fmt="none")
        collection, = errors.lines[2]
        self.assertEqual(tuple(collection.get_colors()[0]), colors.main)
        self.assertEqual(collection.get_segments()[0].tolist(), [[1, 1.5], [1, 2.5]])

    def test_learning_curve_draws_both_series_for_every_run(self):
        pairs = plotting.learning_curve(self.axis, [
            ([1, 2, 3], [3, 2, 1], [5, 3, 2]),
            ([1, 2], [4, 1], [4, 3]),
        ])
        self.assertEqual(len(self.axis.lines), 4)
        self.assertEqual([list(training.get_ydata()) for training, _ in pairs],
                         [[3, 2, 1], [4, 1]])
        self.assertEqual([list(validation.get_xdata()) for _, validation in pairs],
                         [[1, 2, 3], [1, 2]])
        for training, validation in pairs:
            self.assertEqual(training.get_color(), palette.color_pair(220).main)
            self.assertEqual(validation.get_color(), palette.color_pair(140).main)
            self.assertEqual(training.get_linestyle(), validation.get_linestyle())
        self.assertNotEqual(pairs[0][0].get_linestyle(), pairs[1][0].get_linestyle())
        self.assertEqual(self.axis.get_legend_handles_labels()[1], ["Training", "Validation"])

    def test_scatter_and_bubbles_keep_positions_sizes_and_color(self):
        colors = palette.color_pair(300, 2)
        dots = plotting.scatter(self.axis, [1, 2], [3, 4], hue=300, depth=2, s=8)
        self.assertEqual(dots.get_offsets().tolist(), [[1, 3], [2, 4]])
        self.assertEqual(tuple(dots.get_facecolors()[0]), colors.main)

        bubbles = plotting.bubble(self.axis, [1, 2], [3, 4], [30, 60], hue=300, depth=2)
        self.assertEqual(bubbles.get_sizes().tolist(), [30, 60])
        self.assertEqual(tuple(bubbles.get_facecolors()[0]), colors.main)
        self.assertEqual(tuple(bubbles.get_edgecolors()[0]), colors.main)

    def test_heatmap_draws_an_ordered_color_scale_with_numeric_ticks(self):
        image = plotting.heatmap(self.axis, [[0, .5], [.8, 1]],
                                 hue=220, vmin=0, vmax=1, colorbar_label="F1")
        self.assertEqual(image.get_array().tolist(), [[0, .5], [.8, 1]])
        self.assertEqual(image.get_clim(), (0, 1))
        self.assertEqual(image.get_cmap().N, 4)
        self.assertEqual(image.get_cmap()(0.), palette.color_pair(220, 3).main)
        self.assertEqual(image.get_cmap()(1.), palette.color_pair(220, 0).main)
        color_axis = self.figure.axes[1]
        self.assertEqual(color_axis.get_ylabel(), "F1")
        self.assertIn(0., color_axis.get_yticks())
        self.assertIn(1., color_axis.get_yticks())

    def test_replicate_dots_use_one_group_coordinate_and_each_defined_value(self):
        dots = plotting.replicate_dots(self.axis, 2, [4, -1, float("nan"), -1],
                                      hue=180, depth=3, spread=.25)
        self.assertEqual(dots.get_offsets().tolist(), [[1.75, -1], [2, -1], [2.25, 4]])
        self.assertEqual([tuple(color) for color in dots.get_facecolors()],
                         [palette.color_pair(180, 3).main])
        single = plotting.replicate_dots(self.axis, 6, [0], hue=180, depth=3)
        self.assertEqual(single.get_offsets().tolist(), [[6., 0.]])

    def test_replicate_bar_plots_sample_mean_individual_values_and_sample_sd(self):
        bars, dots, errors = plotting.replicate_bar(
            self.axis, 2, [4, float("nan"), -2, 1, float("inf")],
            hue=180, depth=3, width=.5, spread=.1, label="A",
        )
        bar, = bars
        self.assertEqual((bar.get_height(), bars.get_label()), (1, "A"))
        self.assertEqual(dots.get_offsets().tolist(), [[1.9, -2], [2, 1], [2.1, 4]])
        self.assertEqual(bar.get_edgecolor(), palette.color_pair(180, 3).main)
        collection, = errors.lines[2]
        self.assertEqual(collection.get_segments()[0].tolist(), [[2, -2], [2, 4]])
        self.assertEqual(tuple(collection.get_colors()[0]), bar.get_edgecolor())

        single, single_dot, single_error = plotting.replicate_bar(self.axis, 4, [0])
        self.assertEqual(single[0].get_height(), 0)
        self.assertEqual(single_dot.get_offsets().tolist(), [[4, 0]])
        self.assertIsNone(single_error)
        self.assertEqual(plotting.replicate_bar(self.axis, 6, [float("nan")]),
                         (None, None, None))

    def test_invalid_replicate_spread_leaves_no_partial_marks(self):
        with self.assertRaises(ValueError):
            plotting.replicate_bar(self.axis, 2, [1, 3], spread=-1)
        self.assertEqual((len(self.axis.patches), len(self.axis.collections), len(self.axis.containers)),
                         (0, 0, 0))

    def test_focus_y_contains_replicates_uncertainty_and_padding(self):
        lower, upper = plotting.focus_y(self.axis, [[10, 14], [float("nan"), 15]])
        self.assertAlmostEqual(lower, 8.58873, places=5)
        self.assertAlmostEqual(upper, 15.58284, places=5)
        self.assertEqual(self.axis.get_ylim(), (lower, upper))

    def test_focus_y_preserves_reference_and_expands_constant_values(self):
        lower, upper = plotting.focus_y(self.axis, [[-3, -3]], reference=0)
        self.assertLess(lower, -3)
        self.assertGreater(upper, 0)
        self.assertEqual(plotting.focus_y(self.axis, [[5, 5]]), (4.5, 5.5))
        self.assertEqual(plotting.focus_y(self.axis, [[], [float("nan")]]), (4.5, 5.5))

    def test_histogram_retains_bin_edges_and_counts(self):
        counts, edges, bars = plotting.histogram(self.axis, [0, .3, .75, 1],
                                                 bins=[0, .5, .7, 1], hue=60)
        self.assertEqual(counts.tolist(), [2, 0, 2])
        self.assertEqual(edges.tolist(), [0, .5, .7, 1])
        self.assertEqual(bars[0].get_edgecolor(), palette.color_pair(60).main)
        self.assertEqual(bars[0].get_facecolor(), palette.color_pair(60).secondary)
