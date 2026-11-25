"""Exercise the shared plotting contract through its rendered artists."""

import unittest

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from stage2_util import palette, plotting


class PlottingTests(unittest.TestCase):
    def test_palette_coordinates_preserve_rgb_between_roles(self):
        self.assertEqual(palette.HUES, tuple(range(20, 341, 40)))
        self.assertEqual(palette.DEPTHS, tuple(range(5)))
        for hue in palette.HUES:
            colors = [palette.color_pair(hue, depth) for depth in palette.DEPTHS]
            self.assertEqual(len({pair.main for pair in colors}), 5)
            for pair in colors:
                self.assertEqual(pair.main[:3], pair.secondary[:3])
                self.assertEqual((pair.main[3], pair.secondary[3]), (.90, .40))
        for channel, expected in zip(palette.color_pair(20, 0).main[:3], (.432, .176, .048)):
            self.assertAlmostEqual(channel, expected)
        with self.assertRaises(ValueError):
            palette.color_pair(0)
        with self.assertRaises(ValueError):
            palette.color_pair(20, 5)

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
        self.assertEqual(image.get_cmap().N, 5)
        self.assertEqual(image.get_cmap()(0.), palette.color_pair(220, 4).main)
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

    def test_histogram_retains_bin_edges_and_counts(self):
        counts, edges, bars = plotting.histogram(self.axis, [0, .3, .75, 1],
                                                 bins=[0, .5, .7, 1], hue=60)
        self.assertEqual(counts.tolist(), [2, 0, 2])
        self.assertEqual(edges.tolist(), [0, .5, .7, 1])
        self.assertEqual(bars[0].get_edgecolor(), palette.color_pair(60).main)
        self.assertEqual(bars[0].get_facecolor(), palette.color_pair(60).secondary)
