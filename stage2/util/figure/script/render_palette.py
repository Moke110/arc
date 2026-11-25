"""Render the Stage 2 palette and examples of each color role."""

import argparse
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

UTILITY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(UTILITY / "output/api"))

from stage2_util import palette, plotting


def render(output):
    figure = plt.figure(figsize=(12.9, 8.42), facecolor="#fdfdfd")
    ink = "#253044"

    def label(x, y, text, size, **kwargs):
        figure.text(x / 1290, 1 - y / 842, text,
                    color=ink, fontsize=size, va="center", **kwargs)

    label(48, 62, "Stage 2 palette", 23, weight="bold")
    label(48, 103, "Nine hues × five depths, with two opacities per coordinate", 13)
    label(46, 167, "DEPTH / HSL", 10.5, weight="bold")
    for column, hue in enumerate(palette.HUES):
        x = 195 + 116 * column
        label(x + 53, 168, f"{hue}°", 12, weight="bold", ha="center")
        for depth in palette.DEPTHS:
            y = 187 + 76 * depth
            colors = palette.color_pair(hue, depth)
            for side, rgba in enumerate(colors):
                figure.add_artist(Rectangle(
                    ((x + 53 * side) / 1290, 1 - (y + 57) / 842),
                    53 / 1290, 57 / 842, transform=figure.transFigure,
                    facecolor=rgba, edgecolor="none",
                ))
    for depth, (lightness, saturation) in enumerate(palette.DEPTH_LEVELS):
        y = 202 + 76 * depth
        label(46, y, f"Depth {depth + 1}", 11, weight="bold")
        label(46, y + 26, f"L {lightness:.0%} · S {saturation:.0%}", 10)
    label(195, 586, f"Left  {palette.MAIN_OPACITY:.0%} main     "
                    f"Right  {palette.SECONDARY_OPACITY:.0%} secondary", 11.5)
    figure.add_artist(Line2D([48 / 1290, 1243 / 1290], [1 - 618 / 842] * 2,
                             transform=figure.transFigure, color="#d3dbe7", linewidth=.6))
    label(48, 655, "Same coordinate in a chart", 14, weight="bold")

    for index, hue in enumerate((20, 140, 220, 300)):
        x = 195 + 251 * index
        axis = figure.add_axes((x / 1290, 1 - 773 / 842, 173 / 1290, 75 / 842))
        axis.set(xlim=(0, 6), ylim=(0, 1))
        axis.set_axis_off()
        depth = 2
        plotting.bar(axis, [1.1], [.65], hue=hue, depth=depth, width=1.5, linewidth=1.8)
        plotting.errorbar(axis, [1.1], [.47], yerr=[[0], [.37]], fmt="none",
                          hue=hue, depth=depth, capsize=6, capthick=1.6, linewidth=1.6)
        plotting.scatter(axis, [1.1], [.47], hue=hue, depth=depth, s=26)
        epochs = [2.75, 3.75, 4.75, 5.75]
        values = [.22, .53, .46, .78]
        plotting.range_band(axis, epochs, [value - .13 for value in values],
                            [value + .19 for value in values], hue=hue, depth=depth)
        plotting.line(axis, epochs, values, hue=hue, depth=depth, linewidth=1.7)
        label(x + 86.5, 796, f"({hue}°, {depth + 1})", 10.5, ha="center")

    plotting.save_figure(figure, output, dpi=100)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=UTILITY / "figure/palette.png")
    render(parser.parse_args().output)


if __name__ == "__main__":
    main()
