"""Render Step 06.1 architecture aggregates from the accepted metrics CSV."""

import argparse
from itertools import product
import math
from pathlib import Path

from render_shape import FIGURES, LABELS, MEAN_HUE, metric_direction, plt, plotting, read_table


def render(data, output):
    summary = read_table(data, ("shape", "metric"))
    shapes = tuple(dict.fromkeys(shape for shape, _ in summary))
    if (not shapes or shapes[0] != "widen_z-100"
            or set(summary) != set(product(shapes, LABELS))):
        raise ValueError("aggregate plotting CSV requires a complete architecture-by-metric grid")

    positions = range(len(shapes))
    for field, filename in FIGURES.items():
        difference = field == "score_difference"
        figure, plots = plt.subplots(3, 3, figsize=(15, 10), layout="constrained")
        for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
            metric_direction(plot, metric, label)
            means = [summary[shape, metric][field] for shape in shapes]
            plotting.scatter(plot, positions, means, hue=MEAN_HUE, marker="D", s=30, zorder=3)
            for position, mean in zip(positions, means, strict=True):
                if math.isnan(mean):
                    plot.text(position, .02, "NA", ha="center", fontsize=7, transform=plot.get_xaxis_transform())
            if difference:
                plot.axhline(0, color="#666666", linewidth=0.8, zorder=2)
            plot.set_xticks(positions, shapes, rotation=30, ha="right", fontsize=8)
            plot.set_xlim(-0.5, len(shapes) - 0.5)
            plot.grid(axis="y", color="#777777", alpha=0.25, zorder=1)
            plot.spines[["top", "right"]].set_visible(False)
        figure.supylabel("Mean paired difference from the original architecture" if difference
                         else "Equal-fold mean held-out score")
        plotting.save_figure(figure, output / filename, dpi=180)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
