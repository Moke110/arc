"""Render HVG-DEG coverage from the plotted CSV."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

FIGURE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIGURE.parents[2] / "util/output/api"))

from stage2_util.plotting import grouped_bar, save_figure

HVG_HUE = 220


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=FIGURE / "data")
    parser.add_argument("--output", type=Path, default=FIGURE)
    args = parser.parse_args()
    coverage = pd.read_csv(args.data / "disease-DEG1000-coverage_bars.csv")
    groups = coverage.pivot(index="disease", columns="axis_size", values="coverage")
    groups = groups.reindex(coverage.disease.drop_duplicates())
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": "#ededed",
    })

    x = np.arange(len(groups))
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    bars = [(groups[size], HVG_HUE, 5 - size // 1000, f"{size // 1000}k")
            for size in groups.columns]
    grouped_bar(ax, x, bars, width=.60)
    ax.set(title="Disease top-DEG1000 coverage", ylabel="Coverage", ylim=(0, 1.02))
    ax.set_xticks(x, groups.index)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.legend(title="Axis size", frameon=False, ncol=len(groups.columns), loc="upper left")
    save_figure(fig, args.output / "disease-DEG1000-coverage.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
