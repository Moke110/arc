"""Render HVG-DEG coverage from the plotted CSV."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from stage2_util.palette import color
from stage2_util.plotting import save_figure


def main():
    figure_path = Path(__file__).resolve().parents[1]
    coverage = pd.read_csv(figure_path / "data/disease-DEG1000-coverage_bars.csv")
    groups = coverage.pivot(index="disease", columns="axis_size", values="coverage")
    groups = groups.reindex(coverage.disease.drop_duplicates())
    palette = {2000: "orange", 3000: "green", 4000: "purple", 5000: "magenta"}
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": "#ededed",
    })

    x = np.arange(len(groups))
    width = 0.15
    center = (len(groups.columns) - 1) / 2
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    for i, size in enumerate(groups.columns):
        ax.bar(x + (i - center) * width, groups[size], width,
               color=color(palette[size], 2), label=f"{size // 1000}k")
    title = "disease top-DEG1000 coverage"
    ax.set(title=title, ylabel=title, ylim=(0, 1))
    ax.set_xticks(x, groups.index)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.legend(title="Axis size", frameon=False, ncol=len(groups.columns), loc="upper left")
    save_figure(fig, figure_path / "disease-DEG1000-coverage.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
