"""Render tenv diagnostics using the recorded CSVs."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

SCORE_COLUMNS = (
    ("p_tenv", r"$p_{\mathrm{tenv}}$  ·  tenv top-DEG occurrences"),
    ("p_disease", "disease top-DEG1000 occurrences"),
    ("score", r"Score = $p_{\mathrm{tenv}}\times 2^{p_{\mathrm{disease}}}$"),
)


def plot_exclusions(ax, x, y, excluded, population_size, color, excluded_color):
    n = int(excluded.sum())
    ax.scatter(x[~excluded], y[~excluded], s=6, color=color, alpha=0.2, linewidths=0)
    ax.scatter(x[excluded], y[excluded], s=7, color=excluded_color,
               alpha=0.35, linewidths=0, label=f"Excluded: {n:,} ({n / population_size:.1%})")


def variance_boundary(ax, x, y, excluded, k, label, color, excluded_color):
    plot_exclusions(ax, x, y, excluded, int((x > 0).sum()), color, excluded_color)
    ymax = 1.05 * y.max()
    xmax = min(x.max(), ymax / k)
    ax.plot([0, xmax], [0, k * xmax], color=excluded_color,
            label=f"SD = {k:.3f} × mean log₂(1 + CPM)")
    ax.set(xlabel="Mean DMSO log₂(1 + CPM)", ylabel=label, ylim=(-0.05, ymax))
    ax.legend(frameon=False, fontsize=9, loc="upper right")


def render_scores(figure_path, exclusions, histogram, save_figure, palette):
    scores = pd.read_csv(figure_path / "data/diagnostics_scores.csv")
    scores = scores.merge(exclusions, on=["gene_index", "symbol"], validate="one_to_one")
    scores = scores.loc[~scores.excluded]
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4), layout="constrained")
    fig.suptitle(f"Top-DEG presence and weighted score · {len(scores):,} retained genes")
    first_axis = histogram.loc[histogram.axis_size.eq(histogram.axis_size.min())]
    for ax, (key, label), fill in zip(axes, SCORE_COLUMNS, palette, strict=True):
        values = scores[key]
        bins = first_axis.loc[first_axis.metric.eq(key)]
        edges = np.r_[bins.bin_left, bins.bin_right.iloc[-1]]
        ax.hist(values, bins=edges, color=fill)
        if key == "p_disease":
            ax.set_xticks(np.arange(values.max() + 1))
        ax.set(xlabel=label, ylabel="Genes", xscale="linear", yscale="log")
        ax.text(0.96, 0.95, f"Genes at zero: {(values == 0).sum():,}",
                transform=ax.transAxes, ha="right", va="top")
    save_figure(fig, figure_path / "deg-scores.png", dpi=180)
    plt.close(fig)


def render_scores_in_axes(figure_path, histogram, save_figure, light_blue, blue):
    sizes = sorted(histogram.axis_size.unique())
    fig, axes = plt.subplots(len(sizes), len(SCORE_COLUMNS), figsize=(12.8, 12.0),
                             sharex="col", sharey=True, layout="constrained")
    ymax = 1.3 * histogram.retained_genes.max()
    for size, row in zip(sizes, axes, strict=True):
        for ax, (key, label) in zip(row, SCORE_COLUMNS, strict=True):
            bins = histogram.loc[histogram.axis_size.eq(size) & histogram.metric.eq(key)]
            width = bins.bin_right - bins.bin_left
            ax.bar(bins.bin_left, bins.retained_genes, width=width, align="edge", color=light_blue)
            chosen = bins.axis_genes.gt(0)
            ax.bar(bins.bin_left[chosen], bins.axis_genes[chosen], width=width[chosen],
                   bottom=(bins.retained_genes - bins.axis_genes)[chosen],
                   align="edge", color=blue)
            ax.set(yscale="log", ylim=(0.8, ymax))
            if key == "p_disease":
                ax.set_xticks(np.arange(int(bins.bin_right.max() + 0.5)))
            if size == sizes[-1]:
                ax.set_xlabel(label)
        row[0].set_ylabel(f"G{size}\nGenes")
    fig.suptitle("Top-DEG scores by gene-axis membership")
    fig.legend(handles=[Patch(color=light_blue, label="Retained genes"),
                        Patch(color=blue, label="Selected in this axis")],
               loc="outside lower center", ncol=2, frameon=False)
    save_figure(fig, figure_path / "score-in-axis.png", dpi=180)
    plt.close(fig)


def render_coverage(figure_path, save_figure, palette):
    coverage = pd.read_csv(figure_path / "data/disease-DEG1000-coverage.csv")
    groups = coverage.pivot(index="disease", columns="axis_size", values="coverage")
    groups = groups.reindex(coverage.disease.drop_duplicates())
    x = np.arange(len(groups))
    width = 0.15
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    for i, (size, fill) in enumerate(zip(groups.columns, palette, strict=True)):
        ax.bar(x + (i - 2) * width, groups[size], width, color=fill, label=f"{size // 1000}k")
    names = [name.replace("genetic_cardiomyopathy", "Genetic cardiomyopathy") for name in groups.index]
    ax.set(title="disease top-DEG1000 coverage",
           ylabel="disease top-DEG1000 coverage", ylim=(0, 1))
    ax.set_xticks(x, names)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.legend(title="Axis size", frameon=False, ncol=len(groups.columns), loc="upper left")
    save_figure(fig, figure_path / "disease-DEG200-coverage.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--step", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / "stage2/util/output/api"))
    from stage2_util.palette import color
    from stage2_util.plotting import save_figure

    figure_path = args.step / "figure"
    genes = pd.read_csv(figure_path / "data/diagnostics_genes.csv.gz")
    exclusions = pd.read_csv(figure_path / "data/diagnostics_gene_exclusions.csv")
    genes = genes.merge(exclusions, on=["gene_index", "symbol"], validate="one_to_one")
    slopes = pd.read_csv(figure_path / "data/diagnostics_variance_boundaries.csv")
    slopes = slopes.loc[slopes.target.eq(.10)].set_index("metric").k
    s = color("blue", 2)
    t = color("orange", 2)
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": "#ededed",
    })

    fig, axes = plt.subplots(1, 3, figsize=(15.9, 4.7), layout="constrained")
    ax = axes[0]
    plot_exclusions(ax, np.log10(1 + genes.mean_cpm), genes.detection_fraction,
                    genes.low_detection, len(genes), s, t)
    ax.axhline(.10, color=t, label="Nonzero fraction = 0.10")
    ax.set(title="Abundance across nonblank wells", xlabel="log10(1 + mean CPM)",
           ylabel="Fraction of wells with nonzero counts", ylim=(-0.02, 1.02))
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for ax, (metric, title, ylabel) in zip(axes[1:], (
        ("plate_sd", "Between-plate DMSO variation", "Between-plate SD"),
        ("within_plate_sd", "Within-plate DMSO variation", "Within-plate residual SD"),
    ), strict=True):
        ax.set_title(title)
        variance_boundary(ax, genes.dmso_mean_log2_cpm1, genes[metric],
                          genes[metric + "_above_boundary"], slopes[metric], ylabel, s, t)
    save_figure(fig, figure_path / "denoise.png", dpi=180)
    plt.close(fig)

    histogram = pd.read_csv(figure_path / "data/score-in-axis.csv")
    render_scores(figure_path, exclusions, histogram, save_figure, (s, t, color("green", 2)))
    render_scores_in_axes(figure_path, histogram, save_figure, color("blue", 7), s)
    render_coverage(figure_path, save_figure,
                    [color(family, 2) for family in ("blue", "orange", "green", "purple", "magenta")])


if __name__ == "__main__":
    main()
