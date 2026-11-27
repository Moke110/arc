"""Render tenv diagnostics using the recorded CSVs."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

SCORE_COLUMNS = (
    ("p_tenv", r"$p_{\mathrm{tenv}}$  ·  tenv top-DEG occurrences", 220),
    ("p_disease", "disease top-DEG1000 occurrences", 20),
    ("score", r"Score = $p_{\mathrm{tenv}}\times 2^{p_{\mathrm{disease}}}$", 140),
)
STEP = Path(__file__).resolve().parents[2]
RETAINED_HUE, EXCLUDED_HUE = 220, 20
SIZE_HUES = {1000: 20, 2000: 100, 3000: 140, 4000: 220, 5000: 300}
DISEASE_DEPTHS = {"nHF": 0, "ICM": 1, "genetic_cardiomyopathy": 2, "TOX": 3}


def plot_exclusions(ax, values, excluded, bins, cutoff, cutoff_label, drawing):
    centers = (bins[:-1] + bins[1:]) / 2
    width = np.diff(bins) * .44
    counts, _ = np.histogram(values[~excluded], bins)
    drawing.bar(ax, centers - width / 2, counts, width=width, hue=RETAINED_HUE,
                label=f"Pass: {(~excluded).sum():,}")
    counts, _ = np.histogram(values[excluded], bins)
    n = int(excluded.sum())
    drawing.bar(ax, centers + width / 2, counts, width=width, hue=EXCLUDED_HUE,
                label=f"Excluded: {n:,} ({n / len(values):.1%})")
    ax.axvline(cutoff, color="#666666", linestyle="--", linewidth=1, label=cutoff_label)
    ax.set(ylabel="Genes", yscale="log", ylim=(0.8, None))


def variance_boundary(ax, genes, metric, k, drawing):
    measured = genes.dmso_mean_log2_cpm1.gt(0)
    ratios = genes.loc[measured, metric] / genes.loc[measured, "dmso_mean_log2_cpm1"]
    values = np.log2(1 + ratios)
    excluded = genes.loc[measured, metric + "_above_boundary"]
    bins = np.linspace(0, values.max(), 31)
    plot_exclusions(ax, values, excluded, bins, np.log2(1 + k), f"SD / mean = {k:.3f}", drawing)
    ax.set_xlabel("log₂(1 + SD / mean DMSO log₂(1 + CPM))")
    ax.legend(frameon=False, fontsize=9, loc="upper right")


def render_scores(data_path, output, exclusions, histogram, drawing):
    scores = pd.read_csv(data_path / "diagnostics_scores.csv")
    scores = scores.merge(exclusions, on=["gene_index", "symbol"], validate="one_to_one")
    scores = scores.loc[~scores.excluded]
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4), layout="constrained")
    fig.suptitle(f"Top-DEG presence and weighted score · {len(scores):,} retained genes")
    first_axis = histogram.loc[histogram.axis_size.eq(histogram.axis_size.min())]
    for ax, (key, label, hue) in zip(axes, SCORE_COLUMNS, strict=True):
        values = scores[key]
        bins = first_axis.loc[first_axis.metric.eq(key)]
        edges = np.r_[bins.bin_left, bins.bin_right.iloc[-1]]
        drawing.histogram(ax, values, bins=edges, hue=hue)
        if key == "p_disease":
            ax.set_xticks(np.arange(values.max() + 1))
        ax.set(xlabel=label, ylabel="Genes", xscale="linear", yscale="log")
        ax.text(0.96, 0.95, f"Genes at zero: {(values == 0).sum():,}",
                transform=ax.transAxes, ha="right", va="top")
    drawing.save_figure(fig, output / "deg-scores.png", dpi=180)
    plt.close(fig)


def render_scores_in_axes(output, histogram, drawing):
    sizes = sorted(histogram.axis_size.unique())
    fig, axes = plt.subplots(len(sizes), len(SCORE_COLUMNS), figsize=(12.8, 12.0),
                             sharex="col", sharey=True, layout="constrained")
    ymax = 1.3 * histogram.retained_genes.max()
    for size, row in zip(sizes, axes, strict=True):
        for ax, (key, label, _) in zip(row, SCORE_COLUMNS, strict=True):
            bins = histogram.loc[histogram.axis_size.eq(size) & histogram.metric.eq(key)]
            centers = (bins.bin_left + bins.bin_right) / 2
            width = (bins.bin_right - bins.bin_left) * .40
            size_hue = SIZE_HUES[size]
            drawing.bar(ax, centers - width / 2, bins.retained_genes - bins.axis_genes,
                        width=width, hue=size_hue, depth=3, label="Outside this axis")
            drawing.bar(ax, centers + width / 2, bins.axis_genes,
                        width=width, hue=size_hue, depth=0, label="Selected in this axis")
            ax.set(yscale="log", ylim=(0.8, ymax))
            if key == "p_disease":
                ax.set_xticks(np.arange(int(bins.bin_right.max() + 0.5)))
            if size == sizes[-1]:
                ax.set_xlabel(label)
        row[0].set_ylabel(f"G{size}\nGenes")
    fig.suptitle("Top-DEG scores by gene-axis membership")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    drawing.save_figure(fig, output / "score-in-axis.png", dpi=180)
    plt.close(fig)


def render_coverage(data_path, output, drawing):
    coverage = pd.read_csv(data_path / "disease-DEG1000-coverage.csv")
    groups = coverage.pivot(index="disease", columns="axis_size", values="coverage")
    groups = groups.reindex(coverage.disease.drop_duplicates())
    x = np.arange(len(groups))
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    depths = [DISEASE_DEPTHS[disease] for disease in groups.index]
    bars = [(groups[size], SIZE_HUES[size], depths, f"{size // 1000}k")
            for size in groups.columns]
    drawing.grouped_bar(ax, x, bars, width=.75)
    names = [name.replace("genetic_cardiomyopathy", "Genetic cardiomyopathy") for name in groups.index]
    ax.set(title="Denoise-DEG · disease top-DEG1000 coverage", ylabel="Coverage", ylim=(0, 1))
    ax.set_xticks(x, names)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.legend(title="Axis size", frameon=False, ncol=len(groups.columns), loc="upper left")
    drawing.save_figure(fig, output / "disease-DEG1000-coverage.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    repo = next(parent for parent in STEP.parents if (parent / ".git").exists())
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument("--step", type=Path, default=STEP)
    parser.add_argument("--data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / "cpst/util/output/api"))
    from cpst_util import plotting

    figure_path = args.step / "figure"
    data_path = args.data or figure_path / "data"
    output = args.output or figure_path
    genes = pd.read_csv(data_path / "diagnostics_genes.csv.gz")
    exclusions = pd.read_csv(data_path / "diagnostics_gene_exclusions.csv")
    genes = genes.merge(exclusions, on=["gene_index", "symbol"], validate="one_to_one")
    slopes = pd.read_csv(data_path / "diagnostics_variance_boundaries.csv")
    slopes = slopes.loc[slopes.target.eq(.10)].set_index("metric").k
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": "#ededed",
    })

    fig, axes = plt.subplots(1, 3, figsize=(15.9, 4.7), layout="constrained")
    ax = axes[0]
    plot_exclusions(ax, genes.detection_fraction, genes.low_detection,
                    np.linspace(0, 1, 41), .10, "Nonzero fraction = 0.10", plotting)
    ax.set(title="Detection across nonblank wells", xlabel="Fraction of wells with nonzero counts")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for ax, (metric, title) in zip(axes[1:], (
        ("plate_sd", "Between-plate DMSO variation"),
        ("within_plate_sd", "Within-plate DMSO variation"),
    ), strict=True):
        ax.set_title(title)
        variance_boundary(ax, genes, metric, slopes[metric], plotting)
    plotting.save_figure(fig, output / "denoise.png", dpi=180)
    plt.close(fig)

    histogram = pd.read_csv(data_path / "score-in-axis.csv")
    render_scores(data_path, output, exclusions, histogram, plotting)
    render_scores_in_axes(output, histogram, plotting)
    render_coverage(data_path, output, plotting)


if __name__ == "__main__":
    main()
