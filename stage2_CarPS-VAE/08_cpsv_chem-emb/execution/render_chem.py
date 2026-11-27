"""Render Step 08 nine held-out metrics from the published fold and seed score CSV."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

STAGE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STAGE / "util/output/api"))
from stage2_util import palette, plotting

ARMS = ("ECFP4C", "ECFP6C", "FCFP4C", "FCFP6C", "E3FP", "KPGT", "UniMol")
HUES = palette.HUES[:len(ARMS)]
METRICS = (
    ("JAC_top200", "JAC top200"),
    ("JAC_pos100", "JAC pos100"),
    ("JAC_neg100", "JAC neg100"),
    ("Δ PCC", "Δ PCC"),
    ("Δ Spearman", "Δ Spearman"),
    ("Δ RMSE", "Δ RMSE"),
    ("Δ E-distance", "Δ E-distance"),
    ("Δ MMD", "Δ MMD"),
    ("Δ R²", "Δ R²"),
)
LOWER_IS_BETTER = {"Δ RMSE", "Δ E-distance", "Δ MMD"}
FIGURES = {"value": "delta-metrics.png", "score_difference": "delta-metrics-paired-difference.png"}


def read_scores(path):
    scores = {}
    pairs = set()
    for row in plotting.read_rows(path):
        coordinate = int(row["fold"]), int(row["seed"])
        key = (row["arm"], row["metric"], *coordinate)
        values = {field: float(row[field]) if row[field] else math.nan for field in FIGURES}
        if key in scores or any(math.isinf(value) for value in values.values()):
            raise ValueError(f"invalid or duplicate score for {key}")
        scores[key] = values
        pairs.add(coordinate)

    bounded = {(0, 17)}
    formal = set(product(range(3), (17, 42, 2026)))
    if pairs not in (bounded, formal) or set(scores) != {
        (arm, metric, *pair) for arm, (metric, _), pair in product(ARMS, METRICS, pairs)
    }:
        raise ValueError("plotting scores do not cover every registered arm, metric, fold and seed")
    return scores, sorted(pairs)


def render(data, output):
    scores, pairs = read_scores(data / "fold-scores.csv")
    output.mkdir(exist_ok=True, parents=True)
    positions = range(len(ARMS))
    for field, filename in FIGURES.items():
        figure, axes = plt.subplots(3, 3, figsize=(18, 11), layout="constrained")
        paired = field == "score_difference"
        for (metric, label), ax in zip(METRICS, axes.flat, strict=True):
            direction = "↓" if metric in LOWER_IS_BETTER else "↑"
            ax.set_title(f"{label} {direction}", fontweight="bold")
            groups = [[scores[arm, metric, *pair][field] for pair in pairs] for arm in ARMS]
            for index, (values, hue) in enumerate(zip(groups, HUES, strict=True)):
                bars, _, _ = plotting.replicate_bar(ax, index, values, hue=hue, width=.68)
                if bars is None:
                    ax.text(index, .02, "NA", ha="center", fontsize=7, transform=ax.get_xaxis_transform())
            reference = 0 if paired or metric == "Δ R²" else None
            lower, upper = plotting.focus_y(ax, groups, reference=reference)
            if lower <= 0 <= upper:
                ax.axhline(0, color="#666666", linewidth=.8, zorder=1.5)
            ax.set_xticks(positions, ARMS, rotation=30, ha="right", fontsize=8)
            ax.set_xlim(-.5, len(ARMS) - .5)
            ax.grid(axis="y", color="#777777", alpha=.25, zorder=1)
            ax.spines[["top", "right"]].set_visible(False)

        title = "Bars are means · Dots are fold and seed coordinates · Error bars are sample SD"
        if not paired:
            title += "\nY-axis ranges follow the data and may omit zero."
        figure.suptitle(title)
        figure.supylabel("Paired difference from KPGT" if paired else "Held-out score")
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
