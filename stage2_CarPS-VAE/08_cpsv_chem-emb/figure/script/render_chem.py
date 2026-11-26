"""Render Step 08 primary metrics from the published fold and seed score CSV."""

import argparse
import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

STAGE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(STAGE / "util/output/api"))
from stage2_util import palette, plotting

ARMS = ("ECFP4C", "ECFP6C", "FCFP4C", "FCFP6C", "E3FP", "KPGT", "UniMol")
HUES = palette.HUES[:len(ARMS)]
PRIMARY = ("Δ PCC", "Δ RMSE")
FIGURES = {"value": "primary-scores.png", "score_difference": "paired-primary.png"}


def render(data, output):
    with (data / "fold-scores.csv").open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row["metric"] in PRIMARY]
    keys = {(row["arm"], row["fold"], row["seed"], row["metric"]) for row in rows}
    pairs = {(row["fold"], row["seed"]) for row in rows}
    if (not pairs or len(keys) != len(rows)
            or keys != {(arm, fold, seed, metric) for arm in ARMS for fold, seed in pairs for metric in PRIMARY}):
        raise ValueError("plotting scores do not cover every paired metric coordinate")
    output.mkdir(exist_ok=True, parents=True)
    for value, filename in FIGURES.items():
        figure, axes = plt.subplots(1, 2, figsize=(11.2, 3.6), layout="constrained")
        paired = value == "score_difference"
        for ax, metric in zip(axes, PRIMARY, strict=True):
            groups = [[float(row[value]) for row in rows if row["arm"] == arm and row["metric"] == metric
                       and row[value]] for arm in ARMS]
            for index, (scores, hue) in enumerate(zip(groups, HUES, strict=True)):
                plotting.replicate_bar(ax, index, scores, hue=hue, width=.68)
            if paired:
                ax.axhline(0, color="0.4", linewidth=.8, linestyle="--")
            limits = plotting.focus_y(ax, groups, reference=0 if paired else None)
            ax.set(xticks=list(range(len(ARMS))), xticklabels=ARMS,
                   ylabel=f"{metric} minus KPGT" if paired else metric)
            ax.tick_params(axis="x", labelrotation=30)
            if not paired and limits and not limits[0] <= 0 <= limits[1]:
                ax.set_title("Y-axis omits zero", fontsize=8, loc="right")
        plotting.save_figure(figure, output / filename, dpi=160)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
