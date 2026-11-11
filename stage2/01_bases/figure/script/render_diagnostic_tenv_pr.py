#!/usr/bin/env python3
"""Render the published tenv diagnostic PR panel from its plotted-values CSV alone.

This is a publication script, not scientific code: it reads only
``figure/data/diagnostic_tenv_pr.csv`` and writes the PNG. It never opens output
products, matrices, checkpoints or truth, and it recomputes no metric. The CSV
carries the curve coordinates, the fixed eligible-group counts, the displayed
equal-group macro AP values and the plotted title, so this renderer needs no
other input.

    python3 stage2/01_bases/figure/script/render_diagnostic_tenv_pr.py \
        --csv stage2/01_bases/figure/data/diagnostic_tenv_pr.csv \
        --png stage2/01_bases/figure/diagnostic_tenv_pr.png

The committed PNG has sha256
35d14266283b296c1585af3d58f0bacb15a3d7f854cde4972c433375bbdd911f; re-rendering
from the CSV reproduces it byte for byte (recorded by ticket #156).
"""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402  (backend must be set first)
import pandas as pd  # noqa: E402

METRICS = ("AUPRC", "AUPRC_positive", "AUPRC_negative")
RECORDS = ("figure", "panel", "summary", "curve")


def load(path):
    """Read and validate the complete plotted-values table."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("plotted-values CSV is empty")
    kinds = {row["record"] for row in rows}
    if not set(RECORDS) <= kinds:
        raise ValueError(f"plotted-values CSV is incomplete: {sorted(kinds)}")
    figures = [row for row in rows if row["record"] == "figure"]
    if len(figures) != 1 or not figures[0]["title"]:
        raise ValueError("exactly one figure record with a title is required")
    panels = {row["panel"]: row["panel_title"] for row in rows if row["record"] == "panel"}
    if tuple(panels) != METRICS or any(not value for value in panels.values()):
        raise ValueError("panel records must define AUPRC, AUPRC_positive and AUPRC_negative")
    summaries = {metric: {} for metric in METRICS}
    series_order = []
    for row in rows:
        if row["record"] == "summary":
            if row["panel"] not in summaries or not row["series"]:
                raise ValueError("summary record without a panel or series")
            summaries[row["panel"]][row["series"]] = float(row["macro_ap"])
            if row["series"] not in series_order:
                series_order.append(row["series"])
    curves = {}
    for row in rows:
        if row["record"] == "curve":
            curves.setdefault((row["panel"], row["series"]), []).append(row)
    runs = []
    for series in series_order:
        summary = {metric: summaries[metric][series] for metric in METRICS}
        drawn = {}
        for metric in METRICS:
            part = curves.get((metric, series), [])
            if not part:
                raise ValueError(f"no curve values for {series}/{metric}")
            frame = pd.DataFrame(
                [
                    {
                        "recall": float(row["recall"]),
                        "precision": float(row["precision"]),
                        "contributing_groups": int(row["contributing_groups"]),
                        "eligible_groups": int(row["eligible_groups"]),
                    }
                    for row in part
                ]
            )
            drawn[metric] = frame
        runs.append({"label": series, "summary": summary, "curves": drawn})
    return figures[0]["title"], panels, runs


def render(title, panels, runs, output):
    """Draw the three-panel figure with the accepted publication geometry."""
    output = Path(output)
    if output.suffix != ".png":
        raise ValueError("a .png output path is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), sharex=True, sharey=True,
                             layout="constrained")
    colors = plt.get_cmap("tab10").colors
    try:
        for index, run in enumerate(runs):
            label = run["label"]
            for name, ax in zip(METRICS, axes):
                curve, value = run["curves"][name], run["summary"][name]
                ax.step(curve.recall, curve.precision, where="post",
                        color=colors[index % len(colors)], linewidth=1.8,
                        label=f"{label} · macro AP {value:.3f}")
        for name, ax in zip(METRICS, axes):
            ax.set(title=panels[name], xlabel="Recall", xlim=(0, 1.01), ylim=(0, 1.04))
            ax.grid(alpha=.18)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(loc="upper center", bbox_to_anchor=(.5, -.23), fontsize=8,
                          frameon=False)
        axes[0].set_ylabel("Precision")
        fig.suptitle(title, fontsize=13)
        fig.savefig(output, dpi=160)
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, required=True, help="figure/data plotted values")
    parser.add_argument("--png", type=Path, required=True, help="figure PNG to write")
    args = parser.parse_args()
    render(*load(args.csv), args.png)
    print(f"rendered {args.png} from {args.csv}")


if __name__ == "__main__":
    main()
