"""Render the Step 10 learning curves and metrics from plotted CSVs alone."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(STEP / "execution"), str(STEP.parent / "util/output/api")]
from dose_protocol import ARMS, REFERENCE, TRANSIGEN
from cpst_util.comparison import LOSS_LABEL, METRIC_PANELS
from cpst_util.plotting import focus_y, learning_curve, replicate_bar, save_figure

COORDINATE = ["arm", "fold", "seed"]
METRIC_KEY = [*COORDINATE, "metric"]
CONTRAST_ARMS = tuple(arm for arm in ARMS if arm not in (REFERENCE, TRANSIGEN))
STYLES = {
    TRANSIGEN: (260, 2), REFERENCE: (180, 2),
    "Xmol-scaling": (20, 1), "Xmol-FiLM": (60, 1),
    "Zmol-scaling": (20, 3), "Zmol-FiLM": (60, 3),
    "Xdose-num": (140, 1), "Xdose-MLP": (220, 1),
    "Zdose-num": (140, 3), "Zdose-MLP": (220, 3),
}


def read_panels(data):
    curves, scores, paired = (pd.read_csv(data / name, float_precision="round_trip")
                              for name in ("learning-curve.csv", "metrics.csv", "metrics-paired.csv"))
    if (curves.columns.tolist() != [*COORDINATE, "epoch", "train", "validation"]
            or scores.columns.tolist() != [*METRIC_KEY, "value"]
            or paired.columns.tolist() != [*METRIC_KEY, "reference", "arm_value", "reference_value", "value"]
            or set(scores.arm) != set(ARMS) or set(paired.arm) != set(CONTRAST_ARMS)):
        raise ValueError("unexpected Step 10 figure CSV columns or arm groups")
    metric_names = {name for name, _ in METRIC_PANELS}
    if (scores.duplicated(METRIC_KEY).any()
            or any(set(values) != metric_names for _, values in scores.groupby(COORDINATE).metric)
            or not np.isfinite(scores.value.dropna()).all()
            or not np.isfinite(curves[["train", "validation"]]).all().all()
            or set(map(tuple, curves[COORDINATE].to_numpy())) != set(map(tuple, scores[COORDINATE].to_numpy()))):
        raise ValueError("invalid Step 10 absolute metrics or curve coordinates")
    for _, trace in curves.groupby(COORDINATE):
        if trace.epoch.tolist() != list(range(1, len(trace) + 1)):
            raise ValueError("incomplete Step 10 epoch trace")
    expected = scores.loc[scores.arm.isin(CONTRAST_ARMS)].rename(columns={"value": "arm_value"})
    baseline = scores.loc[scores.arm.eq(REFERENCE), [*METRIC_KEY[1:], "value"]].rename(
        columns={"value": "reference_value"})
    expected = expected.merge(baseline, how="left", on=METRIC_KEY[1:], validate="many_to_one")
    expected["reference"] = REFERENCE
    expected["value"] = expected.arm_value - expected.reference_value
    if not paired.equals(expected[paired.columns]):
        raise ValueError("Step 10 contrasts differ from the no-dose fold–seed pairs")
    return curves, scores, paired


def format_axis(axis):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#eeeeee")
    axis.set_axisbelow(True)


def draw_curves(curves):
    figure, axes = plt.subplots(4, 3, figsize=(13.2, 11.2), layout="constrained")
    for axis, arm in zip(axes.flat, ARMS):
        selected = curves.loc[curves.arm.eq(arm)]
        traces = [(trace.epoch, trace.train, trace.validation)
                  for _, trace in selected.groupby(["fold", "seed"])]
        learning_curve(axis, traces)
        axis.set(title=arm, xlabel="Epoch", ylabel=LOSS_LABEL)
        axis.xaxis.set_major_locator(MaxNLocator(integer=True))
        format_axis(axis)
    for axis in list(axes.flat)[len(ARMS):]:
        axis.remove()
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside lower center", ncols=2, frameon=False)
    figure.suptitle("CarPS-VAE dose input · learning curves")
    return figure


def draw_metrics(table, arms, title):
    figure, axes = plt.subplots(2, 4, figsize=(18, 7.4), layout="constrained")
    for axis, (metric, label) in zip(axes.flat, METRIC_PANELS):
        data = table.loc[table.metric.eq(metric)]
        groups = [data.loc[data.arm.eq(arm), "value"].dropna().to_numpy() for arm in arms]
        for position, (arm, values) in enumerate(zip(arms, groups, strict=True)):
            hue, depth = STYLES[arm]
            if len(values):
                replicate_bar(axis, position, values, hue=hue, depth=depth, width=.66)
            else:
                axis.text(position, .04, "NA", ha="center", transform=axis.get_xaxis_transform())
        axis.set_xticks(range(len(arms)), arms, fontsize=8, rotation=32, ha="right", rotation_mode="anchor")
        axis.set_xlim(-.5, len(arms) - .5)
        axis.axhline(0, color="#777777", lw=.6)
        focus_y(axis, groups, reference=0)
        format_axis(axis)
        axis.set_title(label)
    axes[1, 3].remove()
    figure.suptitle(title)
    return figure


def render(data, output):
    curves, scores, paired = read_panels(data)
    figures = {
        "learning-curve.png": draw_curves(curves),
        "metrics.png": draw_metrics(scores, ARMS, "CarPS-VAE dose input · held-out scores · mean ± SD"),
        "metrics-paired.png": draw_metrics(
            paired, CONTRAST_ARMS, "CarPS-VAE dose input · paired differences vs no-dose · mean ± SD"),
    }
    for filename, figure in figures.items():
        save_figure(figure, output / filename, dpi=180)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
