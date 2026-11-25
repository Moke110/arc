"""Render the twelve loss and weight panels from the accepted epoch CSV."""

import argparse
from collections import defaultdict
from itertools import product
import math
from pathlib import Path
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP.parent / "util/output/api"))

from stage2_util.palette import color
from stage2_util.plotting import read_rows, save_figure

LOSSES = {"mse": "MSE", "cos": "Cos", "pcc": "PCC", "l1": "L1"}
ALPHAS = ("1/4", "1/3", "1/2")
FOLDS = (0, 1, 2)
COLORS = tuple(color(family, 3) for family in ("blue", "orange", "green"))
EPOCH_SOURCE = "epoch-objectives.csv"
FIGURE_NAME = "learning-curves.png"


def read_curves(path):
    """Read one continuous, checkpoint-selected trace per formal coordinate."""
    curves = defaultdict(list)
    for row in read_rows(path):
        key = row["loss"], row["alpha"], int(row["fold"])
        identity = row["run_sha256"], int(row["seed"]), int(row["best_epoch"]), int(row["stopping_epoch"])
        epoch = int(row["epoch"])
        train, valid, mse = (float(row[field]) for field in ("train_loss", "validation_loss", "validation_mse"))
        if not all(math.isfinite(value) for value in (train, valid, mse)) or min(train, valid) <= 0:
            raise ValueError(f"invalid objective in {key} epoch {epoch}")
        curves[key].append((identity, epoch, train, valid, mse))

    if curves.keys() != set(product(LOSSES, ALPHAS, FOLDS)):
        raise ValueError("epoch CSV must contain all 36 formal loss, weight and fold coordinates")

    run_digests = set()
    for key, points in curves.items():
        identity = points[0][0]
        run_digest, seed, best, stop = identity
        if ([point[1] for point in points] != list(range(1, stop + 1)) or not 1 <= stop <= 1000):
            raise ValueError(f"missing or padded epochs for {key}")
        if (seed != 42 or not re.fullmatch(r"[0-9a-f]{64}", run_digest) or run_digest in run_digests
                or any(point[0] != identity for point in points)):
            raise ValueError(f"inconsistent run identity for {key}")
        run_digests.add(run_digest)

        best_epoch, best_loss = 0, math.inf
        for _, epoch, _, valid, _ in points:
            if valid < best_loss:
                best_epoch, best_loss = epoch, valid
            if epoch - best_epoch == 30:
                break
        if epoch != stop or best_epoch != best or (stop < 1000 and stop - best != 30):
            raise ValueError(f"checkpoint selection or stopping epoch differs for {key}")
    return curves


def render(data, output):
    curves = read_curves(data / EPOCH_SOURCE)
    fig, plots = plt.subplots(4, 3, figsize=(13, 11), sharex=True, sharey=True,
                             layout="constrained")
    for (loss, alpha), plot in zip(product(LOSSES, ALPHAS), plots.flat, strict=True):
        for fold in FOLDS:
            _, epochs, train, valid, _ = zip(*curves[loss, alpha, fold], strict=True)
            plot.plot(epochs, train, linestyle="--", color=COLORS[fold], linewidth=1.0)
            plot.plot(epochs, valid, color=COLORS[fold], linewidth=1.2)
        plot.set_title(f"{LOSSES[loss]} · α = {alpha}", fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)

    handles = [Line2D([], [], color=COLORS[fold], label=f"Fold {fold}") for fold in FOLDS]
    handles.extend((Line2D([], [], color="#333333", linestyle="--", label="Train"),
                    Line2D([], [], color="#333333", label="Validation")))
    fig.legend(handles=handles, ncols=5, loc="outside upper center", frameon=False)
    fig.supxlabel("Epoch")
    fig.supylabel("Objective per condition (log scale)")
    save_figure(fig, output / FIGURE_NAME, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
