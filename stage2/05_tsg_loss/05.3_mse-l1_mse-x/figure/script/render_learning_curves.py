"""Render the reconstruction MSE weight panels from the accepted epoch CSV."""

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

LOSS = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(LOSS.parent / "util/output/api"), str(LOSS / "figure/script"),
               str(Path(__file__).resolve().parents[2] / "execution")]

from loss_plot_style import FOLD_HUES, curve_legend_handles
from stage2_util import plotting

from sweep_protocol import ALPHAS, SMOKE_ALPHAS

FOLDS = (0, 1, 2)
EPOCH_SOURCE = "epoch-objectives.csv"
FIGURE_NAME = "learning-curves.png"


def read_curves(path, alphas, folds, cap):
    """Read one continuous, checkpoint-selected trace per registered coordinate."""
    curves = defaultdict(list)
    for row in plotting.read_rows(path):
        key = row["alpha"], int(row["fold"])
        identity = row["run_sha256"], int(row["seed"]), int(row["best_epoch"]), int(row["stopping_epoch"])
        epoch = int(row["epoch"])
        train, valid, mse = (float(row[field]) for field in ("train_loss", "validation_loss", "validation_mse"))
        if not all(math.isfinite(value) for value in (train, valid, mse)) or min(train, valid) <= 0:
            raise ValueError(f"invalid objective in {key} epoch {epoch}")
        curves[key].append((identity, epoch, train, valid, mse))

    if curves.keys() != set(product(alphas, folds)):
        raise ValueError("epoch CSV must contain the complete weight and fold roster")

    run_digests = set()
    for key, points in curves.items():
        identity = points[0][0]
        run_digest, seed, best, stop = identity
        if ([point[1] for point in points] != list(range(1, stop + 1)) or not 1 <= stop <= cap):
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
        if epoch != stop or best_epoch != best or (stop < cap and stop - best != 30):
            raise ValueError(f"checkpoint selection or stopping epoch differs for {key}")
    return curves


def render(data, output, *, smoke=False):
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    folds = (0,) if smoke else FOLDS
    curves = read_curves(data / EPOCH_SOURCE, alphas, folds, 2 if smoke else 1000)
    shape = (1, 2) if smoke else (3, 3)
    fig, plots = plt.subplots(*shape, figsize=(9, 4) if smoke else (13, 9),
                              sharex=True, sharey=True, layout="constrained")
    for alpha, plot in zip(alphas, plots.flat, strict=True):
        for fold in folds:
            _, epochs, train, valid, _ = zip(*curves[alpha, fold], strict=True)
            plotting.line(plot, epochs, train, hue=FOLD_HUES[fold], linestyle="--", linewidth=1.0)
            plotting.line(plot, epochs, valid, hue=FOLD_HUES[fold], linewidth=1.2)
        plot.set_title(f"Reconstruction α = {alpha}", fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)

    handles = curve_legend_handles(folds)
    fig.legend(handles=handles, ncols=len(handles), loc="outside upper center", frameon=False)
    fig.supxlabel("Epoch")
    fig.supylabel("Objective per condition (log scale)")
    plotting.save_figure(fig, output / FIGURE_NAME, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    render(args.data, args.output, smoke=args.smoke)


if __name__ == "__main__":
    main()
