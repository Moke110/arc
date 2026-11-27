"""Render equal-fold Delta scores and signed differences against unit reconstruction weight."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

LOSS = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(LOSS.parent / "util/output/api"), str(LOSS / "execution")]

from loss_plot_style import FOLD_SCORES, read_replicates, render_weight_panels
from stage2_util import plotting

from sweep_protocol import ALPHAS, SMOKE_ALPHAS, REFERENCE_ALPHA

LABELS = {
    "JAC_top200": "JAC top200",
    "JAC_pos100": "JAC pos100",
    "JAC_neg100": "JAC neg100",
    "Δ PCC": "Δ PCC",
    "Δ Spearman": "Δ Spearman",
    "Δ RMSE": "Δ RMSE",
    "Δ E-distance": "Δ E-distance",
    "Δ MMD": "Δ MMD",
    "Δ R²": "Δ R²",
}
FIGURES = {"value": "delta-metrics.png",
           "score_difference": "delta-metrics-improvement.png"}


def read_scores(path, alphas):
    scores = {}
    for row in plotting.read_rows(path):
        key = row["alpha"], row["metric"]
        if key in scores:
            raise ValueError(f"duplicate score in complete sweep: {key}")
        record = {name: float(row[name]) if row[name] else math.nan
                  for name in (*FIGURES, "reference_weight_score")}
        if any(math.isinf(value) for value in record.values()):
            raise ValueError(f"infinite score for {key}")
        scores[key] = record
    if set(scores) != set(product(alphas, LABELS)):
        raise ValueError("metric summary must contain the complete sweep")
    for key, record in scores.items():
        baseline = scores[REFERENCE_ALPHA, key[1]]["value"]
        actual = record["reference_weight_score"], record["score_difference"]
        expected = baseline, record["value"] - baseline
        if any(a != b and not (math.isnan(a) and math.isnan(b))
               for a, b in zip(actual, expected, strict=True)):
            raise ValueError(f"signed score difference disagrees with the unit-weight score: {key}")
    return scores


def render(data, output, *, smoke=False):
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    scores = read_scores(data / "delta-metrics.csv", alphas)
    replicates = read_replicates(data / FOLD_SCORES, scores, coordinates=("alpha",),
                                 reference=(REFERENCE_ALPHA,), folds=(0,) if smoke else (0, 1, 2))
    for field, filename in FIGURES.items():
        render_weight_panels(
            replicates, field, output / filename, alphas, LABELS,
            xlabel="Shared reconstruction MSE weight α", reference=REFERENCE_ALPHA,
            rotate=not smoke,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    render(args.data, args.output, smoke=args.smoke)


if __name__ == "__main__":
    main()
