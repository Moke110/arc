"""Reference arms and palette for the five Step 05 comparison publications."""

import argparse
from pathlib import Path
import sys
from typing import NamedTuple

LOSS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(LOSS.parent / "util/output/api"))

from stage2_util.comparison import render_comparison
from stage2_util.palette import HUES

RENDERER = "figure/script/render_comparison.py"


class Study(NamedTuple):
    title: str
    receipt: str
    reference: str
    arms: dict
    similarity: str | None = None


def _weight_arms(weights):
    return {weight: (f"α {weight}", hue, 2) for weight, hue in zip(weights, HUES)}


def _similarity_study(label, metric):
    weights = ("0", "0.1", "0.2", "0.3", "0.4", "0.5")
    return Study(f"TranSiGen Delta {label}", "339-formal.json", "0", _weight_arms(weights), metric)


STUDIES = {
    "05.1_loss_metrics": Study("TranSiGen Delta loss forms", "330-formal.json", "mse:1", {
        f"{form}:{weight}": (f"{form.upper()} × {weight}", hue, depth)
        for form, hue in zip(("mse", "cos", "pcc", "l1"), HUES)
        for depth, weight in enumerate(("0.75", "1", "1.25"), start=1)
    }),
    "05.2_mse_l1": Study("TranSiGen additive Delta L1", "333-formal.json", "0",
                         _weight_arms(("0", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2"))),
    "05.3_mse-l1_mse-x": Study("TranSiGen reconstruction MSE", "336-formal.json", "1",
                              _weight_arms(("0", "0.5", "0.75", "0.9", "1", "1.1", "1.25", "1.5", "2"))),
    "05.4_mse-l1_pcc": _similarity_study("PCC", "pcc"),
    "05.5_mse-l1_cos": _similarity_study("cosine", "cos"),
}


def render_study(step):
    figure = step / "figure"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=figure / "data")
    parser.add_argument("--output", type=Path, default=figure)
    args = parser.parse_args()
    study = STUDIES[step.name]
    render_comparison(args.data, args.output, arms=study.arms, reference=study.reference, title=study.title)
