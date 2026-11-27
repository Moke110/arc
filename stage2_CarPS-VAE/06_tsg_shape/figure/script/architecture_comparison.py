"""Arm styles and references for the four Step 06 comparison publications."""

import argparse
from pathlib import Path
import sys
from typing import NamedTuple

SHAPE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SHAPE.parent / "util/output/api"))

from stage2_util.comparison import render_comparison

RENDERER = "figure/script/render_comparison.py"


class Study(NamedTuple):
    title: str
    ticket: int
    reference: str
    seeds: tuple
    arms: dict


def _styles(colors):
    return {arm: (arm, hue, depth) for arm, (hue, depth) in colors.items()}


STUDIES = {
    "06.1_round1": Study("TranSiGen expression widths and depths", 350, "original", (42,), _styles({
        "original": (20, 2), "widen_z-110": (60, 2), "widen_z-120": (100, 2), "widen_z-130": (140, 2),
        "widen_enc-2k": (180, 3), "widen_enc-4k": (180, 2), "widen_enc-6k": (180, 1),
        "deepen_enc": (220, 2), "deepen_dec": (260, 2), "deepen_both": (300, 2),
    })),
    "06.2_round2": Study("TranSiGen prediction architectures", 353, "original", (42,), _styles({
        "original": (20, 2), "output-delta": (60, 2), "output-delta-nox2": (60, 3),
        "z-interact_MLP": (140, 2), "z-interact_FiLM": (140, 3), "shared-z": (220, 2),
        "z-res": (260, 2), "shared-z-res": (300, 2),
    })),
    "06.3_round3": Study("TranSiGen architecture combinations", 356, "0_original", (17, 42, 2026, 43), _styles({
        "0_original": (20, 2), "1_output-delta": (60, 2), "2_output-delta-no-x2": (60, 3),
        "3_z-interact_FiLM": (140, 2), "4_1+3": (220, 2), "5_2+3": (220, 3),
    })),
    "06.4_round4": Study("TranSiGen measured encoder targets", 359, "XX-XX", (17, 42, 2026), _styles({
        "XX-XX": (20, 2), "XX-XXΔ": (60, 2), "X-XΔ": (140, 2), "XΔ-XΔ": (220, 2), "XΔ-XΔΔ": (300, 2),
    })),
}


def render_study(step):
    figure = step / "figure"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=figure / "data")
    parser.add_argument("--output", type=Path, default=figure)
    args = parser.parse_args()
    study = STUDIES[step.name]
    render_comparison(args.data, args.output, arms=study.arms, reference=study.reference, title=study.title)
