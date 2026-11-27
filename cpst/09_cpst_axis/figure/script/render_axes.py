"""Render the three Step 09 figures using only their plotted CSVs."""

import argparse
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "execution"))
from axis_protocol import AXES, REFERENCE_AXIS
from cpst_util.comparison import render_comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    families = {"hvg": ("HVG-DEG", 220), "denoise": ("Denoise-DEG", 20)}
    arms = {}
    for axis in AXES:
        family, width = axis.split("-")
        label, hue = families[family]
        arms[axis] = (f"{label} {width}", hue, 5 - int(width[0]))
    render_comparison(args.data, args.output, arms=arms, reference=REFERENCE_AXIS,
                      title="CarPS-VAE prediction axes")


if __name__ == "__main__":
    main()
