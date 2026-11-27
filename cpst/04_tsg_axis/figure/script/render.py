"""Render the TranSiGen prediction-axis comparison from its plotted CSVs."""

import argparse
from pathlib import Path
import sys

FIGURE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIGURE.parents[1] / "util/output/api"))

from cpst_util.comparison import render_comparison

ARMS = {f"{family}-{width}k": (f"{family}-{width}k", hue, 5 - width)
        for family, hue in (("hvg", 220), ("denoise", 20)) for width in (3, 4, 5)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=FIGURE / "data")
    parser.add_argument("--output", type=Path, default=FIGURE)
    args = parser.parse_args()
    render_comparison(args.data, args.output, arms=ARMS, reference="denoise-5k",
                      title="TranSiGen prediction axes")


if __name__ == "__main__":
    main()
