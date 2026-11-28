"""Render the TranSiGen module comparison from its plotted CSVs."""

import argparse
from pathlib import Path
import sys

FIGURE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIGURE.parents[1] / "util/output/api"))

from cpst_util.comparison import COMPARISON_MODES, render_comparison

ARMS = {
    "XX-XX": ("XX-XX", 20, 2),
    "XX-XXΔ": ("XX-XXΔ", 60, 2),
    "X-XΔ": ("X-XΔ", 140, 2),
    "XΔ-XΔ": ("XΔ-XΔ", 220, 2),
    "XΔ-XΔΔ": ("XΔ-XΔΔ", 300, 2),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=FIGURE / "data")
    parser.add_argument("--output", type=Path, default=FIGURE)
    parser.add_argument("--mode", choices=COMPARISON_MODES)
    args = parser.parse_args()
    render_comparison(args.data, args.output, arms=ARMS, reference="XX-XX",
                      title="TranSiGen modules", mode=args.mode)


if __name__ == "__main__":
    main()
