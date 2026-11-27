"""Render the canonical Step 08 comparison from its three plot CSVs."""

import argparse
from pathlib import Path
import sys

STAGE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(STAGE / "util/output/api"))

from cpst_util.comparison import render_comparison
from cpst_util.palette import HUES
from cpst_util.representations import ARMS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    arms = {arm: (arm, hue, 0) for arm, hue in zip(ARMS, HUES[:len(ARMS)], strict=True)}
    render_comparison(args.data, args.output, arms=arms, reference="KPGT",
                      title="CarPS-VAE chemical representations")


if __name__ == "__main__":
    main()
