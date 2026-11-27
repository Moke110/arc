"""Render this study's three canonical figures from plotted CSVs."""

from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP.parent / "figure/script"))

from architecture_comparison import render_study


if __name__ == "__main__":
    render_study(STEP)
