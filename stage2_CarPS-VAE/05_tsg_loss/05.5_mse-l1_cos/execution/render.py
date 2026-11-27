"""Render Step 05.5 figures from the accepted plot CSVs."""

from pathlib import Path
import sys

SHARED = Path(__file__).resolve().parents[2] / "similarity"
sys.path.insert(0, str(SHARED))

from render_similarity import main


if __name__ == "__main__":
    main("cos")
