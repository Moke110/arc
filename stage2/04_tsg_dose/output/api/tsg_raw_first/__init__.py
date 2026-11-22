"""Active raw-first input, evaluation and reporting contracts for dose comparison."""

from pathlib import Path
import sys

STAGE2 = Path(__file__).resolve().parents[4]
for api in ("util/output/api", "02_axis_view/02.1_hvg-deg/output/api",
            "03_view/03.1_hvg-deg_views/output/api"):
    path = str(STAGE2 / api)
    if path not in sys.path:
        sys.path.insert(0, path)
