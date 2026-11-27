"""Registered coordinates for the CarPS-Transcriptome prediction-axis comparison."""

from itertools import product
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
REFERENCE_AXIS = "denoise-5k"
AXES = tuple(f"{family}-{width}k" for family in ("hvg", "denoise") for width in range(2, 6))
SEEDS = (17, 42, 2026)
ARM = "UniMol"
DISEASE_SIZE = 3481
DISEASE_OVERLAP = dict(zip(AXES, (1258, 1889, 2525, 3149, 906, 1225, 1375, 1513), strict=True))
SCOPES = ("axis", "disease")

sys.path[:0] = [
    str(STAGE / "04_tsg_axis/output/api"),
    str(STAGE / "08_cpst_chem-emb/execution"),
    str(STAGE / "03_data_train/output/api"),
    str(STAGE / "util/output/api"),
]


def coordinates(bounded=False):
    """Return each axis, cold-drug fold and training seed in registered order."""
    return tuple(product(AXES, (0,) if bounded else range(3), (17,) if bounded else SEEDS))
