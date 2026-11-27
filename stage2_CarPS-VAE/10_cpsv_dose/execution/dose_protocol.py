"""Registered arms and coordinates for the CarPS-VAE dose-input comparison."""

from itertools import product
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
REFERENCE = "no-dose"
TRANSIGEN = "TranSiGen"
ARMS = (TRANSIGEN, REFERENCE, "Xmol-scaling", "Xmol-FiLM", "Zmol-scaling", "Zmol-FiLM",
        "Xdose-num", "Xdose-MLP", "Zdose-num", "Zdose-MLP")
SEEDS = (17, 42, 2026)
SCENARIO = "cold_dose_drug"
HELD_DOSES = (0.2, 1.0, 10.0)
AXIS = "denoise-5k"

DEFAULT_STEP03 = STAGE / "03_data_train"
DEFAULT_REPRESENTATIONS = STAGE / "util/output/representations"
DEFAULT_SOURCE = STAGE.parent / "pert_model/TranSiGen/source"


def coordinates(bounded=False):
    """Return dose arm, cold-dose-drug fold and training seed in registered order."""
    return tuple(product(ARMS, (0,) if bounded else range(3), (17,) if bounded else SEEDS))


def representation(arm):
    if arm not in ARMS:
        raise ValueError(f"unknown dose arm: {arm}")
    return "KPGT" if arm == TRANSIGEN else "UniMol"
