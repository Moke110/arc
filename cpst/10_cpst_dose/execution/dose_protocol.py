"""Registered arms and coordinates for the CarPS-Transcriptome dose-input comparison."""

from itertools import product
from functools import cache
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
REFERENCE = "no-dose"
TRANSIGEN = "TranSiGen"
ARMS = (TRANSIGEN, REFERENCE, "Xmol-scaling", "Xmol-FiLM", "Zmol-scaling", "Zmol-FiLM",
        "Xdose-num", "Xdose-MLP", "Zdose-num", "Zdose-MLP")
sys.path.insert(0, str(STAGE / "util/output/api"))
from cpst_util.tenv_folds import COLD_DOSE_DRUG, HELD_DOSES, open_tenv_folds

SCENARIO = COLD_DOSE_DRUG


@cache
def protocol():
    return open_tenv_folds().protocol(COLD_DOSE_DRUG)


AXIS = "denoise-5k"

DEFAULT_STEP03 = STAGE / "03_data_train"
DEFAULT_REPRESENTATIONS = STAGE / "util/output/representations"
DEFAULT_SOURCE = STAGE.parent / "pert_model/TranSiGen/source"


def fold_numbers():
    return tuple(fold.number for fold in protocol().folds)


def coordinates(bounded=False):
    """Return dose arm, cold-dose-drug fold and training seed in registered order."""
    folds, seeds = fold_numbers(), protocol().seeds
    # Group 0 at 1 µM exercises a crossed fold absent from the historical study.
    return tuple(product(ARMS, (1,) if bounded else folds, seeds[:1] if bounded else seeds))


def fold_identity(fold):
    split = protocol().folds[fold]
    return {"parent_group": split.parent_group, "held_dose_um": split.held_dose_um,
            "protocol_sha256": protocol().release_sha256}


def representation(arm):
    if arm not in ARMS:
        raise ValueError(f"unknown dose arm: {arm}")
    return "KPGT" if arm == TRANSIGEN else "UniMol"
