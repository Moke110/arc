"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if parent.name == "cpst")))
from util.training.architecture_roster import (
    BASE_OBJECTIVE,
    BOTH,
    CONTROL,
    DIRECT_DELTA,
    FILM,
    HEAD_PARAMETERS,
    MLP,
    NO_X2,
    NO_X2_OBJECTIVE,
    RESIDUAL,
    SHAPES,
    SHARED,
    coordinates,
    objective_label,
)
