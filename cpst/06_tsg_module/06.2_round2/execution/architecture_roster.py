"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
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
