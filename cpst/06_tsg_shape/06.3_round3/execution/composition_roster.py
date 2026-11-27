"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from util.training.composition_roster import (
    COMBINATIONS,
    CONTROL,
    DIRECT_DELTA,
    DIRECT_DELTA_FILM,
    FILM,
    HEAD_PARAMETERS,
    NORMALIZED_FILM_SHAPES,
    NO_X2,
    NO_X2_FILM,
    NO_X2_SHAPES,
    REUSED_SEED,
    ROUND_TWO_SHAPES,
    SEEDS,
    SHAPES,
    SMOKE_FOLDS,
    SMOKE_SEEDS,
    coordinates,
    fit_coordinates,
    objective_label,
)
