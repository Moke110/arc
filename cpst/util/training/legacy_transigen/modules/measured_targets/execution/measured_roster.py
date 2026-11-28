"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if parent.name == "cpst")))
from util.training.measured_roster import (
    CONTROL,
    CONTROL_ONLY,
    DIRECT_DELTA,
    MEASURED_DELTA,
    MEASURED_SEPARATE,
    MEASURED_SHARED,
    SEEDS,
    SHAPES,
    SMOKE_FOLDS,
    SMOKE_SEEDS,
    SOURCE_SHAPES,
    coordinates,
    fit_coordinates,
    objective_label,
)
