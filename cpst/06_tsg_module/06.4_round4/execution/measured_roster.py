"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
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
