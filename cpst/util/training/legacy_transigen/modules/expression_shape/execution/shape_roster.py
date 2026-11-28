"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if parent.name == "cpst")))
from util.training.shape_roster import (
    CONTROL,
    NEW_WIDTHS,
    PRIOR_NAMES,
    SHAPES,
    SMOKE_SHAPES,
    coordinates,
    dimensions,
)
