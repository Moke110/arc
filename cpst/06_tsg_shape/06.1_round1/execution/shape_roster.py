"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from util.training.shape_roster import (
    CONTROL,
    NEW_WIDTHS,
    PRIOR_NAMES,
    SHAPES,
    SMOKE_SHAPES,
    coordinates,
    dimensions,
)
