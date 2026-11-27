"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from util.training.native import (
    Adapter,
    DoseView,
    SCRIPT_HASH,
    SOURCE_DIGEST,
    TRANSFORM_KIND,
    WIDTH,
    WIDTHS,
    adapter,
    configuration,
    configure_backend,
    fit_transform,
    positive_doses,
    transform_doses,
)
