"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if parent.name == "cpst")))
from util.training.delta_objective import (
    ALPHAS,
    ALPHA_LABELS,
    COMPONENTS,
    ComponentTrace,
    EPSILON,
    FORMS,
    NATIVE_ALPHA,
    Objective,
    delta_loss,
)
