"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
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
