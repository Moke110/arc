"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from util.training.additive_objective import (
    ALPHAS,
    ALPHA_LABELS,
    AdditiveObjective,
    COMPONENTS,
    ComponentTrace,
    LOSS_COMPARISON,
    SMOKE_ALPHAS,
)


sys.path.append(str(LOSS_COMPARISON))
