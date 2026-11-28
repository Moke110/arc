"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if parent.name == "cpst")))

STAGE = next(parent for parent in Path(__file__).resolve().parents if parent.name == "cpst")
sys.path[:0] = [str(STAGE / "util/training/legacy_transigen/loss/additive_delta/execution"), str(STAGE / "util/training/legacy_transigen/loss/loss_metrics/execution")]
from run_loss import axis_run
from util.training.composition_model import (
    CompositionAdapter,
    NormalizedFiLMPrediction,
    ROUND_TWO,
    STEP,
    configure_objective,
    objective_identity,
)
