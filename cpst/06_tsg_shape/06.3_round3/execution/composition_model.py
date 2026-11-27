"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

STAGE = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(STAGE / "05_tsg_loss/05.2_mse_l1/execution"), str(STAGE / "05_tsg_loss/05.1_loss_metrics/execution")]
from run_loss import axis_run
from util.training.composition_model import (
    CompositionAdapter,
    NormalizedFiLMPrediction,
    ROUND_TWO,
    STEP,
    configure_objective,
    objective_identity,
)
