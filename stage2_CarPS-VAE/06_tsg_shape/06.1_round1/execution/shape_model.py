"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parents[1]
LOSS = STAGE / "05_tsg_loss"
sys.path[:0] = [str(STAGE / "05_tsg_loss/05.2_mse_l1/execution"), str(STAGE / "05_tsg_loss/05.1_loss_metrics/execution")]
from run_loss import axis_run
from util.training.shape_model import (
    ShapeAdapter,
)
