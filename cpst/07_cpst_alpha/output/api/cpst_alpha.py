"""CarPS-Transcriptome alpha model and its fixed KPGT default."""

from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[2]
ROUND_FOUR = STEP.with_name("06_tsg_shape") / "06.4_round4"
sys.path[:0] = [str(STEP.parent), str(ROUND_FOUR / "execution")]

from measured_model import axis_run
from util.training.cpst_model import AlphaAdapter, KPGT_WIDTH, MEASURED_SHARED, PROJECTION_WIDTH, objective
