"""Compatibility imports for the shared Stage 2 implementation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from util.training.similarity_objective import (
    ALPHAS,
    ALPHA_LABELS,
    METRICS,
    METRIC_LABELS,
    RECONSTRUCTION_EXECUTION,
    REFERENCE_ALPHA,
    SMOKE_ALPHAS,
    SimilarityObjective,
)
