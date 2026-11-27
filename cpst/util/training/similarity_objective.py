"""Add one bounded Delta similarity term to the unit-weight reconstruction loss."""

from fractions import Fraction
from pathlib import Path
import sys

from torch.nn.functional import normalize

RECONSTRUCTION_EXECUTION = Path(__file__).resolve().parents[2] / "05_tsg_loss/05.3_mse-l1_mse-x/execution"
sys.path.append(str(RECONSTRUCTION_EXECUTION))

from .reconstruction_objective import ReconstructionObjective
from .delta_objective import EPSILON

METRICS = ("pcc", "cos")
ALPHAS = ("0", "0.1", "0.2", "0.3", "0.4", "0.5")
ALPHA_LABELS = {Fraction(label): label for label in ALPHAS}
REFERENCE_ALPHA = "0"
SMOKE_ALPHAS = (REFERENCE_ALPHA, "0.5")
METRIC_LABELS = {"pcc": "PCC", "cos": "cosine"}


class SimilarityObjective:
    """Penalize row-wise dissimilarity between paired matched-control Delta vectors."""

    def __init__(self, metric, alpha):
        if metric not in METRICS or alpha not in ALPHA_LABELS:
            raise ValueError(f"unsupported Delta similarity coordinate: {metric}, {alpha}")
        self.metric = metric
        self.alpha = alpha
        self.predecessor = ReconstructionObjective(Fraction(1))

    def identity(self):
        return {**self.predecessor.identity(), "alpha": ALPHA_LABELS[self.alpha],
                "similarity_metric": self.metric, "similarity_weight": float(self.alpha),
                "similarity_epsilon": EPSILON}

    def __call__(self, native_loss, control, control_rec, mu1, logvar1,
                 measured, treatment_rec, mu2, logvar2, predicted, mu_pred, logvar_pred):
        base = self.predecessor(native_loss, control, control_rec, mu1, logvar1,
                                measured, treatment_rec, mu2, logvar2, predicted, mu_pred, logvar_pred)
        if not self.alpha:
            return base
        predicted_delta, measured_delta = predicted - control, measured - control
        if self.metric == "pcc":
            predicted_delta = predicted_delta - predicted_delta.mean(dim=-1, keepdim=True)
            measured_delta = measured_delta - measured_delta.mean(dim=-1, keepdim=True)
        predicted_unit = normalize(predicted_delta, dim=-1, eps=EPSILON)
        measured_unit = normalize(measured_delta, dim=-1, eps=EPSILON)
        similarity = (predicted_unit * measured_unit).sum(dim=-1).clamp(-1, 1)
        penalty = ((1 - similarity) / 2).mean()
        return base[0] + float(self.alpha) * penalty, *base[1:]
