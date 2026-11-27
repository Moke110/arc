"""Keep native TranSiGen arithmetic and add the matched-control Delta L1 term."""

from fractions import Fraction
import math
from pathlib import Path
import sys

import torch

LOSS_COMPARISON = Path(__file__).resolve().parents[2] / "05_tsg_loss/05.1_loss_metrics/execution"

from .delta_objective import delta_loss

ALPHAS = ("0", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2")
ALPHA_LABELS = {Fraction(label): label for label in ALPHAS}
SMOKE_ALPHAS = ("0", "1")
COMPONENTS = ("control_mse", "treatment_mse", "delta_mse", "control_kl",
              "treatment_kl", "prediction_kl", "delta_l1")


class AdditiveObjective:
    def __init__(self, alpha):
        if alpha not in ALPHA_LABELS:
            raise ValueError(f"unsupported additive Delta L1 weight: {alpha}")
        self.alpha = alpha

    def identity(self):
        return {"alpha": ALPHA_LABELS[self.alpha], "delta_l1_weight": float(self.alpha),
                "delta_mse_weight": 1.0, "reconstruction_weight": 1.0, "kl_weight": 0.1}

    def __call__(self, native_loss, control, control_rec, mu1, logvar1,
                 measured, treatment_rec, mu2, logvar2, predicted, mu_pred, logvar_pred):
        total, *native_parts = native_loss(control, control_rec, mu1, logvar1,
                                           measured, treatment_rec, mu2, logvar2,
                                           predicted, mu_pred, logvar_pred)
        l1_delta = delta_loss("l1", predicted - control, measured - control)
        if self.alpha:
            total = total + float(self.alpha) * l1_delta
        return total, *native_parts, l1_delta


class ComponentTrace:
    """Accumulate detached raw batch sums and export per-condition means."""

    def __init__(self, loss):
        self.loss = loss
        self.count = 0
        self.sums = None

    def __call__(self, *inputs):
        result = self.loss(*inputs)
        parts = torch.stack([value.detach() for value in result[1:]]).to(torch.float64)
        self.sums = parts if self.sums is None else self.sums + parts
        self.count += len(inputs[0])
        return result

    def per_condition(self, partition):
        values = (self.sums / self.count).tolist()
        self.count, self.sums = 0, None
        if any(not math.isfinite(value) for value in values):
            raise ValueError(f"nonfinite {partition} objective component")
        return {f"{partition}_{name}": value for name, value in zip(COMPONENTS, values, strict=True)}
