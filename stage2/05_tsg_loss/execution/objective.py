"""TranSiGen objectives for the matched-plate Delta loss comparison."""

from fractions import Fraction
import math

import torch
import torch.nn.functional as F

FORMS = ("mse", "cos", "pcc", "l1")
ALPHAS = (Fraction(1, 4), Fraction(1, 3), Fraction(1, 2))
EPSILON = 1e-8
COMPONENTS = ("control_mse", "treatment_mse", "delta", "control_kl", "treatment_kl", "prediction_kl")


def delta_loss(form, predicted, measured):
    """Reduce Delta across conditions and the ordered gene axis."""
    if form == "mse":
        return F.mse_loss(predicted, measured, reduction="sum")
    if form == "l1":
        return F.l1_loss(predicted, measured, reduction="sum")
    if form == "pcc":
        predicted = predicted - predicted.mean(dim=-1, keepdim=True)
        measured = measured - measured.mean(dim=-1, keepdim=True)
    elif form != "cos":
        raise ValueError(f"unknown Delta loss: {form}")
    predicted = F.normalize(predicted, dim=-1, eps=EPSILON)
    measured = F.normalize(measured, dim=-1, eps=EPSILON)
    similarity = (predicted * measured).sum(dim=-1)
    return (1 - similarity).sum() * predicted.shape[-1]


class Objective:
    """Use native reconstructions and KL terms with the registered Delta term."""

    def __init__(self, form, alpha):
        if form not in FORMS or alpha not in ALPHAS:
            raise ValueError(f"unsupported Delta loss coordinate: {form}, {alpha}")
        self.form = form
        self.alpha = alpha

    def identity(self):
        return {
            "form": self.form, "alpha": str(self.alpha),
            "delta_weight": float(3 * self.alpha),
            "reconstruction_weight": float(Fraction(3, 2) * (1 - self.alpha)),
            "kl_weight": 0.1,
            "similarity_epsilon": EPSILON,
        }

    def __call__(self, native_loss, control, control_rec, mu1, logvar1,
                 measured, treatment_rec, mu2, logvar2, predicted, mu_pred, logvar_pred):
        native = native_loss(control, control_rec, mu1, logvar1, measured, treatment_rec,
                             mu2, logvar2, predicted, mu_pred, logvar_pred)
        if self.form == "mse" and self.alpha == Fraction(1, 3):
            return native
        _, control_mse, treatment_mse, mse_delta, kl1, kl2, kl_prediction = native
        delta = (mse_delta if self.form == "mse"
                 else delta_loss(self.form, predicted - control, measured - control))
        reconstruction_weight = float(Fraction(3, 2) * (1 - self.alpha))
        delta_weight = float(3 * self.alpha)
        loss = (0.1 * (kl1 + kl2 + kl_prediction)
                + delta_weight * delta + reconstruction_weight * (control_mse + treatment_mse))
        return loss, control_mse, treatment_mse, delta, kl1, kl2, kl_prediction


class ComponentTrace:
    """Accumulate the raw objective components without changing the fit's RNG."""

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
