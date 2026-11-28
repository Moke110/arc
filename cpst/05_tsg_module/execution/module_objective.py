"""Common Step 05 loss for the five accepted TranSiGen modules."""

from functools import partial

import torch
from torch.nn import functional as F

from util.training.measured_roster import CONTROL_ONLY, MEASURED_DELTA, SHAPES


def objective_identity(shape):
    if shape not in SHAPES:
        raise ValueError(f"Unknown TranSiGen module: {shape}")
    measured = shape != CONTROL_ONLY
    return {
        "control_mse_weight": 1.0,
        "control_kl_weight": 0.1,
        "measured_target": "signed Delta" if measured else None,
        "measured_mse_weight": float(measured),
        "measured_kl_weight": 0.1 if measured else 0.0,
        "prediction_kl_weight": 0.1 if measured else 0.0,
        "delta_mse_weight": 1.0,
        "delta_l1_weight": 0.0,
    }


def gaussian_kl(mean, logvar, target_mean, target_logvar):
    return 0.5 * torch.sum(
        target_logvar - logvar - 1
        + (logvar.exp() + (mean - target_mean).square()) * (-target_logvar).exp()
    )


def module_loss(model, shape, control, control_rec, mu1, logvar1, treated, treated_rec,
                mu2, logvar2, predicted, mu_pred, logvar_pred):
    observed_delta = treated - control
    control_mse = F.mse_loss(control_rec, control, reduction="sum")
    prediction_mse = F.mse_loss(predicted - control, observed_delta, reduction="sum")
    zero = control_mse.new_zeros(())
    control_kl = gaussian_kl(mu1, logvar1, zero, zero)
    measured_mse = measured_kl = alignment_kl = zero
    if shape != CONTROL_ONLY:
        if shape in MEASURED_DELTA:
            reconstructed_delta, mu2, logvar2 = model.reconstruct_measured_delta(observed_delta)
        else:
            reconstructed_delta = treated_rec - control
        measured_mse = F.mse_loss(reconstructed_delta, observed_delta, reduction="sum")
        measured_kl = gaussian_kl(mu2, logvar2, zero, zero)
        alignment_kl = gaussian_kl(mu_pred, logvar_pred, mu2, logvar2)
    loss = 0.1 * (control_kl + measured_kl + alignment_kl) + control_mse + measured_mse + prediction_mse
    return loss, control_mse, measured_mse, prediction_mse, control_kl, measured_kl, alignment_kl, zero


def configure_module_objective(model, shape):
    model.loss = partial(module_loss, model, shape)
