"""TranSiGen measured encoder targets and signed Delta readouts for ADR 0098."""

from functools import partial
from pathlib import Path
import sys

import torch
from torch.nn import functional as F

from .measured_roster import CONTROL_ONLY, DIRECT_DELTA, MEASURED_DELTA, MEASURED_SHARED, SHAPES, SOURCE_SHAPES

STEP = Path(__file__).resolve().parents[2] / "06_tsg_shape/06.4_round4"
ROUND_THREE = STEP.with_name("06.3_round3")
sys.path.insert(0, str(ROUND_THREE / "execution"))

from .composition_model import CompositionAdapter
from .composition_model import configure_objective as configure_round_three
from .composition_model import objective_identity as round_three_objective


class MeasuredAdapter(CompositionAdapter):
    def __init__(self, source, transform, shape):
        if shape not in SHAPES:
            raise ValueError(f"unregistered Step 06.4 arm: {shape}")
        super().__init__(source, transform, SOURCE_SHAPES.get(shape, SOURCE_SHAPES[DIRECT_DELTA]))
        self.measured_shape = shape

    def network(self):
        shape = self.measured_shape
        original = super().network()
        if shape in SOURCE_SHAPES:
            return original

        class MeasuredNetwork(original):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                if shape in MEASURED_DELTA:
                    self.encoder_delta = self.encoder_x2
                    self.mu_delta = self.mu_z2
                    self.logvar_delta = self.logvar_z2
                    signed_measurement = self.decoder_x2[:-1]
                    if shape == MEASURED_SHARED:
                        self.decoder_delta = signed_measurement
                    else:
                        self.decoder_measured_delta = signed_measurement
                del self.encoder_x2, self.mu_z2, self.logvar_z2, self.decoder_x2

            @property
            def prediction_decoder(self):
                return self.decoder_delta

            @property
            def measurement_decoder(self):
                return self.decoder_delta if shape == MEASURED_SHARED else self.decoder_measured_delta

            def encode_x2(self, measured):
                # The shared loops supply measured X2 here and matched X1 to the objective.
                # New measured paths execute in that objective, where both targets are available.
                return None, None, None

            def decode_x2(self, latent):
                if latent is not None:
                    raise ValueError("the measured treated-expression decoder is absent")
                return None

            def reconstruct_measured_delta(self, delta):
                hidden = self.encoder_delta(delta)
                mean, logvar = self.mu_delta(hidden), self.logvar_delta(hidden)
                return self.measurement_decoder(self.sample_latent(mean, logvar)), mean, logvar

        return MeasuredNetwork


def objective_identity(shape):
    if shape in SOURCE_SHAPES:
        return round_three_objective(SOURCE_SHAPES[shape])
    return {
        "alpha": "1.75", "control_mse_weight": 1.0, "control_kl_weight": 0.1,
        "measured_target": None if shape == CONTROL_ONLY else "signed Delta",
        "measured_mse_weight": 0.0 if shape == CONTROL_ONLY else 1.0,
        "measured_kl_weight": 0.0 if shape == CONTROL_ONLY else 0.1,
        "prediction_kl_weight": 0.0 if shape == CONTROL_ONLY else 0.1,
        "delta_mse_weight": 2.0 if shape == CONTROL_ONLY else 1.0, "delta_l1_weight": 1.75,
    }


def standard_normal_kl(mean, logvar):
    return -0.5 * torch.sum(1. + logvar - mean.square() - logvar.exp())


def measured_objective(model, shape, control, control_rec, mu1, logvar1, measured, treatment_rec,
                       mu2, logvar2, predicted, mu_pred, logvar_pred):
    if any(value is not None for value in (treatment_rec, mu2, logvar2)):
        raise ValueError("unexpected measured treated-expression path")
    delta_hat, delta = predicted - control, measured - control
    control_mse = F.mse_loss(control_rec, control, reduction="sum")
    delta_mse = F.mse_loss(delta_hat, delta, reduction="sum")
    delta_l1 = F.l1_loss(delta_hat, delta, reduction="sum")
    k1 = standard_normal_kl(mu1, logvar1)
    if shape == CONTROL_ONLY:
        measured_mse = k2 = kp = control_mse.new_zeros(())
        loss = 0.1 * k1 + control_mse + 2 * delta_mse + 1.75 * delta_l1
    else:
        delta_rec, mu_delta, logvar_delta = model.reconstruct_measured_delta(delta)
        measured_mse = F.mse_loss(delta_rec, delta, reduction="sum")
        k2 = standard_normal_kl(mu_delta, logvar_delta)
        kp = -0.5 * torch.sum(1 + (logvar_pred - logvar_delta)
                              - ((mu_pred - mu_delta).square() + logvar_pred.exp()) / logvar_delta.exp())
        loss = 0.1 * (k1 + k2 + kp) + control_mse + measured_mse + delta_mse + 1.75 * delta_l1
    return loss, control_mse, measured_mse, delta_mse, k1, k2, kp, delta_l1


def configure_objective(model, shape):
    if shape in SOURCE_SHAPES:
        configure_round_three(model, SOURCE_SHAPES[shape])
    else:
        model.loss = partial(measured_objective, model, shape)
