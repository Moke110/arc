"""TranSiGen prediction architectures on the fixed Step 06.1 expression shape."""

from fractions import Fraction
from functools import partial
from pathlib import Path
import sys

import torch
from torch import nn
from torch.nn import functional as F

STEP = Path(__file__).resolve().parents[2] / "06_tsg_shape/06.2_round2"
ROUND_ONE = STEP.with_name("06.1_round1")
sys.path.insert(0, str(ROUND_ONE / "execution"))

from . import native
from .shape_model import ShapeAdapter
from .additive_objective import AdditiveObjective
from .architecture_roster import BOTH, CONTROL, DIRECT_DELTA, FILM, MLP, NO_X2, RESIDUAL, SHARED, SHAPES

BASE_OBJECTIVE = AdditiveObjective(Fraction("1.75"))


class FiLMPrediction(nn.Module):
    def __init__(self):
        super().__init__()
        self.conditioner = nn.Sequential(nn.Linear(432, 136), nn.ReLU(), nn.Linear(136, 200))
        nn.init.zeros_(self.conditioner[-1].weight)
        nn.init.zeros_(self.conditioner[-1].bias)
        self.projection = nn.Linear(100, 200)

    def forward(self, control, drug_dose):
        scale, shift = self.conditioner(drug_dose).chunk(2, dim=1)
        return self.projection((1 + scale) * control + shift)


class ArchitectureAdapter(ShapeAdapter):
    def network(self):
        if self.shape not in SHAPES:
            raise ValueError(f"unregistered Step 06.2 architecture: {self.shape}")
        original = super().network() if self.shape == CONTROL else native.Adapter.network(self)
        if self.shape == CONTROL:
            return original
        shape = self.shape

        class Architecture(original):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                if shape in (SHARED, BOTH):
                    self.encoder_x2 = self.encoder_x1
                if shape in (MLP, FILM):
                    del self.mu_z2Fz1, self.logvar_z2Fz1
                    self.prediction_head = (FiLMPrediction() if shape == FILM else
                                            nn.Sequential(nn.Linear(532, 145), nn.ReLU(), nn.Linear(145, 200)))
                if shape == DIRECT_DELTA:
                    hidden = self.n_de_hidden[0]
                    self.decoder_delta = nn.Sequential(
                        nn.Linear(self.n_latent, hidden), nn.BatchNorm1d(hidden), nn.ReLU(),
                        nn.Dropout(self.dropout), nn.Linear(hidden, self.n_genes),
                    )
                    self.decoder_delta.apply(self._init_weights)
                if shape == NO_X2:
                    self.decoder_x2 = self.decoder_x2[:-1]

            def encode_x2(self, measured):
                latent, mean, logvar = super().encode_x2(measured)
                # Preserve the native measured-encoder and sampling order.
                return None if shape == NO_X2 else latent, mean, logvar

            def decode_x2(self, latent):
                # The shared training/validation loop has no reconstruction
                # input for this arm, so the treated decoder is never run here.
                return None if latent is None else super().decode_x2(latent)

            def predictive_distribution(self, z1, drug_dose):
                if shape == FILM:
                    return self.prediction_head(z1, drug_dose).chunk(2, dim=1)
                combined = torch.cat((z1, drug_dose), dim=1)
                if shape == MLP:
                    return self.prediction_head(combined).chunk(2, dim=1)
                mu = self.mu_z2Fz1(combined)
                if shape in (RESIDUAL, BOTH):
                    mu = mu + z1
                return mu, self.logvar_z2Fz1(combined)

            def forward(self, x1, features):
                width = native.WIDTH
                if features.ndim != 2 or features.shape[1] != width + 1:
                    raise ValueError("dose transport requires KPGT2304 plus one scalar")
                z1, mu1, logvar1 = self.encode_x1(x1)
                x1_rec = self.decode_x1(z1)
                drug_dose = torch.cat((self.feat_embeddings(features[:, :width]),
                                      features[:, width:].repeat(1, 32)), dim=1)
                mu_pred, logvar_pred = self.predictive_distribution(z1, drug_dose)
                z_pred = self.sample_latent(mu_pred, logvar_pred)
                if shape == DIRECT_DELTA:
                    predicted = x1 + self.decoder_delta(z_pred)
                elif shape == NO_X2:
                    predicted = x1 + self.decoder_x2(z_pred)
                else:
                    predicted = self.decode_x2(z_pred)
                return x1_rec, mu1, logvar1, predicted, mu_pred, logvar_pred, z_pred

        return Architecture


def objective_identity(shape):
    identity = BASE_OBJECTIVE.identity()
    if shape == NO_X2:
        del identity["reconstruction_weight"]
        identity.update(control_mse_weight=1., treatment_mse_weight=0., delta_mse_weight=2.)
    return identity


def no_treatment_objective(control, control_rec, mu1, logvar1, measured, treatment_rec,
                           mu2, logvar2, predicted, mu_pred, logvar_pred):
    """Keep native KL components with an omitted treatment autoencoder MSE."""
    if treatment_rec is not None:
        raise ValueError("output-delta-nox2 has no treatment reconstruction")
    control_mse = F.mse_loss(control_rec, control, reduction="sum")
    delta_hat, delta = predicted - control, measured - control
    delta_mse = F.mse_loss(delta_hat, delta, reduction="sum")
    delta_l1 = F.l1_loss(delta_hat, delta, reduction="sum")
    k1 = -0.5 * torch.sum(1. + logvar1 - mu1.pow(2) - logvar1.exp())
    k2 = -0.5 * torch.sum(1. + logvar2 - mu2.pow(2) - logvar2.exp())
    kp = -0.5 * torch.sum(1 + (logvar_pred - logvar2)
                          - ((mu_pred - mu2).pow(2) + logvar_pred.exp()) / logvar2.exp())
    loss = 0.1 * (k1 + k2 + kp) + control_mse + 2 * delta_mse + 1.75 * delta_l1
    return loss, control_mse, control_mse.new_zeros(()), delta_mse, k1, k2, kp, delta_l1


def configure_objective(model, shape):
    model.loss = (no_treatment_objective if shape == NO_X2 else partial(BASE_OBJECTIVE, model.loss))
