"""Native TranSiGen with the accepted dose transforms and a KPGT input layer."""

import torch
from torch import nn

from tsg_dose.model import (
    Adapter,
    DoseAdapter,
    ScalingAdapter,
    configure_backend,
    fit_transform,
    SOURCE_DIGEST,
)
from tsg_round1.model import configuration as round1_configuration

from .inputs import WIDTH


def configuration(method, transform):
    base = round1_configuration(method, transform)
    adaptations = [
        value
        for value in base['adaptations']
        if 'Morgan2048' not in value
        and 'no pretraining' not in value
        and value != 'unstandardized fingerprint scaling'
    ]
    value = dict(
        base,
        method_identity='tsg-tenv-kpgt-' + method + '-f32-v1',
        molecule_feature='KPGT',
        representation_width=WIDTH,
        features_dim=WIDTH + (32 if method.startswith('input-') else 0),
        adaptations=adaptations
        + ['frozen KPGT2304 input', 'random response-model initialization'],
    )
    if method == 'scaling':
        value['dose_scaling'] = 'KPGT × log10(1+dose/1µM); unstandardized'
    return value


class NativeWidth:
    def build(self, seed, device):
        torch.manual_seed(seed)
        self.device = torch.device(device)
        model = self.network()(
            n_genes=5000,
            n_latent=100,
            n_en_hidden=[1200],
            n_de_hidden=[800],
            features_dim=self.input_width,
            features_embed_dim=[400],
            init_w=True,
            beta=0.1,
            dropout=0.1,
            random_seed=seed,
            device=self.device,
        )
        return model.to(device=self.device, dtype=torch.float32)

    def network(self):
        return self.native

    input_width = WIDTH


class NoDoseAdapter(NativeWidth, Adapter):
    pass


class ScaledAdapter(NativeWidth, ScalingAdapter):
    pass


class ConditionedAdapter(NativeWidth, DoseAdapter):
    @property
    def input_width(self):
        return WIDTH + (32 if self.method.startswith('input-') else 0)

    def network(self):
        method = self.method

        class Conditioned(self.native):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                if method.startswith('latent-'):
                    self.mu_z2Fz1 = nn.Sequential(nn.Linear(532, 100))
                    self.logvar_z2Fz1 = nn.Sequential(nn.Linear(532, 100))
                if method.endswith('embedding'):
                    self.dose_embedding = nn.Sequential(
                        nn.Linear(1, 32), nn.ReLU(), nn.Linear(32, 32)
                    )
                if method == 'film':
                    self.film = nn.Sequential(
                        nn.Linear(1, 32), nn.ReLU(), nn.Linear(32, 800)
                    )
                    nn.init.zeros_(self.film[-1].weight)
                    nn.init.zeros_(self.film[-1].bias)

            def forward(self, x1, features):
                if features.ndim != 2 or features.shape[1] != WIDTH + 1:
                    raise ValueError('dose transport requires KPGT2304 plus one scalar')
                kpgt, scalar = features[:, :WIDTH], features[:, WIDTH:]
                if method != 'film':
                    dose = (
                        self.dose_embedding(scalar)
                        if method.endswith('embedding')
                        else scalar.repeat(1, 32)
                    )
                    if method.startswith('input-'):
                        return super().forward(x1, torch.cat([kpgt, dose], dim=1))
                z1, mu1, logvar1 = self.encode_x1(x1)
                x1_rec = self.decode_x1(z1)
                chemical = self.feat_embeddings(kpgt)
                if method == 'film':
                    gamma, beta = self.film(scalar).chunk(2, dim=1)
                    conditioned = torch.cat([z1, (1 + gamma) * chemical + beta], dim=1)
                else:
                    conditioned = torch.cat([z1, chemical, dose], dim=1)
                mu_pred = self.mu_z2Fz1(conditioned)
                logvar_pred = self.logvar_z2Fz1(conditioned)
                z2_pred = self.sample_latent(mu_pred, logvar_pred)
                return (
                    x1_rec,
                    mu1,
                    logvar1,
                    self.decode_x2(z2_pred),
                    mu_pred,
                    logvar_pred,
                    z2_pred,
                )

        return Conditioned


def adapter(source, method, transform):
    configuration(method, transform)
    if method == 'no-dose':
        return NoDoseAdapter(source)
    if method == 'scaling':
        return ScaledAdapter(source)
    return ConditionedAdapter(source, method, transform)
