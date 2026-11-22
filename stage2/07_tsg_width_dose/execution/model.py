"""Width adapters for the immutable native TranSiGen and Step 04 dose rules."""

import torch
from torch import nn

from run import METHODS, WIDTHS
from tsg_kpgt import model as step04
from tsg_kpgt.inputs import WIDTH as KPGT_WIDTH


def dimensions(width):
    numerator = WIDTHS[width]
    return {name: base * numerator // 2 for name, base in
            dict(latent=100, encoder=1200, decoder=800, chemical=400, dose=32).items()}


class ExpandedWidth:
    def build(self, seed, device):
        torch.manual_seed(seed)
        self.device = torch.device(device)
        dims = self.dims
        model = self.network()(
            n_genes=5000,
            n_latent=dims['latent'],
            n_en_hidden=[dims['encoder']],
            n_de_hidden=[dims['decoder']],
            features_dim=KPGT_WIDTH,
            features_embed_dim=[dims['chemical']],
            init_w=True,
            beta=0.1,
            dropout=0.1,
            random_seed=seed,
            device=self.device,
        )
        return model.to(device=self.device, dtype=torch.float32)


class ExpandedNoDose(ExpandedWidth, step04.NoDoseAdapter):
    pass


class ExpandedConditioned(ExpandedWidth, step04.ConditionedAdapter):
    def network(self):
        method = self.method
        latent, chemical, dose = (self.dims[key] for key in ('latent', 'chemical', 'dose'))

        class Conditioned(self.native):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                if method == 'film':
                    self.film = nn.Sequential(nn.Linear(1, dose), nn.ReLU(),
                                              nn.Linear(dose, 2 * chemical))
                    nn.init.zeros_(self.film[-1].weight)
                    nn.init.zeros_(self.film[-1].bias)
                else:
                    size = latent + chemical + dose
                    self.mu_z2Fz1 = nn.Sequential(nn.Linear(size, latent))
                    self.logvar_z2Fz1 = nn.Sequential(nn.Linear(size, latent))
                    if method == 'latent-embedding':
                        self.dose_embedding = nn.Sequential(
                            nn.Linear(1, dose), nn.ReLU(), nn.Linear(dose, dose)
                        )

            def forward(self, x1, features):
                if features.ndim != 2 or features.shape[1] != KPGT_WIDTH + 1:
                    raise ValueError('Dose transport requires KPGT2304 and one scalar')
                kpgt, scalar = features[:, :KPGT_WIDTH], features[:, KPGT_WIDTH:]
                z1, mu1, logvar1 = self.encode_x1(x1)
                x1_rec = self.decode_x1(z1)
                embedding = self.feat_embeddings(kpgt)
                if method == 'film':
                    gamma, beta = self.film(scalar).chunk(2, dim=1)
                    conditioned = torch.cat([z1, (1 + gamma) * embedding + beta], dim=1)
                else:
                    dose_features = (self.dose_embedding(scalar) if method == 'latent-embedding'
                                     else scalar.repeat(1, dose))
                    conditioned = torch.cat([z1, embedding, dose_features], dim=1)
                mu_pred = self.mu_z2Fz1(conditioned)
                logvar_pred = self.logvar_z2Fz1(conditioned)
                z2_pred = self.sample_latent(mu_pred, logvar_pred)
                return x1_rec, mu1, logvar1, self.decode_x2(z2_pred), mu_pred, logvar_pred, z2_pred

        return Conditioned


def adapter(source, width, method, transform):
    dose_method = METHODS[method]
    if width == 'original':
        return step04.adapter(source, dose_method, transform)
    step04.configuration(dose_method, transform)
    native = (ExpandedNoDose(source) if dose_method == 'no-dose'
              else ExpandedConditioned(source, dose_method, transform))
    native.dims = dimensions(width)
    return native
