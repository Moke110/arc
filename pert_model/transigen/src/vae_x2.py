from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class VAE_x2(torch.nn.Module):
    def __init__(self, n_genes, n_latent, n_en_hidden, n_de_hidden, BCE_pos_weight=None, **kwargs):
        super().__init__()
        self.n_genes = n_genes
        self.n_latent = n_latent
        self.n_en_hidden = n_en_hidden
        self.n_de_hidden = n_de_hidden
        self.init_w = kwargs.get("init_w", False)
        self.beta = kwargs.get("beta", 0.05)
        self.dev = kwargs.get("device", torch.device("cpu"))
        self.dropout = kwargs.get("dropout", 0.3)
        self.BCE_pos_weight = torch.tensor([BCE_pos_weight], device=self.dev) if BCE_pos_weight is not None else None

        encoder = [nn.Linear(n_genes, n_en_hidden[0]), nn.BatchNorm1d(n_en_hidden[0]), nn.ReLU(), nn.Dropout(self.dropout)]
        for i in range(len(n_en_hidden) - 1):
            encoder.extend([nn.Linear(n_en_hidden[i], n_en_hidden[i + 1]), nn.BatchNorm1d(n_en_hidden[i + 1]), nn.ReLU(), nn.Dropout(self.dropout)])
        self.encoder_x2 = nn.Sequential(*encoder)
        self.mu_z2 = nn.Linear(n_en_hidden[-1], n_latent)
        self.logvar_z2 = nn.Linear(n_en_hidden[-1], n_latent)

        if len(n_de_hidden) == 0:
            decoder = [nn.Linear(n_latent, n_genes)]
        else:
            decoder = [nn.Linear(n_latent, n_de_hidden[0]), nn.BatchNorm1d(n_de_hidden[0]), nn.ReLU(), nn.Dropout(self.dropout)]
            for i in range(len(n_de_hidden) - 1):
                decoder.extend([nn.Linear(n_de_hidden[i], n_de_hidden[i + 1]), nn.BatchNorm1d(n_de_hidden[i + 1]), nn.ReLU(), nn.Dropout(self.dropout)])
            decoder.extend([nn.Linear(n_de_hidden[-1], n_genes), nn.ReLU()])
        self.decoder_x2 = nn.Sequential(*decoder)

        if self.init_w:
            self.encoder_x2.apply(self._init_weights)
            self.decoder_x2.apply(self._init_weights)

    @staticmethod
    def _init_weights(layer):
        if isinstance(layer, nn.Linear):
            torch.nn.init.xavier_uniform_(layer.weight)

    def encode(self, x):
        y = self.encoder_x2(x)
        mu, logvar = self.mu_z2(y), self.logvar_z2(y)
        return self.sample_latent(mu, logvar), mu, logvar

    def decode(self, z):
        return self.decoder_x2(z)

    def sample_latent(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std, device=std.device)
        return eps * std + mu

    def to_latent(self, x):
        z, _, _ = self.encode(x)
        return z

    def to_mu_logvar(self, x):
        y = self.encoder_x2(x)
        return self.mu_z2(y), self.logvar_z2(y)

    def forward(self, x2):
        z2, z2_mu, z2_logvar = self.encode(x2)
        rec_x2 = self.decode(z2)
        return rec_x2, z2_mu, z2_logvar

    def vae_loss(self, y_pred, y_true, mu, logvar):
        kld = -0.5 * torch.sum(1.0 + logvar - mu.pow(2) - logvar.exp())
        mse = F.mse_loss(y_pred, y_true, reduction="sum")
        return mse + self.beta * kld, mse, kld
