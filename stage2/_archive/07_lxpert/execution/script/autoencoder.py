"""Full-axis deterministic transcriptome autoencoder."""
from math import log

import torch
from torch import nn
from torch.nn import functional as F


LATENT_TOKENS = 512
WIDTH = 128
HEADS = 8


def architecture(gene_count):
    return {"gene_count": gene_count, "latent_tokens": LATENT_TOKENS,
            "width": WIDTH, "heads": HEADS}


class CrossAttention(nn.Module):
    def __init__(self):
        super().__init__()
        self.query = nn.Linear(WIDTH, WIDTH)
        self.key = nn.Linear(WIDTH, WIDTH)
        self.value = nn.Linear(WIDTH, WIDTH)
        self.output = nn.Linear(WIDTH, WIDTH)

    def forward(self, queries, tokens):
        batch = tokens.shape[0]
        q = self.query(queries).unflatten(-1, (HEADS, WIDTH // HEADS)).transpose(1, 2)
        k = self.key(tokens).unflatten(-1, (HEADS, WIDTH // HEADS)).transpose(1, 2)
        v = self.value(tokens).unflatten(-1, (HEADS, WIDTH // HEADS)).transpose(1, 2)
        attended = F.scaled_dot_product_attention(q.expand(batch, -1, -1, -1), k, v)
        return self.output(attended.transpose(1, 2).flatten(-2))


class TranscriptomeAutoencoder(nn.Module):
    def __init__(self, gene_count):
        super().__init__()
        frequencies = torch.exp(torch.arange(0, WIDTH, 2) * (-log(10000.0) / WIDTH))
        angles = torch.arange(gene_count)[:, None] * frequencies
        positions = torch.stack((angles.sin(), angles.cos()), dim=-1).flatten(-2)
        self.register_buffer("gene_positions", positions.unsqueeze(0), persistent=False)
        self.latent_queries = nn.Parameter(torch.randn(1, LATENT_TOKENS, WIDTH) / WIDTH**0.5)
        self.expression = nn.Linear(1, WIDTH)
        self.read_expression = CrossAttention()
        self.read_latent = CrossAttention()
        self.latent_norm = nn.LayerNorm(WIDTH)
        self.reconstruct = nn.Linear(WIDTH, 1)

    def encode(self, expression):
        values = self.expression(expression.unsqueeze(-1))
        tokens = values + self.gene_positions.to(dtype=values.dtype)
        return self.latent_norm(self.read_expression(self.latent_queries, tokens)
                                + self.latent_queries)

    def decode(self, latent):
        attended = self.read_latent(self.gene_positions, latent)
        return self.reconstruct(attended + self.gene_positions.to(dtype=attended.dtype)).squeeze(-1)

    def forward(self, expression):
        return self.decode(self.encode(expression))


def reconstruction_mse(prediction, expression, observed):
    """Return one observed-target MSE per profile, before averaging profiles."""
    error = (prediction.float() - expression).square()
    return torch.where(observed, error, 0).sum(dim=-1) / observed.sum(dim=-1)
