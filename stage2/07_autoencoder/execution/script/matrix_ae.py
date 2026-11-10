"""Full-axis autoencoder with a deterministic 64 by 8 latent matrix."""

from functools import partial

import torch
from torch import nn

from aae import CrossAttention, HEADS, WIDTH, attention_context, recipe_for
from benchmark import MATRIX_AE_CHANNELS, MATRIX_AE_MODEL_ID, MATRIX_AE_TOKENS
from neural import fit_neural, load_neural


LATENT_BLOCKS = 2
FEEDFORWARD_WIDTH = 4 * WIDTH


class LatentTransformerBlock(nn.Module):
    def __init__(self):
        super().__init__()
        self.attention_norm = nn.LayerNorm(WIDTH)
        self.attention = CrossAttention()
        self.feedforward_norm = nn.LayerNorm(WIDTH)
        self.feedforward = nn.Sequential(
            nn.Linear(WIDTH, FEEDFORWARD_WIDTH),
            nn.GELU(),
            nn.Linear(FEEDFORWARD_WIDTH, WIDTH),
        )

    def forward(self, tokens):
        normalized = self.attention_norm(tokens)
        tokens = tokens + self.attention(normalized, normalized)
        return tokens + self.feedforward(self.feedforward_norm(tokens))


class MatrixAutoencoder(nn.Module):
    def __init__(self, gene_count):
        super().__init__()
        self.gene_identities = nn.Parameter(torch.randn(1, gene_count, WIDTH) / WIDTH**0.5)
        self.latent_queries = nn.Parameter(torch.randn(1, MATRIX_AE_TOKENS, WIDTH) / WIDTH**0.5)
        self.expression = nn.Linear(1, WIDTH)
        self.read_expression = CrossAttention()
        self.encoder = nn.Sequential(*(LatentTransformerBlock() for _ in range(LATENT_BLOCKS)))
        self.compress = nn.Linear(WIDTH, MATRIX_AE_CHANNELS)
        self.expand = nn.Linear(MATRIX_AE_CHANNELS, WIDTH)
        self.decoder = nn.Sequential(*(LatentTransformerBlock() for _ in range(LATENT_BLOCKS)))
        self.read_latent = CrossAttention()
        self.reconstruct = nn.Linear(WIDTH, 1)

    def encode(self, expression):
        values = self.expression(expression.unsqueeze(-1))
        tokens = values + self.gene_identities.to(dtype=values.dtype)
        pooled = self.read_expression(self.latent_queries, tokens)
        pooled = pooled + self.latent_queries.to(dtype=pooled.dtype)
        return self.compress(self.encoder(pooled))

    def decode(self, latent):
        tokens = self.decoder(self.expand(latent))
        attended = self.read_latent(self.gene_identities, tokens)
        return self.reconstruct(attended + self.gene_identities.to(dtype=attended.dtype)).squeeze(-1)

    def forward(self, expression):
        return self.decode(self.encode(expression))


def architecture(gene_count):
    return {
        "gene_count": gene_count,
        "width": WIDTH,
        "heads": HEADS,
        "latent_tokens": MATRIX_AE_TOKENS,
        "latent_channels": MATRIX_AE_CHANNELS,
        "encoder_blocks": LATENT_BLOCKS,
        "decoder_blocks": LATENT_BLOCKS,
        "feedforward_width": FEEDFORWARD_WIDTH,
    }


def fit_matrix_ae(corpus, split, output, batch_size, device):
    recipe = recipe_for(corpus.fixture, device, model_name=MATRIX_AE_MODEL_ID)
    if batch_size is not None and batch_size != recipe["train_batch_size"]:
        raise ValueError(f"{MATRIX_AE_MODEL_ID} fit requires the frozen batch size "
                         f"{recipe['train_batch_size']}")
    return fit_neural(
        corpus, split, output, MATRIX_AE_MODEL_ID, lambda: MatrixAutoencoder(corpus.gene_count),
        recipe, device,
        {"architecture": architecture(corpus.gene_count), "latent_tokens": MATRIX_AE_TOKENS,
         "latent_scalars": MATRIX_AE_TOKENS * MATRIX_AE_CHANNELS},
        context=partial(attention_context, device), deterministic=device == "cpu",
    )


def matrix_ae_predictor(path, record, fixture, batch_size, device):
    recipe = recipe_for(fixture, device, model_name=MATRIX_AE_MODEL_ID)
    metadata = record["metadata"]
    gene_count = record["gene_count"]
    if (metadata["architecture"] != architecture(gene_count)
            or metadata["latent_tokens"] != MATRIX_AE_TOKENS
            or metadata["latent_scalars"] != MATRIX_AE_TOKENS * MATRIX_AE_CHANNELS
            or metadata["recipe"] != recipe
            or set(record["files"]) != {"weights.pt", "epochs.json"}):
        raise ValueError(f"invalid {MATRIX_AE_MODEL_ID} artifact or incompatible replay runtime")
    score_batch_size = recipe["validation_batch_size"]
    if batch_size is not None and batch_size != score_batch_size:
        raise ValueError(f"{MATRIX_AE_MODEL_ID} replay requires the frozen batch size "
                         f"{score_batch_size}")
    predict = load_neural(
        path, MatrixAutoencoder(gene_count), device,
        context=partial(attention_context, device), deterministic=device == "cpu",
    )
    return predict, score_batch_size
