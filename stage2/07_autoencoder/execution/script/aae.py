"""Full-axis Attention Autoencoder and its frozen Step 07 fit recipe."""

from contextlib import contextmanager
from functools import partial
from math import log

import torch
from torch import nn
from torch.nn import functional as F
from torch.nn.attention import SDPBackend, sdpa_kernel

from benchmark import AAE_MODEL_ID, SEED
from neural import fit_neural, load_neural


LATENT_TOKENS = 512
WIDTH = 128
HEADS = 8
RECIPE = {
    "seed": SEED,
    "epochs": 20,
    "train_batch_size": 32,
    "validation_batch_size": 32,
    "max_seconds": 3600,
    "optimizer": "AdamW",
    "learning_rate": 0.001,
    "weight_decay": 0.01,
    "objective": "full-axis reconstruction MSE",
    "train_order": "PCG64(seed, epoch) shuffle of the joint train membership",
}


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


@contextmanager
def attention_context(device):
    if device == "cuda":
        # Require the bounded fused kernel instead of allocating a gene-by-token matrix.
        with sdpa_kernel(SDPBackend.FLASH_ATTENTION), torch.autocast("cuda", dtype=torch.bfloat16):
            yield
    else:
        yield


def recipe_for(fixture, device):
    if device == "cuda":
        if not torch.cuda.is_bf16_supported(including_emulation=False):
            raise ValueError("AAE requires CUDA with native BF16 support")
        recipe = {**RECIPE, "precision": "bfloat16", "attention_backend": "FlashAttention"}
    elif fixture:
        recipe = {**RECIPE, "precision": "float32", "attention_backend": "PyTorch CPU SDPA"}
    else:
        raise ValueError("full-axis AAE fitting and replay require CUDA")
    if fixture:
        recipe.update(epochs=3, train_batch_size=2, validation_batch_size=2)
    return recipe


def fit_aae(corpus, split, output, batch_size, device):
    recipe = recipe_for(corpus.fixture, device)
    if batch_size is not None and batch_size != recipe["train_batch_size"]:
        raise ValueError(f"AAE fit requires the frozen batch size {recipe['train_batch_size']}")
    return fit_neural(
        corpus, split, output, AAE_MODEL_ID, lambda: TranscriptomeAutoencoder(corpus.gene_count),
        recipe, device,
        {"architecture": architecture(corpus.gene_count), "latent_tokens": LATENT_TOKENS,
         "latent_scalars": LATENT_TOKENS * WIDTH},
        context=partial(attention_context, device), deterministic=device == "cpu",
    )


def aae_predictor(path, record, fixture, batch_size, device):
    recipe = recipe_for(fixture, device)
    replay_recipes = [recipe]
    if not fixture:
        # Keep the originally published four-epoch checkpoint replayable.
        replay_recipes.append({**recipe, "epochs": 4, "max_seconds": 1800})
    metadata = record["metadata"]
    gene_count = record["gene_count"]
    if (metadata["architecture"] != architecture(gene_count)
            or metadata["latent_tokens"] != LATENT_TOKENS
            or metadata["latent_scalars"] != LATENT_TOKENS * WIDTH
            or metadata["recipe"] not in replay_recipes
            or set(record["files"]) != {"weights.pt", "epochs.json"}):
        raise ValueError("invalid AAE artifact or incompatible replay runtime")
    score_batch_size = recipe["validation_batch_size"]
    if batch_size is not None and batch_size != score_batch_size:
        raise ValueError(f"AAE replay requires the frozen batch size {score_batch_size}")
    predict = load_neural(
        path, TranscriptomeAutoencoder(gene_count), device,
        context=partial(attention_context, device), deterministic=device == "cpu",
    )
    return predict, score_batch_size
