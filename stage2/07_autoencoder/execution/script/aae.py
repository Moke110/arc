"""Full-axis Attention Autoencoder and its frozen Step 07 fit recipe."""

from contextlib import contextmanager
from functools import partial
from math import log

import torch
from torch import nn
from torch.nn import functional as F
from torch.nn.attention import SDPBackend, sdpa_kernel

from benchmark import (
    AAE_MODEL_ID, AAE_NUISANCE_MODEL_ID, AAE_TRANSFORMER_MODEL_ID, AAE_TRANSFORMER_MODEL_IDS, SEED,
)
from neural import check_selected_checkpoint, fit_neural, load_neural


LATENT_TOKENS = 512
NUISANCE_TOKENS = 16
WIDTH = 128
HEADS = 8
LATENT_BLOCKS = 2
FEEDFORWARD_WIDTH = 4 * WIDTH
RECIPE = {
    "seed": SEED,
    "max_optimizer_steps": 10_000,
    "validation_patience": 3,
    "train_batch_size": 32,
    "validation_batch_size": 32,
    "max_seconds": 1800,
    "optimizer": "AdamW",
    "learning_rate": 0.001,
    "weight_decay": 0.01,
    "objective": "full-axis reconstruction MSE",
    "train_order": "PCG64(seed, epoch) shuffle of the joint train membership",
}


def architecture(gene_count, model_id=AAE_MODEL_ID):
    result = {
        "gene_count": gene_count, "latent_tokens": LATENT_TOKENS, "width": WIDTH, "heads": HEADS,
    }
    if model_id in AAE_TRANSFORMER_MODEL_IDS:
        result.update(encoder_blocks=LATENT_BLOCKS, feedforward_width=FEEDFORWARD_WIDTH,
                      decoder_tokens=LATENT_TOKENS)
    if model_id == AAE_NUISANCE_MODEL_ID:
        result["latent_tokens"] += NUISANCE_TOKENS
    return result


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


class TransformerTranscriptomeAutoencoder(TranscriptomeAutoencoder):
    def __init__(self, gene_count, *, nuisance_tokens=0):
        super().__init__(gene_count)
        self.encoder = nn.Sequential(*(LatentTransformerBlock() for _ in range(LATENT_BLOCKS)))
        if nuisance_tokens:
            # Draw extra queries last to preserve every shared parameter at the same seed.
            extra_queries = torch.randn(1, nuisance_tokens, WIDTH) / WIDTH**0.5
            self.latent_queries = nn.Parameter(
                torch.cat((self.latent_queries.detach(), extra_queries), dim=1),
            )

    def encode(self, expression):
        return self.encoder(super().encode(expression))

    def decode(self, latent):
        return super().decode(latent[:, :LATENT_TOKENS])


MODEL_CLASSES = {
    AAE_MODEL_ID: TranscriptomeAutoencoder,
    AAE_TRANSFORMER_MODEL_ID: TransformerTranscriptomeAutoencoder,
    AAE_NUISANCE_MODEL_ID: partial(TransformerTranscriptomeAutoencoder, nuisance_tokens=NUISANCE_TOKENS),
}


@contextmanager
def attention_context(device):
    if device == "cuda":
        # Require the bounded fused kernel instead of allocating a gene-by-token matrix.
        with sdpa_kernel(SDPBackend.FLASH_ATTENTION), torch.autocast("cuda", dtype=torch.bfloat16):
            yield
    else:
        yield


def recipe_for(fixture, device, *, model_name="AAE"):
    if device == "cuda":
        if not torch.cuda.is_bf16_supported(including_emulation=False):
            raise ValueError(f"{model_name} requires CUDA with native BF16 support")
        recipe = {**RECIPE, "precision": "bfloat16", "attention_backend": "FlashAttention"}
    elif fixture:
        recipe = {**RECIPE, "precision": "float32", "attention_backend": "PyTorch CPU SDPA"}
    else:
        raise ValueError(f"full-axis {model_name} fitting and replay require CUDA")
    if fixture:
        recipe.update(max_optimizer_steps=5, train_batch_size=2, validation_batch_size=2)
    return recipe


def fit_aae(corpus, split, output, batch_size, device, model_id=AAE_MODEL_ID):
    recipe = recipe_for(corpus.fixture, device)
    model_class = MODEL_CLASSES[model_id]
    design = architecture(corpus.gene_count, model_id)
    latent_tokens = design["latent_tokens"]
    if batch_size is not None and batch_size != recipe["train_batch_size"]:
        raise ValueError(f"AAE fit requires the frozen batch size {recipe['train_batch_size']}")
    return fit_neural(
        corpus, split, output, model_id, lambda: model_class(corpus.gene_count),
        recipe, device,
        {"architecture": design, "latent_tokens": latent_tokens, "latent_scalars": latent_tokens * WIDTH},
        context=partial(attention_context, device), deterministic=device == "cpu",
    )


def aae_predictor(path, record, fixture, batch_size, device):
    recipe = recipe_for(fixture, device)
    replay_recipes = [recipe]
    model_id = record["model_id"]
    if model_id == AAE_MODEL_ID:
        # Keep the originally published fixed-epoch checkpoints replayable.
        legacy_recipe = {**recipe, "epochs": 3 if fixture else 4}
        del legacy_recipe["max_optimizer_steps"], legacy_recipe["validation_patience"]
        replay_recipes.append(legacy_recipe)
    metadata = record["metadata"]
    gene_count = record["gene_count"]
    design = architecture(gene_count, model_id)
    latent_tokens = design["latent_tokens"]
    if (metadata["architecture"] != design
            or metadata["latent_tokens"] != latent_tokens
            or metadata["latent_scalars"] != latent_tokens * WIDTH
            or metadata["recipe"] not in replay_recipes
            or set(record["files"]) != {"weights.pt", "epochs.json"}):
        raise ValueError("invalid AAE artifact or incompatible replay runtime")
    score_batch_size = recipe["validation_batch_size"]
    if batch_size is not None and batch_size != score_batch_size:
        raise ValueError(f"AAE replay requires the frozen batch size {score_batch_size}")
    if model_id in AAE_TRANSFORMER_MODEL_IDS:
        check_selected_checkpoint(path, metadata, model_id, device)
    predict = load_neural(
        path, MODEL_CLASSES[model_id](gene_count), device,
        context=partial(attention_context, device), deterministic=device == "cpu",
    )
    return predict, score_batch_size
