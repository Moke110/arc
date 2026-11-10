"""Channel-shared token mixing for full-axis matrix MLP reconstruction."""

from contextlib import nullcontext
from functools import partial

import torch
from torch import nn

from aae import GenePositionedAutoencoder, LATENT_TOKENS, WIDTH, recipe_for as aae_recipe_for
from benchmark import MATRIX_MLP_MODEL_ID
from neural import check_selected_checkpoint, fit_neural, load_neural


TOKEN_MLP_HIDDEN = 128


def token_mlp(input_tokens, output_tokens):
    return nn.Sequential(
        nn.Linear(input_tokens, TOKEN_MLP_HIDDEN),
        nn.GELU(),
        nn.Linear(TOKEN_MLP_HIDDEN, output_tokens),
    )


class MatrixMLPAutoencoder(GenePositionedAutoencoder):
    def __init__(self, gene_count):
        super().__init__(gene_count)
        self.expression = nn.Linear(1, WIDTH)
        self.encode_tokens = token_mlp(gene_count, LATENT_TOKENS)
        self.decode_tokens = token_mlp(LATENT_TOKENS, gene_count)
        self.latent_norm = nn.LayerNorm(WIDTH)
        self.reconstruct = nn.Linear(WIDTH, 1)

    def encode(self, expression):
        tokens = self.embed_expression(expression)
        # Linear maps the final axis and shares weights over profiles and channels.
        mixed = self.encode_tokens(tokens.transpose(1, 2).contiguous())
        return self.latent_norm(mixed.transpose(1, 2))

    def decode(self, latent):
        tokens = self.decode_tokens(latent.transpose(1, 2).contiguous())
        return self.read_expression_tokens(tokens.transpose(1, 2).contiguous())


def design_metadata(gene_count):
    with torch.device("meta"):
        model = MatrixMLPAutoencoder(gene_count)
    return {
        "architecture": {
            "gene_count": gene_count, "latent_tokens": LATENT_TOKENS, "width": WIDTH,
            "token_mlp_hidden": TOKEN_MLP_HIDDEN, "token_mlp_activation": "GELU",
        },
        "latent_tokens": LATENT_TOKENS,
        "latent_scalars": LATENT_TOKENS * WIDTH,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
    }


def recipe_for(fixture, device, *, check_device=True):
    recipe = aae_recipe_for(fixture, device, model_name=MATRIX_MLP_MODEL_ID, check_device=check_device)
    del recipe["attention_backend"]
    return recipe


def precision_context(device):
    if device == "cuda":
        return partial(torch.autocast, "cuda", dtype=torch.bfloat16)
    return nullcontext


def fit_matrix_mlp(corpus, split, output, batch_size, device):
    recipe = recipe_for(corpus.fixture, device)
    if batch_size is not None and batch_size != recipe["train_batch_size"]:
        raise ValueError(f"{MATRIX_MLP_MODEL_ID} fit requires the frozen batch size "
                         f"{recipe['train_batch_size']}")
    return fit_neural(
        corpus, split, output, MATRIX_MLP_MODEL_ID, lambda: MatrixMLPAutoencoder(corpus.gene_count),
        recipe, device, design_metadata(corpus.gene_count),
        context=precision_context(device), deterministic=device == "cpu",
    )


def check_matrix_mlp_artifact(path, record, fixture, device, *, check_device=False):
    recipe = recipe_for(fixture, device, check_device=check_device)
    metadata = record["metadata"]
    expected = design_metadata(record["gene_count"])
    if (any(metadata.get(field) != value for field, value in expected.items())
            or metadata["recipe"] != recipe
            or set(record["files"]) != {"weights.pt", "epochs.json"}):
        raise ValueError(f"{MATRIX_MLP_MODEL_ID}: invalid artifact or incompatible replay runtime")
    check_selected_checkpoint(path, metadata, MATRIX_MLP_MODEL_ID, device)
    return recipe["validation_batch_size"]


def matrix_mlp_predictor(path, record, fixture, batch_size, device):
    score_batch_size = check_matrix_mlp_artifact(path, record, fixture, device, check_device=True)
    if batch_size is not None and batch_size != score_batch_size:
        raise ValueError(f"{MATRIX_MLP_MODEL_ID} replay requires the frozen batch size {score_batch_size}")
    predict = load_neural(
        path, MatrixMLPAutoencoder(record["gene_count"]), device,
        context=precision_context(device), deterministic=device == "cpu",
    )
    return predict, score_batch_size
