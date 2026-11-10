"""Neural token maps for the frozen Step 07 expression interface."""

from functools import partial
import json
from math import ceil
from pathlib import Path
import time

import torch
from torch import nn

from aae import RECIPE, attention_context
from benchmark import SOURCES, digest
from neural import (
    check_selected_checkpoint, configure_runtime, fit_neural, joint_training_rows, load_neural,
    neural_predictor, train_epoch, training_windows,
)
from shared_expression import load_interface


NEURAL_IDS = ("MLP_1024", "AAEc_1024", "AAEcs_1024", "AAEcs_1024+64")
LATENT_TOKENS = 1024
EXTRA_TOKENS = 64
HEADS = 8
MAX_PROFILES_PER_BATCH = 32
BUDGET_FIELDS = {"train_batch_size", "validation_batch_size", "max_optimizer_steps"}
FIXTURE_BUDGET = {"train_batch_size": 2, "validation_batch_size": 2, "max_optimizer_steps": 5}


def token_shape(model_id, gene_count, fixture):
    tokens = min(3, gene_count) if fixture else LATENT_TOKENS
    extras = (2 if fixture else EXTRA_TOKENS) if model_id == "AAEcs_1024+64" else 0
    return tokens, extras


def token_map(input_tokens, output_tokens, hidden_tokens):
    return nn.Sequential(
        nn.Linear(input_tokens, hidden_tokens), nn.GELU(),
        nn.Linear(hidden_tokens, output_tokens),
    )


class TokenMLP(nn.Module):
    def __init__(self, gene_count, tokens):
        super().__init__()
        self.encoder = token_map(gene_count, tokens, tokens)
        self.decoder = token_map(tokens, gene_count, tokens)

    def encode(self, matrix):
        return self.encoder(matrix.transpose(1, 2).contiguous()).transpose(1, 2)

    def decode(self, latent, positions):
        matrix = self.decoder(latent.transpose(1, 2).contiguous())
        return matrix.transpose(1, 2).contiguous() + positions


class LatentAttention(nn.Module):
    def __init__(self, channels):
        super().__init__()
        self.norm = nn.LayerNorm(channels)
        self.attention = nn.MultiheadAttention(channels, HEADS, batch_first=True)

    def forward(self, latent):
        normalized = self.norm(latent)
        return latent + self.attention(normalized, normalized, normalized, need_weights=False)[0]


class MatrixAttention(nn.Module):
    def __init__(self, channels, tokens, extras, self_attention):
        super().__init__()
        self.reconstruction_tokens = tokens
        self.latent_queries = nn.Parameter(torch.randn(1, tokens, channels) / channels**0.5)
        self.encoder_cross = nn.MultiheadAttention(channels, HEADS, batch_first=True)
        self.decoder_cross = nn.MultiheadAttention(channels, HEADS, batch_first=True)
        self.encoder_norm = nn.LayerNorm(channels)
        self.encoder_self = LatentAttention(channels) if self_attention else nn.Identity()
        self.decoder_self = LatentAttention(channels) if self_attention else nn.Identity()
        if extras:
            # Draw the extra queries after every parameter shared by the pair.
            extra_queries = torch.randn(1, extras, channels) / channels**0.5
            self.latent_queries = nn.Parameter(
                torch.cat((self.latent_queries.detach(), extra_queries), dim=1),
            )

    def encode(self, matrix):
        queries = self.latent_queries.expand(matrix.shape[0], -1, -1)
        attended = self.encoder_cross(queries, matrix, matrix, need_weights=False)[0]
        latent = self.encoder_self(self.encoder_norm(attended + queries))
        return latent[:, :self.reconstruction_tokens].contiguous()

    def decode(self, latent, positions):
        if latent.shape[1] != self.reconstruction_tokens:
            raise ValueError("decoder requires the fixed reconstruction token count")
        latent = self.decoder_self(latent)
        queries = positions.unsqueeze(0).expand(latent.shape[0], -1, -1)
        return self.decoder_cross(queries, latent, latent, need_weights=False)[0] + positions


def build_matrix_model(model_id, gene_count, channels, *, fixture):
    if model_id not in NEURAL_IDS:
        raise ValueError(f"unsupported shared neural model {model_id}")
    tokens, extras = token_shape(model_id, gene_count, fixture)
    if model_id == "MLP_1024":
        return TokenMLP(gene_count, tokens)
    return MatrixAttention(channels, tokens, extras, self_attention=model_id != "AAEc_1024")


class SharedReconstructor(nn.Module):
    def __init__(self, matrix_model, interface):
        super().__init__()
        self.matrix_model = matrix_model
        # The verified interface lives in its own artifact. Keep its methods
        # outside the checkpoint and optimizer, and retain the lookup as a buffer.
        self.generate = interface.generate
        self.restore = interface.restore
        self.register_buffer("positions", interface.positions, persistent=False)

    def forward(self, expression):
        # Use the same frozen float32 interface as PCA and the round-trip diagnostic.
        with torch.autocast(expression.device.type, enabled=False):
            matrix = self.generate(expression)
        reconstructed = self.matrix_model.decode(self.matrix_model.encode(matrix), self.positions)
        with torch.autocast(expression.device.type, enabled=False):
            return self.restore(reconstructed.float())


def read_budget(path, fixture):
    if path is None:
        if fixture:
            return FIXTURE_BUDGET.copy()
        raise ValueError("shared neural fitting requires a common --budget")
    return json.loads(Path(path).read_text())


def recipe_for(budget, fixture, device):
    if (not isinstance(budget, dict) or set(budget) != BUDGET_FIELDS
            or any(type(value) is not int or value <= 0 for value in budget.values())
            or max(budget["train_batch_size"], budget["validation_batch_size"]) > MAX_PROFILES_PER_BATCH):
        raise ValueError("invalid common shared neural budget")
    if device != "cuda" and not fixture:
        raise ValueError("real shared neural fitting and replay require CUDA")
    return {
        **RECIPE, **budget,
        "max_seconds": 14_400,
        "validation_patience": budget["max_optimizer_steps"] + 1,
        "validation_interval_steps": ceil(budget["max_optimizer_steps"] / 3),
        "precision": "float32 interface, CUDA BF16 token maps" if device == "cuda" else "float32",
        "attention_backend": "FlashAttention" if device == "cuda" else "PyTorch CPU SDPA",
    }


def design_metadata(model_id, gene_count, channels, fixture):
    tokens, extras = token_shape(model_id, gene_count, fixture)
    with torch.device("meta"):
        model = build_matrix_model(model_id, gene_count, channels, fixture=fixture)
    return {
        "latent_shape": [tokens, channels],
        "encoder_shape": [tokens + extras, channels],
        "trainable_parameters": sum(parameter.numel() for parameter in model.parameters()),
    }


def smoke_shared_neural(corpus, split, model_id, build_model, recipe, device, context):
    train_rows = joint_training_rows(split)
    batch_size = recipe["train_batch_size"]
    steps = recipe["max_optimizer_steps"]
    if steps > min(64, ceil(len(train_rows) / batch_size)):
        raise ValueError("shared neural resource smoke requires at most 64 steps within one epoch")
    configure_runtime(device, deterministic=device == "cpu")
    torch.manual_seed(recipe["seed"])
    if device == "cuda":
        torch.cuda.reset_peak_memory_stats()
    model = build_model().to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=recipe["learning_rate"],
                                  weight_decay=recipe["weight_decay"])
    start = time.perf_counter()
    record = train_epoch(corpus, model, optimizer, train_rows, recipe, 1, device, context,
                         start + recipe["max_seconds"], batch_limit=steps)
    train_seconds = time.perf_counter() - start
    predict = neural_predictor(model, device, context)
    validation_profiles = {}
    for source in SOURCES:
        rows = split.rows[source]["validation"][:recipe["validation_batch_size"]]
        values = corpus.read(source, rows)
        prediction = torch.from_numpy(predict(values))
        if prediction.shape != values.shape or not torch.isfinite(prediction).all():
            raise ValueError(f"{model_id} smoke: invalid {source} prediction")
        validation_profiles[source] = len(rows)
    return {
        **record,
        "model_id": model_id,
        "smoke": True,
        "train_batch_size": batch_size,
        "validation_batch_size": recipe["validation_batch_size"],
        "validation_profiles": validation_profiles,
        "optimizer_steps": steps,
        "train_seconds": train_seconds,
        "optimizer_steps_per_second": steps / train_seconds,
        "train_profiles_per_second": min(len(train_rows), steps * batch_size) / train_seconds,
        "peak_cuda_allocated_bytes": torch.cuda.max_memory_allocated() if device == "cuda" else None,
        "peak_cuda_reserved_bytes": torch.cuda.max_memory_reserved() if device == "cuda" else None,
    }


def fit_shared_neural(corpus, split, output, model_id, device, interface_path, budget_path, smoke=False):
    budget = read_budget(budget_path, corpus.fixture)
    recipe = recipe_for(budget, corpus.fixture, device)
    interface_record, interface = load_interface(interface_path, corpus, split, device)
    if device == "cuda" and not torch.cuda.is_bf16_supported(including_emulation=False):
        raise ValueError("shared neural models require CUDA with native BF16 support")
    design = design_metadata(model_id, corpus.gene_count, interface.channels, corpus.fixture)

    def build_model():
        return SharedReconstructor(
            build_matrix_model(model_id, corpus.gene_count, interface.channels, fixture=corpus.fixture),
            interface,
        )

    context = partial(attention_context, device)
    if smoke:
        return smoke_shared_neural(corpus, split, model_id, build_model, recipe, device, context)
    return fit_neural(
        corpus, split, output, model_id, build_model, recipe, device,
        {**design, "interface_manifest_sha256": interface_record["artifact_manifest_sha256"],
         "budget_sha256": digest(budget)},
        context=context, deterministic=device == "cpu",
    )


def check_shared_neural_artifact(path, record, fixture, device):
    model_id = record["model_id"]
    metadata = record["metadata"]
    recipe = metadata["recipe"]
    budget = {key: recipe[key] for key in BUDGET_FIELDS}
    design = design_metadata(model_id, record["gene_count"], 16 if fixture else 256, fixture)
    if (recipe != recipe_for(budget, fixture, device)
            or metadata["budget_sha256"] != digest(budget)
            or any(metadata.get(field) != value for field, value in design.items())
            or set(record["files"]) != {"weights.pt", "epochs.json"}):
        raise ValueError(f"{model_id}: incompatible shared neural artifact")
    check_selected_checkpoint(path, metadata, model_id, device)
    history = json.loads((Path(path) / "epochs.json").read_text())
    training_profiles = sum(metadata["fit_profiles"][source] for source in SOURCES)
    batch_size = recipe["train_batch_size"]
    windows = training_windows(ceil(training_profiles / batch_size), recipe["max_optimizer_steps"],
                               recipe["validation_interval_steps"])
    steps = 0
    for check, (entry, (epoch, offset, size)) in enumerate(zip(history, windows, strict=True), 1):
        steps += size
        expected = (check, epoch, offset, steps,
                    min(training_profiles - offset * batch_size, size * batch_size))
        actual = tuple(entry[key] for key in (
            "epoch", "training_epoch", "batch_offset", "optimizer_steps", "training_profiles",
        ))
        if actual != expected:
            raise ValueError(f"{model_id}: training history differs from the common budget")
    if steps != recipe["max_optimizer_steps"] or metadata["stop_reason"] != "optimizer step limit":
        raise ValueError(f"{model_id}: training did not reach the common update cap")
    return recipe["validation_batch_size"]


def shared_neural_predictor(path, record, corpus, split, device, interface_path):
    batch_size = check_shared_neural_artifact(path, record, corpus.fixture, device)
    interface_record, interface = load_interface(interface_path, corpus, split, device)
    if record["metadata"]["interface_manifest_sha256"] != interface_record["artifact_manifest_sha256"]:
        raise ValueError(f"{record['model_id']}: shared expression interface digest mismatch")
    matrix_model = build_matrix_model(record["model_id"], corpus.gene_count, interface.channels,
                                      fixture=corpus.fixture)
    model = SharedReconstructor(matrix_model, interface)
    return load_neural(path, model, device, context=partial(attention_context, device),
                       deterministic=device == "cpu"), batch_size
