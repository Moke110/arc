"""Bounded joint training and validation checkpoint selection for Step 07."""

from contextlib import nullcontext
import hashlib
import json
from math import ceil
import os
from pathlib import Path
import sys
import time

import numpy as np
import torch
from torch.nn import functional as F

from benchmark import (
    SOURCES, batches, check_checkpoint_selection, json_bytes, publish_artifact, verify_artifact,
)


def check_deadline(deadline):
    if time.perf_counter() >= deadline:
        raise ValueError("fit exceeded its frozen wall time ceiling")


def neural_predictor(model, device, context=nullcontext):
    model.to(device).eval()

    @torch.inference_mode()
    def predict(values):
        expression = torch.from_numpy(values).to(device)
        with context():
            prediction = model(expression)
        return prediction.float().cpu().numpy()

    return predict


def validation_mse(corpus, split, predict, batch_size, deadline):
    scores = {}
    for source in SOURCES:
        sse = 0.0
        for values in batches(corpus, split, source, "validation", batch_size):
            check_deadline(deadline)
            prediction = predict(values)
            if prediction.shape != values.shape or not np.isfinite(prediction).all():
                raise ValueError(f"{source}: invalid validation prediction")
            residual = np.subtract(prediction, values, dtype=np.float64)
            sse += float(np.square(residual, out=residual).sum())
        scores[source] = sse / (len(split.rows[source]["validation"]) * corpus.gene_count)
        if not np.isfinite(scores[source]):
            raise ValueError(f"{source}: nonfinite validation MSE")
    return scores


def train_epoch(corpus, model, optimizer, train_rows, recipe, epoch, device, context, deadline,
                batch_limit=None):
    order = train_rows.copy()
    np.random.Generator(np.random.PCG64([recipe["seed"], epoch])).shuffle(order)
    batch_size = recipe["train_batch_size"]
    if batch_limit is not None:
        order = order[:batch_limit * batch_size]
    model.train()
    total = torch.zeros((), device=device, dtype=torch.float64)
    for start in range(0, len(order), batch_size):
        check_deadline(deadline)
        members = order[start:start + batch_size]
        values = np.empty((len(members), corpus.gene_count), dtype=np.float32)
        for index, source in enumerate(SOURCES):
            source_mask = members[:, 0] == index
            if source_mask.any():
                values[source_mask] = corpus.read(source, members[source_mask, 1])
        expression = torch.from_numpy(values).to(device)
        optimizer.zero_grad(set_to_none=True)
        with context():
            prediction = model(expression)
        loss = F.mse_loss(prediction.float(), expression)
        loss.backward()
        optimizer.step()
        total += loss.detach().double() * len(members)
    mse = float(total) / len(order)
    if not np.isfinite(mse):
        raise ValueError("nonfinite training MSE")
    return {
        "epoch": epoch,
        "training_mse": mse,
        "sample_order_sha256": hashlib.sha256(order.astype("<i8", copy=False).tobytes()).hexdigest(),
    }


def configure_runtime(device, deterministic):
    if deterministic:
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA device is unavailable")
    torch.set_float32_matmul_precision("highest")
    torch.use_deterministic_algorithms(deterministic)


def fit_neural(corpus, split, output, model_id, build_model, recipe, device, metadata,
               *, context=nullcontext, deterministic=True):
    """Fit a module accepting complete profiles and seal its best inference state."""
    configure_runtime(device, deterministic)
    start = time.perf_counter()
    deadline = start + recipe["max_seconds"]
    torch.manual_seed(recipe["seed"])
    train_rows = np.concatenate([
        np.column_stack((np.full(len(split.rows[source]["train"]), index),
                         split.rows[source]["train"]))
        for index, source in enumerate(SOURCES)
    ])
    metadata = {
        **metadata,
        "fit_role": "train",
        "fit_profiles": {source: len(split.rows[source]["train"]) for source in SOURCES},
        "selection_role": "validation",
        "selection_metric": "tenv MSE, earliest epoch on ties",
        "validation_profiles": {source: len(split.rows[source]["validation"]) for source in SOURCES},
        "recipe": recipe,
        "runtime": {
            "device": device,
            "torch_version": str(torch.__version__),
            "cuda_version": torch.version.cuda if device == "cuda" else None,
            "gpu": torch.cuda.get_device_name() if device == "cuda" else None,
            "float32_matmul_precision": torch.get_float32_matmul_precision(),
            "deterministic_algorithms": deterministic,
            "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG")
            if device == "cuda" else None,
        },
    }
    with publish_artifact(output, corpus, split, model_id, metadata) as directory:
        if device == "cuda":
            torch.cuda.reset_peak_memory_stats()
        model = build_model().to(device)
        optimizer = torch.optim.AdamW(
            model.parameters(), lr=recipe["learning_rate"], weight_decay=recipe["weight_decay"],
        )
        predict = neural_predictor(model, device, context)
        best_mse = float("inf")
        history = []
        max_steps = recipe.get("max_optimizer_steps")
        steps_per_epoch = ceil(len(train_rows) / recipe["train_batch_size"])
        epochs = recipe["epochs"] if max_steps is None else ceil(max_steps / steps_per_epoch)
        steps = 0
        stale_checks = 0
        for epoch in range(1, epochs + 1):
            batch_limit = None if max_steps is None else min(steps_per_epoch, max_steps - steps)
            record = train_epoch(
                corpus, model, optimizer, train_rows, recipe, epoch, device, context, deadline,
                batch_limit=batch_limit,
            )
            if batch_limit is not None:
                steps += batch_limit
                record["optimizer_steps"] = steps
                record["training_profiles"] = min(len(train_rows),
                                                  batch_limit * recipe["train_batch_size"])
            model.eval()
            scores = validation_mse(
                corpus, split, predict, recipe["validation_batch_size"], deadline,
            )
            record["validation_mse"] = scores
            history.append(record)
            if scores["tenv"] < best_mse:
                best_mse = scores["tenv"]
                stale_checks = 0
                metadata["selected_epoch"] = epoch
                torch.save(
                    {name: tensor.detach().cpu() for name, tensor in model.state_dict().items()},
                    directory / "weights.pt",
                )
            else:
                stale_checks += 1
            print(json_bytes({"model_id": model_id, **record}).decode(), end="",
                  file=sys.stderr, flush=True)
            if max_steps is not None:
                if steps == max_steps:
                    metadata["stop_reason"] = "optimizer step limit"
                    break
                if stale_checks >= recipe["validation_patience"]:
                    metadata["stop_reason"] = "validation patience"
                    break
        check_deadline(deadline)
        (directory / "epochs.json").write_bytes(json_bytes(history))
        metadata["selected_validation_mse"] = history[metadata["selected_epoch"] - 1]["validation_mse"]
        metadata["fit_elapsed_seconds"] = time.perf_counter() - start
        if device == "cuda":
            metadata["runtime"]["peak_cuda_allocated_bytes"] = torch.cuda.max_memory_allocated()
            metadata["runtime"]["peak_cuda_reserved_bytes"] = torch.cuda.max_memory_reserved()
    return verify_artifact(output, corpus, split)


def check_selected_checkpoint(path, metadata, model_id, device):
    if (metadata["runtime"]["device"] != device
            or metadata["selection_role"] != "validation"
            or metadata["selection_metric"] != "tenv MSE, earliest epoch on ties"):
        raise ValueError(f"{model_id}: incompatible checkpoint selection or fit device")
    history = json.loads((Path(path) / "epochs.json").read_text())
    check_checkpoint_selection(history, metadata, model_id)


def load_neural(path, model, device, *, context=nullcontext, deterministic=True):
    configure_runtime(device, deterministic)
    state = torch.load(Path(path) / "weights.pt", map_location="cpu", weights_only=True)
    model.load_state_dict(state, strict=True)
    return neural_predictor(model, device, context)
