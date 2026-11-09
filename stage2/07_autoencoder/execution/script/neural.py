"""Bounded joint training and validation checkpoint selection for Step 07."""

from contextlib import nullcontext
import hashlib
import os
from pathlib import Path
import sys
import time

import numpy as np
import torch
from torch.nn import functional as F

from benchmark import SOURCES, batches, json_bytes, publish_artifact, verify_artifact


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


def train_epoch(corpus, model, optimizer, train_rows, recipe, epoch, device, context, deadline):
    order = train_rows.copy()
    np.random.Generator(np.random.PCG64([recipe["seed"], epoch])).shuffle(order)
    model.train()
    total = torch.zeros((), device=device, dtype=torch.float64)
    batch_size = recipe["train_batch_size"]
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
        for epoch in range(1, recipe["epochs"] + 1):
            record = train_epoch(
                corpus, model, optimizer, train_rows, recipe, epoch, device, context, deadline,
            )
            model.eval()
            scores = validation_mse(
                corpus, split, predict, recipe["validation_batch_size"], deadline,
            )
            record["validation_mse"] = scores
            history.append(record)
            if scores["tenv"] < best_mse:
                best_mse = scores["tenv"]
                metadata["selected_epoch"] = epoch
                torch.save(
                    {name: tensor.detach().cpu() for name, tensor in model.state_dict().items()},
                    directory / "weights.pt",
                )
            print(json_bytes({"model_id": model_id, **record}).decode(), end="",
                  file=sys.stderr, flush=True)
        check_deadline(deadline)
        (directory / "epochs.json").write_bytes(json_bytes(history))
        metadata["selected_validation_mse"] = history[metadata["selected_epoch"] - 1]["validation_mse"]
        metadata["fit_elapsed_seconds"] = time.perf_counter() - start
        if device == "cuda":
            metadata["runtime"]["peak_cuda_allocated_bytes"] = torch.cuda.max_memory_allocated()
            metadata["runtime"]["peak_cuda_reserved_bytes"] = torch.cuda.max_memory_reserved()
    return verify_artifact(output, corpus, split)


def load_neural(path, model, device, *, context=nullcontext, deterministic=True):
    configure_runtime(device, deterministic)
    state = torch.load(Path(path) / "weights.pt", map_location="cpu", weights_only=True)
    model.load_state_dict(state, strict=True)
    return neural_predictor(model, device, context)
