"""Calibrate and freeze the shared continuous expression representation."""

from pathlib import Path

import numpy as np
import torch
from torch import nn

from benchmark import SEED, SOURCES, publish_artifact, verify_artifact


INTERFACE_ID = "shared_expression"
SCALAR_BATCH = 8192
CALIBRATION_STEPS = 1500
SCORE_BATCH = 16


class ExpressionInterface(nn.Module):
    def __init__(self, gene_count, channels):
        super().__init__()
        self.gene_count = gene_count
        self.channels = channels
        self.generator = nn.Sequential(nn.Linear(1, 64), nn.GELU(), nn.Linear(64, channels))
        self.restorer = nn.Sequential(nn.Linear(channels, 64), nn.GELU(), nn.Linear(64, 1))
        angles = torch.arange(gene_count, dtype=torch.float32)[:, None] * (
            10000 ** (-torch.arange(0, channels, 2, dtype=torch.float32) / channels))
        positions = torch.stack((angles.sin(), angles.cos()), dim=-1).flatten(1)
        self.register_buffer("positions", positions, persistent=False)

    def embed(self, values):
        return self.generator(values.unsqueeze(-1))

    def generate(self, values):
        return self.embed(values) + self.positions

    def readout(self, embeddings):
        return self.restorer(embeddings).squeeze(-1)

    def restore(self, representation):
        return self.readout(representation - self.positions)


def dimensions(corpus):
    return (min(3, corpus.gene_count), 16) if corpus.fixture else (1024, 256)


def check_device(corpus, device):
    if not corpus.fixture and device != "cuda":
        raise ValueError("real shared expression commands require --device cuda")
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA is unavailable for shared expression commands")


def fit_interface(corpus, split, output, device):
    check_device(corpus, device)
    if Path(output).exists():
        raise ValueError(f"artifact already exists: {output}")
    channels = dimensions(corpus)[1]
    torch.set_num_threads(2)
    torch.manual_seed(SEED)
    torch.set_float32_matmul_precision("highest")
    model = ExpressionInterface(corpus.gene_count, channels).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.002)
    rng = np.random.default_rng(SEED)
    # Sampling profiles with replacement from their joint membership and using
    # each selected gene with the same probability gives a uniform scalar draw.
    membership = [split.rows[source]["train"] for source in SOURCES]
    counts = np.array([len(rows) for rows in membership])
    steps, scalar_batch = (80, 32) if corpus.fixture else (CALIBRATION_STEPS, SCALAR_BATCH)
    rows_per_step = 4 if corpus.fixture else 32
    for _ in range(steps):
        sampled = rng.integers(counts.sum(), size=rows_per_step)
        source_indices = np.searchsorted(counts.cumsum(), sampled, side="right")
        values = []
        for index, source in enumerate(SOURCES):
            rows = sampled[source_indices == index] - counts[:index].sum()
            if len(rows):
                profiles = corpus.read(source, membership[index][rows])
                columns = rng.integers(corpus.gene_count, size=(len(rows), scalar_batch // rows_per_step))
                values.append(np.take_along_axis(profiles, columns, axis=1).reshape(-1))
        target = torch.from_numpy(np.concatenate(values)).to(device)
        reconstructed = model.readout(model.embed(target))
        loss = nn.functional.mse_loss(reconstructed, target)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
    if not all(torch.isfinite(p).all().item() for p in model.parameters()):
        raise ValueError("nonfinite shared expression weights")
    metadata = {
        "fit_role": "train",
        "fit_profiles": dict(zip(SOURCES, map(int, counts))),
        "seed": SEED,
        "steps": steps,
        "scalar_batch": scalar_batch,
        "profiles_per_step": rows_per_step,
        "sampling": "joint train profiles and genes sampled uniformly with replacement",
        "generator": [1, 64, channels],
        "restorer": [channels, 64, 1],
        "activation": "GELU",
        "positions": "fixed sinusoidal, in release gene order",
        "optimizer": {"name": "AdamW", "learning_rate": 0.002, "weight_decay": 0.01},
        "objective": "scalar round-trip MSE",
        "trainable_parameters": sum(p.numel() for p in model.parameters()),
        "fit_device": device,
        "torch_version": str(torch.__version__),
    }
    with publish_artifact(output, corpus, split, INTERFACE_ID, metadata) as directory:
        torch.save(model.cpu().state_dict(), directory / "weights.pt")
    return verify_artifact(output, corpus, split)


def load_interface(path, corpus, split, device):
    check_device(corpus, device)
    path = Path(path)
    record = verify_artifact(path, corpus, split)
    channels = dimensions(corpus)[1]
    if (record["model_id"] != INTERFACE_ID or set(record["files"]) != {"weights.pt"}
            or record["metadata"]["generator"] != [1, 64, channels]
            or record["metadata"]["restorer"] != [channels, 64, 1]
            or record["metadata"]["fit_role"] != "train"):
        raise ValueError("invalid shared expression artifact")
    torch.set_num_threads(2)
    model = ExpressionInterface(corpus.gene_count, channels)
    model.load_state_dict(torch.load(path / "weights.pt", weights_only=True, map_location="cpu"))
    return record, model.requires_grad_(False).eval().to(device)


def interface_predictor(model, device):
    torch.set_float32_matmul_precision("highest")

    @torch.inference_mode()
    def predict(values):
        tensor = torch.from_numpy(values).to(device)
        return model.restore(model.generate(tensor)).cpu().numpy()

    return predict, SCORE_BATCH
