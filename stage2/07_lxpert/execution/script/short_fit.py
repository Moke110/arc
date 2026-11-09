"""Train and independently replay a bounded Step 07 autoencoder checkpoint."""
from contextlib import contextmanager
import csv
import json
from pathlib import Path
import resource
import tempfile
import time

import numpy as np
import torch
from torch.nn.attention import SDPBackend, sdpa_kernel

from autoencoder import TranscriptomeAutoencoder, architecture, reconstruction_mse
from inputs import SOURCES, load_release, release_files, sha256, verify


ATOL = 1e-6
RTOL = 1e-5
CHECKPOINT_SCHEMA = 1


class ExpressionRelease:
    def __init__(self, path, fixture, source_pins):
        manifest = load_release(path)
        if manifest["fixture"] is not fixture:
            raise ValueError("fixture release requires --fixture, and --fixture requires a fixture release")
        manifest_sha = sha256(path / "manifest.json")
        self.manifest = verify(path, None, source_pins)
        if sha256(path / "manifest.json") != manifest_sha:
            raise ValueError("input manifest changed during verification")
        self.contract = {
            "input_manifest_sha256": manifest_sha,
            "axis_sha256": self.manifest["axis_sha256"],
            "normalization": self.manifest["normalization"],
            "architecture": architecture(self.manifest["gene_count"]),
        }
        self.expressions, self.masks, self.profiles = {}, {}, {}
        for source in SOURCES:
            expression, mask, profiles = release_files(path, source)
            self.expressions[source] = np.load(expression, mmap_mode="r", allow_pickle=False)
            self.masks[source] = np.load(mask, allow_pickle=False)
            with profiles.open(newline="") as file:
                self.profiles[source] = list(csv.DictReader(file))

    def rows(self, source, role):
        return np.array([row for row, profile in enumerate(self.profiles[source])
                         if profile["development"] == role], dtype=np.int64)

    def batch(self, rows, device):
        values = np.stack([self.expressions[source][row] for source, row in rows])
        masks = np.stack([self.masks[source] for source, _ in rows])
        return torch.from_numpy(values).to(device), torch.from_numpy(masks).to(device)


def runtime(device, fixture):
    if device == "cpu":
        if not fixture:
            raise ValueError("full-axis fitting and replay require CUDA")
    elif not torch.cuda.is_available() or not torch.cuda.is_bf16_supported():
        raise ValueError("CUDA with native BF16 support is required")
    else:
        torch.cuda.reset_peak_memory_stats()
    return {"device": device, "precision": "bfloat16" if device == "cuda" else "float32"}


@contextmanager
def attention_runtime(device):
    if device == "cuda":
        # An unsupported CUDA configuration fails instead of materializing a gene-by-token matrix.
        with sdpa_kernel(SDPBackend.FLASH_ATTENTION), torch.autocast("cuda", dtype=torch.bfloat16):
            yield
    else:
        yield


def measure(start, device):
    return {
        "elapsed_seconds": time.perf_counter() - start,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "peak_cuda_allocated_bytes": torch.cuda.max_memory_allocated() if device == "cuda" else 0,
        "peak_cuda_reserved_bytes": torch.cuda.max_memory_reserved() if device == "cuda" else 0,
    }


def probe(release, model, rows, device):
    model.eval()
    expression, observed = release.batch(rows, device)
    with torch.inference_mode(), attention_runtime(device):
        latent = model.encode(expression)
        prediction = model.decode(latent).float()
        losses = reconstruction_mse(prediction, expression, observed)
    if not torch.isfinite(latent).all() or not torch.isfinite(prediction).all():
        raise ValueError("checkpoint probe produced non-finite values")
    return prediction.cpu(), {
        "probe_rows": rows,
        "latent_shape": list(latent.shape),
        "reconstruction_shape": list(prediction.shape),
        "finite": True,
        "probe_mse": {source: float(loss) for (source, _), loss in zip(rows, losses, strict=True)},
    }


def smoke_batches(release, steps, batch_size, seed):
    if steps < 1 or batch_size < 2:
        raise ValueError("short fit requires at least one step and batch size of at least two")
    rng = np.random.default_rng(seed)
    by_source = {}
    for index, source in enumerate(SOURCES):
        count = batch_size // 2 + (index == 0 and batch_size % 2)
        rows = release.rows(source, "train")
        if count * steps > len(rows):
            raise ValueError(f"{source}: insufficient development training rows for this short fit")
        by_source[source] = rng.permutation(rows)[:count * steps].reshape(steps, count)
    return [[(source, int(row)) for source in SOURCES for row in by_source[source][step]]
            for step in range(steps)]


def fit(args, source_pins):
    start = time.perf_counter()
    if args.output.exists():
        raise ValueError(f"fit output already exists: {args.output}")
    execution = runtime(args.device, args.fixture)
    release = ExpressionRelease(args.release, args.fixture, source_pins)
    batches = smoke_batches(release, args.steps, args.batch_size, args.seed)
    torch.manual_seed(args.seed)
    model = TranscriptomeAutoencoder(release.manifest["gene_count"]).to(args.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    step_seconds, losses = [], []
    for rows in batches:
        step_start = time.perf_counter()
        expression, observed = release.batch(rows, args.device)
        optimizer.zero_grad(set_to_none=True)
        with attention_runtime(args.device):
            loss = reconstruction_mse(model(expression), expression, observed).mean()
        loss_value = float(loss.detach())
        if not np.isfinite(loss_value):
            raise ValueError("non-finite reconstruction training loss")
        loss.backward()
        optimizer.step()
        if args.device == "cuda":
            torch.cuda.synchronize()
        step_seconds.append(time.perf_counter() - step_start)
        losses.append(loss_value)
    probe_rows = [(source, int(release.rows(source, "validation")[0])) for source in SOURCES]
    prediction, probe_report = probe(release, model, probe_rows, args.device)
    metadata = {"schema_version": CHECKPOINT_SCHEMA, **release.contract, **execution}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.output.parent, prefix=".fit-") as staging:
        directory = Path(staging)
        checkpoint = directory / "checkpoint.pt"
        torch.save({"metadata": metadata, "model": model.state_dict(), "probe_rows": probe_rows,
                    "probe_prediction": prediction}, checkpoint)
        report = {
            **metadata, **probe_report, **measure(start, args.device),
            "checkpoint_sha256": sha256(checkpoint),
            "seed": args.seed, "steps": args.steps, "batch_size": args.batch_size,
            "step_seconds": step_seconds, "step_loss": losses,
            "training_profiles": {source: sum(row_source == source for rows in batches
                                                for row_source, _ in rows) for source in SOURCES},
            "observed_genes": {source: int(release.masks[source].sum()) for source in SOURCES},
        }
        (directory / "fit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        directory.rename(args.output)
    return report


def load_checkpoint(checkpoint, contract, device):
    payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
    metadata = payload["metadata"]
    for key, value in {"schema_version": CHECKPOINT_SCHEMA, **contract}.items():
        if metadata[key] != value:
            raise ValueError(f"checkpoint metadata mismatch: {key}")
    model = TranscriptomeAutoencoder(contract["architecture"]["gene_count"])
    model.load_state_dict(payload["model"], strict=True)
    return model.to(device), payload


def reload(args, source_pins):
    start = time.perf_counter()
    execution = runtime(args.device, args.fixture)
    release = ExpressionRelease(args.release, args.fixture, source_pins)
    model, payload = load_checkpoint(args.checkpoint, release.contract, args.device)
    if any(payload["metadata"][key] != value for key, value in execution.items()):
        raise ValueError("checkpoint probe replay requires the fit device and precision")
    prediction, report = probe(release, model, payload["probe_rows"], args.device)
    reference = payload["probe_prediction"]
    if prediction.shape != reference.shape or not torch.allclose(prediction, reference, atol=ATOL, rtol=RTOL):
        raise ValueError("checkpoint probe reconstruction exceeded replay tolerance")
    return {**release.contract, **execution, **report, **measure(start, args.device),
            "checkpoint_sha256": sha256(args.checkpoint), "atol": ATOL, "rtol": RTOL,
            "max_abs_error": float((prediction - reference).abs().max())}
