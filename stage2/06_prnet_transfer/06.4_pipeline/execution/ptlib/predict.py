"""Target-free inference: predicted expression, then Δ = prediction - input baseline.

The prediction path deliberately has **no access to treated expression**: for the
CM fold it reads the frozen obs and the fold's fixed c* only; for Tahoe it reads
the obs and the matched ``(plate, cell_line)`` control means only (a control
input, which is a model input, never a scored target).  A prediction directory
therefore cannot leak a held-out response, and the evaluator is the only
component that reads truth.

Determinism: ``model.eval()``, zero noise, fixed batch order.  ``--repeats`` runs
the same inference twice and records whether the two output arrays are
byte-identical.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import torch

from . import PIPELINE_VERSION
from . import checkpoint as checkpoint_mod
from . import config as config_mod
from . import data as data_mod
from . import manifest as manifest_mod
from . import model as model_mod
from . import paths
from . import runtime

__all__ = ["PredictOptions", "run_predict", "build_predict_source", "predict_rows"]

log = logging.getLogger("ptlib.predict")


@dataclass
class PredictOptions:
    stage: str
    checkpoint: Path
    out_dir: Path
    data_dir: Path = paths.DEFAULT_FROZEN_INPUT
    role: str = "test"
    fold: int | None = None
    device: str = "cuda"
    batch_size: int = 512
    repeats: int = 2
    hash_cache: Path | None = None
    table_cache: Path | None = None
    command: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.out_dir = Path(self.out_dir)
        self.data_dir = Path(self.data_dir)
        self.checkpoint = Path(self.checkpoint)
        if self.stage not in ("pretrain", "finetune"):
            raise SystemExit(f"--stage must be pretrain or finetune, got {self.stage!r}")
        if self.stage == "finetune" and self.fold is None:
            raise SystemExit("--fold is required for --stage finetune")
        if self.role not in ("train", "valid", "test"):
            raise SystemExit(f"--role must be train/valid/test, got {self.role!r}")
        if int(self.batch_size) < 1:
            raise SystemExit("--batch-size must be positive")


def build_predict_source(frozen: data_mod.FrozenInput, options: PredictOptions):
    structures = frozen.structures()
    if options.stage == "finetune":
        folds = data_mod.load_folds(frozen)
        return data_mod.CMPredictSource(frozen, folds, int(options.fold), structures)
    return data_mod.TahoePredictSource(frozen, structures)


@torch.no_grad()
def predict_rows(
    model: torch.nn.Module,
    source: Any,
    rows: np.ndarray,
    n_dim: int,
    batch_size: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(predicted_expression, delta)`` for *rows* (fixed order, zero noise)."""
    n_genes = int(source.frozen.n_genes)
    pred = np.empty((rows.size, n_genes), dtype=np.float32)
    delta = np.empty((rows.size, n_genes), dtype=np.float32)
    zero_noise = torch.zeros(int(batch_size), int(n_dim), device=device)
    for start in range(0, rows.size, int(batch_size)):
        stop = min(start + int(batch_size), rows.size)
        batch_rows = np.asarray(rows[start:stop], dtype=np.int64)
        baseline = source.baseline(batch_rows)
        label = source.label(batch_rows)
        noise = zero_noise[: batch_rows.size]
        mean, _var = model(
            torch.from_numpy(np.ascontiguousarray(baseline)).to(device),
            torch.from_numpy(np.ascontiguousarray(label)).to(device),
            noise,
        )
        values = mean.detach().cpu().numpy().astype(np.float32)
        pred[start:stop] = values
        delta[start:stop] = values - baseline
    return pred, delta


def run_predict(options: PredictOptions) -> dict[str, Any]:
    started = time.perf_counter()
    started_utc = manifest_mod.now_utc()
    options.out_dir.mkdir(parents=True, exist_ok=True)
    frozen = data_mod.open_frozen_input(
        options.data_dir, hash_cache=options.hash_cache, table_cache=options.table_cache
    )

    payload = checkpoint_mod.read_checkpoint(options.checkpoint)
    struct = checkpoint_mod.structure_from_checkpoint(payload)
    checkpoint_mod.validate_checkpoint(payload, frozen.n_genes, struct)
    device, device_name = runtime.resolve_device(options.device)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    model = model_mod.build_model(struct)
    checkpoint_mod.load_into_model(model, payload)
    model.to(device)
    # ``check_architecture`` leaves the model in train mode; set eval() after it so
    # dropout/BatchNorm statistics cannot make inference nondeterministic.
    model_mod.check_architecture(model, struct)
    model.eval()

    source = build_predict_source(frozen, options)
    rows = np.asarray(source.rows_for(options.role), dtype=np.int64)
    if rows.size == 0:
        raise data_mod.DataError(f"no rows for stage={options.stage} role={options.role} fold={options.fold}")

    repeats = max(1, int(options.repeats))
    t0 = time.perf_counter()
    pred_expr, delta = predict_rows(
        model, source, rows, int(struct["n_dim"]), int(options.batch_size), device
    )
    first_pass_s = time.perf_counter() - t0
    identical = True
    repeat_s = 0.0
    for _ in range(repeats - 1):
        t0 = time.perf_counter()
        pred_again, delta_again = predict_rows(
            model, source, rows, int(struct["n_dim"]), int(options.batch_size), device
        )
        repeat_s += time.perf_counter() - t0
        if not (
            np.array_equal(pred_again, pred_expr)
            and np.array_equal(delta_again, delta)
            and pred_again.tobytes() == pred_expr.tobytes()
        ):
            identical = False
    if not identical:
        raise data_mod.DataError(
            "repeated inference is not byte-identical: the deterministic prediction "
            "policy is violated (this is a pipeline bug, not a warning)"
        )

    meta = source.meta(rows)
    np.save(options.out_dir / "pred_expr.npy", pred_expr)
    np.save(options.out_dir / "delta.npy", delta)
    meta.to_csv(options.out_dir / "meta.csv", index=False)
    predictions_sha = manifest_mod.sha256_bytes(pred_expr.tobytes())

    manifest = manifest_mod.build_manifest(
        "predict",
        stage=options.stage,
        status="completed",
        pipeline_version=PIPELINE_VERSION,
        command=list(options.command or []),
        git=manifest_mod.git_info(paths.ROOT),
        data=frozen.record(),
        fold=options.fold,
        role=options.role,
        config=payload.get("config", {}),
        config_sha256=config_mod.sha256_payload(payload.get("config", {})),
        structure=struct,
        checkpoint={
            "path": str(options.checkpoint),
            "sha256": checkpoint_mod.checkpoint_sha256(options.checkpoint),
            "epoch": int(payload.get("epoch", -1)),
            "best_valid_auprc": float(payload.get("best_valid_auprc", float("nan"))),
            "seed": int(payload.get("seed", -1)),
            "n_genes": int(payload.get("n_genes", -1)),
        },
        device=manifest_mod.device_record(options.device, str(device), device_name),
        noise_policy={
            "validation": "zero noise, model.eval(), fixed batch order",
            "inference": "zero noise, model.eval(), fixed batch order",
        },
        rows=int(rows.size),
        batch_size=int(options.batch_size),
        repeats=int(repeats),
        repeated_inference_identical=bool(identical),
        predictions={
            "pred_expr_sha256": predictions_sha,
            "pred_expr_shape": [int(x) for x in pred_expr.shape],
            "meta": "meta.csv",
        },
        timings={
            "started_utc": started_utc,
            "finished_utc": manifest_mod.now_utc(),
            "wall_s": time.perf_counter() - started,
            "first_pass_s": first_pass_s,
            "repeat_s": repeat_s,
        },
        peak_rss_mb=manifest_mod.peak_rss_mb(),
        warnings=[],
    )
    manifest_mod.write_json(options.out_dir / "manifest.json", manifest)
    log.info(
        "predicted %s/%s fold=%s: %d rows x %d genes in %.1fs "
        "(repeated inference byte-identical: %s)",
        options.stage, options.role, options.fold, rows.size, frozen.n_genes,
        first_pass_s, identical,
    )
    return {
        "out_dir": str(options.out_dir),
        "rows": int(rows.size),
        "repeated_inference_identical": identical,
        "first_pass_s": first_pass_s,
    }
