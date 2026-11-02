"""PRnet checkpoint payloads, strict warm start and the gene-axis guard.

Both stages save the **identical** payload::

    {"model": "prnet", "state_dict": ..., "config": ..., "seed": ...,
     "n_genes": ..., "best_valid_auprc": ..., "epoch": ...}

Warm start (``--init pretrained --init-checkpoint PATH``) is strict by design:

* key names must match exactly (no positional load, no partial load, no remap —
  a checkpoint outside this pipeline is rejected, never repaired silently);
* ``n_genes`` must equal the resolved gene axis;
* the checkpoint's structure must be shape-compatible (``structure_compatible``);
* the checkpoint sha256 and source path are recorded in the run manifest.

Any mismatch raises :class:`CheckpointError` before training starts.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

import torch

from . import config as config_mod
from . import model as model_mod

__all__ = [
    "CheckpointError",
    "PAYLOAD_KEYS",
    "CHECKPOINT_RELPATH",
    "save_checkpoint",
    "read_checkpoint",
    "checkpoint_sha256",
    "structure_from_checkpoint",
    "validate_checkpoint",
    "load_into_model",
    "build_from_checkpoint",
]

PAYLOAD_KEYS = (
    "model",
    "state_dict",
    "config",
    "seed",
    "n_genes",
    "best_valid_auprc",
    "epoch",
)

#: canonical run-relative location of the best checkpoint
CHECKPOINT_RELPATH = "checkpoint/model.pt"


class CheckpointError(RuntimeError):
    """A checkpoint is missing, malformed, or incompatible with the run."""


def checkpoint_sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(64 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_checkpoint(
    path: Path | str,
    model: torch.nn.Module,
    config: Mapping[str, Any],
    seed: int,
    n_genes: int,
    best_valid_auprc: float,
    epoch: int,
) -> None:
    """Write the standard payload (state dict on CPU)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    state = {key: value.detach().cpu() for key, value in model.state_dict().items()}
    payload = {
        "model": "prnet",
        "state_dict": state,
        "config": dict(config),
        "seed": int(seed),
        "n_genes": int(n_genes),
        "best_valid_auprc": float(best_valid_auprc),
        "epoch": int(epoch),
    }
    tmp = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, tmp)
    tmp.replace(path)


def read_checkpoint(path: Path | str, map_location: str | torch.device = "cpu") -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise CheckpointError(f"checkpoint not found: {path}")
    try:
        payload = torch.load(path, map_location=map_location, weights_only=False)
    except Exception as exc:  # noqa: BLE001 - torch raises many types
        raise CheckpointError(f"cannot read checkpoint {path}: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise CheckpointError(f"checkpoint {path} is not a dict payload")
    missing = [key for key in PAYLOAD_KEYS if key not in payload]
    if missing:
        raise CheckpointError(
            f"checkpoint {path} is not a pipeline checkpoint: missing keys {missing} "
            f"(present: {sorted(payload)})"
        )
    return dict(payload)


def structure_from_checkpoint(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Recover the canonical structure a checkpoint was trained with."""
    config = payload.get("config")
    if not isinstance(config, Mapping):
        raise CheckpointError("checkpoint payload has no config object")
    raw = config.get("structure", config)
    try:
        struct = config_mod.normalize_structure(raw)
    except config_mod.ConfigError as exc:
        raise CheckpointError(f"checkpoint config structure is unreadable: {exc}") from exc
    struct.setdefault("n_genes", payload.get("n_genes"))
    try:
        config_mod.validate_structure(struct)
    except config_mod.ConfigError as exc:
        raise CheckpointError(f"checkpoint structure is invalid: {exc}") from exc
    return struct


def validate_checkpoint(
    payload: Mapping[str, Any],
    expected_n_genes: int,
    expected_structure: Mapping[str, Any],
) -> dict[str, Any]:
    """Fail loudly on an axis/structure mismatch; return the checkpoint structure."""
    checkpoint_struct = structure_from_checkpoint(payload)
    n_genes = payload.get("n_genes")
    if n_genes is None or int(n_genes) != int(expected_n_genes):
        raise CheckpointError(
            f"gene-axis mismatch: checkpoint n_genes={n_genes!r}, run expects {expected_n_genes}"
        )
    if int(checkpoint_struct.get("n_genes", n_genes)) != int(expected_n_genes):
        raise CheckpointError(
            "gene-axis mismatch inside the checkpoint config: "
            f"{checkpoint_struct.get('n_genes')!r} != {expected_n_genes}"
        )
    ok, mismatches = config_mod.structure_compatible(checkpoint_struct, expected_structure)
    if not ok:
        raise CheckpointError(
            "checkpoint structure is incompatible with the run structure: "
            + "; ".join(mismatches)
        )
    if payload.get("model") != "prnet":
        raise CheckpointError(f"checkpoint model is {payload.get('model')!r}, expected 'prnet'")
    if not isinstance(payload.get("state_dict"), Mapping) or not payload["state_dict"]:
        raise CheckpointError("checkpoint state_dict is missing or empty")
    return checkpoint_struct


def load_into_model(model: torch.nn.Module, payload: Mapping[str, Any]) -> None:
    """Strict ``load_state_dict``: no missing key, no unexpected key, no re-init."""
    try:
        incompatible = model.load_state_dict(payload["state_dict"], strict=True)
    except RuntimeError as exc:
        raise CheckpointError(f"strict state_dict load failed: {exc}") from exc
    missing = getattr(incompatible, "missing_keys", [])
    unexpected = getattr(incompatible, "unexpected_keys", [])
    if missing or unexpected:
        raise CheckpointError(
            f"strict state_dict load left missing={missing} unexpected={unexpected}"
        )


def build_from_checkpoint(
    payload: Mapping[str, Any],
    expected_n_genes: int,
    expected_structure: Mapping[str, Any],
    device: str | torch.device = "cpu",
) -> tuple[torch.nn.Module, dict[str, Any]]:
    """Validate, build the model from the checkpoint structure and load strictly."""
    checkpoint_struct = validate_checkpoint(payload, expected_n_genes, expected_structure)
    model = model_mod.build_model(checkpoint_struct)
    load_into_model(model, payload)
    model.to(device)
    return model, checkpoint_struct
