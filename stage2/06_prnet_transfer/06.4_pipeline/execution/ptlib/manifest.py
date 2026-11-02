"""Run manifests, artifact hashing and the reproducibility record.

Every training/prediction/evaluation invocation writes one manifest.  It is the
machine-readable record the brief requires: device, wall time, peak RSS, stop
reason, git SHA, data hashes, checkpoint identity (including the warm-start
source and its sha256), the noise policy, and the resolved config.

The frozen input's ``manifest.json`` is verified **before** any artifact is
read: a missing entry or a sha256 mismatch fails closed.  A per-run hash cache
(``--hash-cache``) can memoise the 16 GB Tahoe hash by ``(size, mtime_ns)``; its
use is recorded in the manifest so the record never claims a fresh read it did
not perform.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import resource
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

__all__ = [
    "ManifestError",
    "REQUIRED_ARTIFACTS",
    "SCHEMA",
    "now_utc",
    "sha256_file",
    "sha256_bytes",
    "normalize_artifact_manifest",
    "verify_artifacts",
    "git_info",
    "peak_rss_mb",
    "read_json",
    "write_json",
    "build_manifest",
    "validate_manifest",
    "manifest_path",
    "unit_complete",
    "hash_cache_path",
]

SCHEMA = "ptlib.run_manifest/v1"

#: Files ticket #3's frozen input contract guarantees.
REQUIRED_ARTIFACTS = (
    "cm_transfer.h5ad",
    "folds.parquet",
    "consensus_controls.npz",
    "consensus_membership.parquet",
    "tahoe_transfer.h5ad",
    "tahoe_control_groups.parquet",
    "tahoe_pairs.npz",
    "locked_config.json",
)

_KIND_REQUIRED = {
    "train": (
        "kind", "stage", "status", "created_utc", "pipeline_version", "command",
        "git", "data", "locked_config", "config", "config_sha256", "structure",
        "seed", "device", "init", "checkpoint", "noise_policy", "timings",
        "peak_rss_mb", "epochs_run", "stop_reason",
    ),
    "predict": (
        "kind", "stage", "status", "created_utc", "pipeline_version", "command",
        "git", "data", "config", "config_sha256", "checkpoint", "device",
        "noise_policy", "timings", "peak_rss_mb", "rows",
    ),
    "evaluate": (
        "kind", "stage", "status", "created_utc", "pipeline_version", "command",
        "git", "data", "timings", "peak_rss_mb", "predictions", "metric",
    ),
}


class ManifestError(RuntimeError):
    """A manifest or an artifact hash does not match its record."""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path | str, chunk_size: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def normalize_artifact_manifest(raw: Mapping[str, Any]) -> dict[str, str]:
    """Extract ``{relative name: sha256}`` from the frozen input's manifest.

    Ticket #3 may encode the manifest as a flat mapping, under ``artifacts`` or
    ``files`` as a mapping or a list of records.  All four shapes are accepted;
    anything else fails closed.  Values may be a bare hash or a record with
    ``sha256``/``content_sha256``.
    """
    if not isinstance(raw, Mapping):
        raise ManifestError("frozen manifest must be a JSON object")

    def _entries(block: Any) -> dict[str, str]:
        out: dict[str, str] = {}
        if isinstance(block, Mapping):
            items = block.items()
        elif isinstance(block, list):
            items = []
            for record in block:
                if not isinstance(record, Mapping):
                    raise ManifestError("artifact list entries must be objects")
                name = record.get("name") or record.get("path") or record.get("file")
                if name is None:
                    raise ManifestError("artifact list entry has no name/path/file")
                items.append((name, record))
        else:
            raise ManifestError(f"unsupported artifact manifest block: {type(block).__name__}")
        for name, record in items:
            if isinstance(record, str):
                digest = record
            elif isinstance(record, Mapping):
                digest = (
                    record.get("sha256")
                    or record.get("content_sha256")
                    or (record.get("hashes") or {}).get("sha256")
                )
                if digest is None and "md5" in record:
                    continue
            else:
                raise ManifestError(f"artifact {name!r} has no usable hash record")
            if not isinstance(digest, str) or len(digest) != 64:
                raise ManifestError(f"artifact {name!r} has no sha256 in the manifest")
            out[str(os.path.basename(str(name)))] = digest.lower()
        return out

    for key in ("artifacts", "files", "sha256", "hashes"):
        if key in raw and isinstance(raw[key], (Mapping, list)):
            block = raw[key]
            if key == "sha256" and isinstance(block, Mapping):
                entries = _entries(block)
                if entries:
                    return entries
            entries = _entries(block)
            if entries:
                return entries
    entries = _entries(raw)
    if not entries:
        raise ManifestError("frozen manifest contains no artifact hashes")
    return entries


def hash_cache_path(path: Path | str | None) -> Path | None:
    if path is None:
        return None
    return Path(path)


def _load_hash_cache(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _save_hash_cache(path: Path | None, payload: Mapping[str, Any]) -> None:
    """Persist the frozen-artifact hash memo; safe for concurrent queue workers.

    The memo path is shared by every unit of a parallel queue, so a fixed
    ``<name>.tmp`` makes two simultaneous completions race: both write the same
    temp file, the first ``replace`` moves it away, and the second raises
    ``FileNotFoundError`` — killing that training unit with a data-path error that
    has nothing to do with the data.  A unique temp file in the same directory
    keeps the write atomic without a cross-process lock; the last writer wins and
    every writer succeeds.
    """
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(tmp_name, path)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


def verify_artifacts(
    data_dir: Path | str,
    artifact_hashes: Mapping[str, str],
    required: Iterable[str] = REQUIRED_ARTIFACTS,
    hash_cache: Path | str | None = None,
) -> dict[str, Any]:
    """Verify every required artifact against its recorded sha256; fail closed.

    Returns a record ``{"verified": {name: sha}, "cached": [...],
    "bytes": total, "manifest_complete": bool}``.  Raises :class:`ManifestError`
    listing every missing artifact and every mismatch at once.
    """
    data_dir = Path(data_dir)
    cache_path = hash_cache_path(hash_cache)
    cache = _load_hash_cache(cache_path)
    problems: list[str] = []
    verified: dict[str, str] = {}
    cached: list[str] = []
    total_bytes = 0

    for name in required:
        if name not in artifact_hashes:
            problems.append(f"{name}: not listed in the frozen manifest")
            continue
        path = data_dir / name
        if not path.exists():
            problems.append(f"{name}: file is missing at {path}")
            continue
        stat = path.stat()
        total_bytes += stat.st_size
        expected = artifact_hashes[name].lower()
        entry = cache.get(str(path.resolve()))
        digest: str | None = None
        if (
            isinstance(entry, Mapping)
            and entry.get("size") == stat.st_size
            and entry.get("mtime_ns") == stat.st_mtime_ns
            and entry.get("sha256") == expected
        ):
            digest = expected
            cached.append(name)
        if digest is None:
            digest = sha256_file(path)
            if digest != expected:
                problems.append(
                    f"{name}: sha256 mismatch (manifest {expected[:12]}…, file {digest[:12]}…)"
                )
                continue
            cache[str(path.resolve())] = {
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": digest,
            }
        verified[name] = digest

    if problems:
        raise ManifestError(
            f"frozen input verification failed for {data_dir}:\n  - " + "\n  - ".join(problems)
        )
    _save_hash_cache(cache_path, cache)
    return {
        "verified": verified,
        "cached": cached,
        "bytes": int(total_bytes),
        "manifest_complete": True,
    }


def git_info(root: Path | str | None = None) -> dict[str, Any]:
    root = Path(root) if root is not None else Path.cwd()
    try:
        sha = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        branch = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True, text=True, check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return {"sha": None, "branch": None, "dirty": None}
    return {"sha": sha, "branch": branch, "dirty": bool(status.strip())}


def peak_rss_mb() -> float:
    """Peak resident set size of this process, MiB (Linux: ru_maxrss is KiB)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def read_json(path: Path | str) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _json_default(value: Any) -> Any:
    """Make numpy scalars/arrays, paths and sets serialisable."""
    import numpy as _np

    if isinstance(value, _np.integer):
        return int(value)
    if isinstance(value, _np.floating):
        return float(value)
    if isinstance(value, _np.bool_):
        return bool(value)
    if isinstance(value, _np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (set, tuple)):
        return list(value)
    raise TypeError(f"object of type {type(value).__name__} is not JSON serialisable")


def write_json(path: Path | str, payload: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=False, default=_json_default) + "\n",
        encoding="utf-8",
    )
    tmp.replace(path)


def device_record(requested: str, resolved: str, name: str) -> dict[str, Any]:
    import torch

    return {
        "requested": requested,
        "resolved": resolved,
        "name": name,
        "torch_version": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_version": torch.version.cuda,
        "python": platform.python_version(),
        "cudnn_deterministic": bool(torch.backends.cudnn.deterministic),
        "cudnn_benchmark": bool(torch.backends.cudnn.benchmark),
    }


def build_manifest(kind: str, **payload: Any) -> dict[str, Any]:
    """Assemble a run manifest; ``kind`` selects the required schema block."""
    manifest = {
        "schema": SCHEMA,
        "kind": kind,
        "created_utc": now_utc(),
        **payload,
    }
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    """Validate the manifest schema; raise :class:`ManifestError` on any gap."""
    if not isinstance(manifest, Mapping):
        raise ManifestError("run manifest must be a JSON object")
    if manifest.get("schema") != SCHEMA:
        raise ManifestError(f"run manifest schema must be {SCHEMA!r}, got {manifest.get('schema')!r}")
    kind = manifest.get("kind")
    if kind not in _KIND_REQUIRED:
        raise ManifestError(f"unknown manifest kind {kind!r}")
    missing = [key for key in _KIND_REQUIRED[kind] if key not in manifest]
    if missing:
        raise ManifestError(f"{kind} manifest is missing required field(s): {missing}")
    if manifest.get("status") not in ("completed", "failed"):
        raise ManifestError("manifest.status must be 'completed' or 'failed'")
    git = manifest.get("git")
    if not isinstance(git, Mapping) or "sha" not in git:
        raise ManifestError("manifest.git must be an object with a 'sha'")
    data = manifest.get("data")
    if not isinstance(data, Mapping) or "dir" not in data or "artifacts" not in data:
        raise ManifestError("manifest.data must hold 'dir' and 'artifacts'")
    timings = manifest.get("timings")
    if not isinstance(timings, Mapping) or "wall_s" not in timings:
        raise ManifestError("manifest.timings must hold 'wall_s'")
    if kind in ("train", "predict", "evaluate"):
        noise = manifest.get("noise_policy")
        if not isinstance(noise, Mapping):
            raise ManifestError("manifest.noise_policy must be an object")
    if kind == "train":
        if not isinstance(manifest.get("config"), Mapping):
            raise ManifestError("train manifest.config must be an object")
        if not isinstance(manifest.get("checkpoint"), Mapping):
            raise ManifestError("train manifest.checkpoint must be an object")
        init = manifest.get("init")
        if not isinstance(init, Mapping) or init.get("mode") not in ("random", "pretrained"):
            raise ManifestError("train manifest.init.mode must be 'random' or 'pretrained'")
        for key in ("epochs_run", "stop_reason", "seed", "structure"):
            if manifest.get(key) is None:
                raise ManifestError(f"train manifest.{key} is required")
    if kind == "predict":
        rows = manifest.get("rows")
        if not isinstance(rows, (int, list)):
            raise ManifestError("predict manifest.rows must be an int or a list of row ids")


def manifest_path(run_dir: Path | str) -> Path:
    return Path(run_dir) / "manifest.json"


def unit_complete(
    run_dir: Path | str,
    expected_config_sha256: str | None = None,
    require_checkpoint: bool = True,
) -> tuple[bool, str]:
    """Resume check for one queue unit: is this run directory a finished unit?"""
    run_dir = Path(run_dir)
    manifest_file = manifest_path(run_dir)
    if not manifest_file.exists():
        return False, "no manifest"
    try:
        manifest = read_json(manifest_file)
        validate_manifest(manifest)
    except (ManifestError, json.JSONDecodeError, OSError) as exc:
        return False, f"unreadable manifest: {exc}"
    if manifest.get("status") != "completed":
        return False, f"status={manifest.get('status')}"
    if expected_config_sha256 and manifest.get("config_sha256") != expected_config_sha256:
        return False, "config fingerprint differs"
    if require_checkpoint:
        checkpoint = manifest.get("checkpoint") or {}
        path = checkpoint.get("path")
        if not path or not (run_dir / path).exists():
            return False, "checkpoint missing"
        recorded = checkpoint.get("sha256")
        if recorded:
            actual = sha256_file(run_dir / path)
            if actual != recorded:
                return False, "checkpoint sha256 mismatch"
    return True, "complete"


def is_finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False
