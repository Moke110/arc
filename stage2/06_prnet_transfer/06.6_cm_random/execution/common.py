"""Shared helpers for the ticket-#6 (06.6_cm_random) driver scripts.

Repo-root discovery per ADR 0030 §4 (walk up to the directory containing
``.git``; never count ``parents[N]``), plus the small hashing/JSON/CSV helpers the
queue wrapper and the record writer share.  The frozen-input verification and the
config fingerprint are **imported from the ticket-#4 pipeline**, never
re-implemented: ``ptlib.manifest`` (artifact normalisation/hashing) and
``ptlib.config`` (``config_fingerprint``) are the study's single definitions.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Mapping

__all__ = [
    "RepoRootError",
    "find_repo_root",
    "REPO_ROOT",
    "PIPELINE_EXEC",
    "SUBSTEP",
    "DEFAULT_DATA_DIR",
    "DEFAULT_OUTPUT",
    "pt_config",
    "pt_manifest",
    "sha256_file",
    "sha256_bytes",
    "canonical_json",
    "content_sha256",
    "load_json",
    "write_json_atomic",
    "write_text_atomic",
    "number",
    "frozen_identity",
    "data_identity_sha256",
    "read_locked_config",
]


class RepoRootError(RuntimeError):
    """No ``.git`` directory exists at or above the start path."""


def find_repo_root(start: Path | str | None = None) -> Path:
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RepoRootError(f"no .git directory at or above {path}")


REPO_ROOT = find_repo_root()
#: the study's only training/evaluation code (ticket #4)
PIPELINE_EXEC = REPO_ROOT / "stage2" / "06_prnet_transfer" / "06.4_pipeline" / "execution"
#: this substep
SUBSTEP = REPO_ROOT / "stage2" / "06_prnet_transfer" / "06.6_cm_random"
#: ticket #3's frozen input
DEFAULT_DATA_DIR = REPO_ROOT / "stage2" / "06_prnet_transfer" / "06.3_input" / "output"
#: this substep's gitignored compute output
DEFAULT_OUTPUT = SUBSTEP / "output"

if str(PIPELINE_EXEC) not in sys.path:
    sys.path.insert(0, str(PIPELINE_EXEC))

from ptlib import config as pt_config  # noqa: E402
from ptlib import manifest as pt_manifest  # noqa: E402


def sha256_file(path: Path | str, chunk_size: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json(payload: Any) -> bytes:
    """Deterministic JSON bytes used for every identity hash in this substep."""
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def content_sha256(payload: Any) -> str:
    return sha256_bytes(canonical_json(payload))


def load_json(path: Path | str) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_text_atomic(path: Path | str, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def write_json_atomic(path: Path | str, payload: Any) -> None:
    write_text_atomic(path, json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def number(value: Any) -> str:
    """Exact round-trip text for a float column cell (``repr`` round-trips JSON floats)."""
    if value is None:
        return ""
    value = float(value)
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return repr(value)


def read_locked_config(data_dir: Path | str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Read + normalise the frozen locked config (ticket #3); returns (config, notes)."""
    return pt_config.read_locked_config(Path(data_dir) / "locked_config.json")


def frozen_identity(data_dir: Path | str) -> dict[str, Any]:
    """Identity of the frozen input: manifest hash, locked-config hash, artifact hashes.

    Cheap by design: it re-reads the two small files and trusts the frozen
    manifest's recorded artifact sha256 values.  Every training run re-verified
    those artifacts (and records the result in its own manifest); this substep's
    ``--check`` cross-checks the run manifests against this identity.
    """
    data_dir = Path(data_dir)
    manifest_file = data_dir / "manifest.json"
    locked_file = data_dir / "locked_config.json"
    if not manifest_file.exists():
        raise FileNotFoundError(f"frozen input manifest is missing: {manifest_file}")
    if not locked_file.exists():
        raise FileNotFoundError(f"frozen locked config is missing: {locked_file}")
    artifacts = pt_manifest.normalize_artifact_manifest(load_json(manifest_file))
    return {
        "dir": str(data_dir),
        "manifest_sha256": sha256_file(manifest_file),
        "locked_config_sha256": sha256_file(locked_file),
        "artifacts": dict(sorted(artifacts.items())),
    }


def data_identity_sha256(frozen: Mapping[str, Any]) -> str:
    """One digest for "which frozen bytes were these candidates trained on"."""
    return content_sha256(
        {
            "manifest_sha256": frozen["manifest_sha256"],
            "locked_config_sha256": frozen["locked_config_sha256"],
            "artifacts": frozen["artifacts"],
        }
    )
