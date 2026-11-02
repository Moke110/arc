"""Shared helpers for `06.9_evaluation` — roots, identity, hashes, atomic records.

This substep owns **no** metric, prediction or training code. Everything scientific comes
from `06.4_pipeline/execution` (`ptlib`), invoked through its CLIs or, for the model-card
writer, imported directly. This module only adds the record-level identity helpers ticket
#9 needs:

* the one **data identity** recipe used by #6/#7 (so the 30 units can be cross-checked
  against the checkpoints' own records, byte for byte);
* checkpoint-root search — the selected checkpoints are gitignored local weights, so a
  worktree record resolves against the worktree root and the main checkout;
* deterministic JSON/text writing for the tracked record.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "EvaluationError",
    "ROOT",
    "TRANSFER",
    "EVAL_DIR",
    "EXECUTION",
    "PIPELINE_EXECUTION",
    "DEFAULT_DATA_DIR",
    "DEFAULT_RANDOM_RECORD",
    "DEFAULT_TAHOE_RECORD",
    "DEFAULT_BASELINES_DIR",
    "DEFAULT_OUT",
    "DEFAULT_TABLE_CACHE",
    "DEFAULT_GPU_PY",
    "find_repo_root",
    "repo_roots",
    "resolve_checkpoint_path",
    "pipeline_modules",
    "sha256_file",
    "sha256_bytes",
    "canonical_json",
    "content_sha256",
    "load_json",
    "write_json",
    "write_text_atomic",
    "now_utc",
    "peak_rss_mb",
    "frozen_identity",
    "data_identity_sha256",
    "mem_available_gib",
    "format_lr_tag",
    "read_table",
]


class EvaluationError(RuntimeError):
    """A unit manifest, record or evaluation artifact is invalid (fail closed)."""


def find_repo_root(start: Path | str | None = None) -> Path:
    """Return the closest ancestor of *start* that contains ``.git`` (file or dir)."""
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise EvaluationError(f"no .git directory at or above {path}")


ROOT = find_repo_root()
#: ``<repo>/stage2/06_prnet_transfer/``
TRANSFER = ROOT / "stage2" / "06_prnet_transfer"
#: ``<repo>/stage2/06_prnet_transfer/06.9_evaluation/``
EVAL_DIR = TRANSFER / "06.9_evaluation"
EXECUTION = EVAL_DIR / "execution"
#: the only place the study's pipeline library lives (never re-implemented here)
PIPELINE_EXECUTION = TRANSFER / "06.4_pipeline" / "execution"

DEFAULT_DATA_DIR = TRANSFER / "06.3_input" / "output"
DEFAULT_RANDOM_RECORD = (
    TRANSFER / "06.6_cm_random" / "output" / "record" / "selected_checkpoints.csv"
)
DEFAULT_TAHOE_RECORD = (
    TRANSFER / "06.7_cm_tahoe" / "output" / "record" / "selected_checkpoints.csv"
)
DEFAULT_BASELINES_DIR = TRANSFER / "06.8_baselines" / "output"
DEFAULT_OUT = EVAL_DIR / "output"
#: derived CSV cache of the frozen parquet tables (built by 06.4; pinned env has no engine)
DEFAULT_TABLE_CACHE = TRANSFER / "06.4_pipeline" / "output" / "input_cache"
DEFAULT_GPU_PY = "/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python"

#: abort a heavy run if MemAvailable (GiB) is below this — CM evaluation is far smaller
#: than the 8.5 GiB full Tahoe-valid evaluation, but leave the repo-wide headroom rule.
MIN_AVAILABLE_GIB = 6.0


# --------------------------------------------------------------------------- #
# Checkpoint roots (selected checkpoints are gitignored local weights)
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def repo_roots() -> tuple[Path, ...]:
    """Places where a worktree-relative record path may resolve.

    Order: ``$ARC_CHECKPOINT_ROOTS`` (colon-separated), this worktree, then every
    worktree listed by git (the main checkout comes first in ``git worktree list``).
    """
    roots: list[Path] = []
    override = os.environ.get("ARC_CHECKPOINT_ROOTS", "").strip()
    if override:
        roots.extend(Path(item) for item in override.split(os.pathsep) if item)
    roots.append(ROOT)
    try:
        done = subprocess.run(
            ["git", "-C", str(ROOT), "worktree", "list", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        if done.returncode == 0:
            for line in done.stdout.splitlines():
                if line.startswith("worktree "):
                    roots.append(Path(line[len("worktree ") :].strip()))
    except (OSError, subprocess.SubprocessError):  # pragma: no cover - odd git setups
        pass
    unique: list[Path] = []
    for root in roots:
        resolved = root.resolve()
        if resolved not in unique:
            unique.append(resolved)
    return tuple(unique)


def resolve_checkpoint_path(relpath: str | Path) -> tuple[Path, Path]:
    """Resolve a record's checkpoint path against the candidate roots.

    Returns ``(path, root)``; raises :class:`EvaluationError` listing every root tried.
    """
    candidate = Path(relpath)
    if candidate.is_absolute() and candidate.exists():
        return candidate, Path(candidate.anchor)
    tried = [candidate]
    for root in repo_roots():
        path = root / candidate
        tried.append(path)
        if path.exists():
            return path, root
    raise EvaluationError(
        "selected checkpoint is missing; tried:\n  - "
        + "\n  - ".join(str(item) for item in tried)
        + "\n(the checkpoints are gitignored local weights: run from a tree that holds "
        "06.6_cm_random/output/runs and 06.7_cm_tahoe/output/runs, or set "
        "ARC_CHECKPOINT_ROOTS)"
    )


# --------------------------------------------------------------------------- #
# The pipeline library (single source of truth for digests / manifests)
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def pipeline_modules() -> dict[str, Any]:
    """Import the ticket-#4 helpers this substep is allowed to share.

    Only stdlib-only modules are imported here (`ptlib.config`,
    `ptlib.manifest`); the evaluator/predictor/model-card code is invoked through its
    CLIs or imported lazily where a run directory is actually read.
    """
    if str(PIPELINE_EXECUTION) not in sys.path:
        sys.path.insert(0, str(PIPELINE_EXECUTION))
    from ptlib import config as pt_config  # noqa: WPS433 - late import by design
    from ptlib import manifest as pt_manifest

    return {"config": pt_config, "manifest": pt_manifest}


# --------------------------------------------------------------------------- #
# Hashes and JSON
# --------------------------------------------------------------------------- #
def sha256_file(path: Path | str, chunk_size: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json(payload: Any) -> bytes:
    """Deterministic JSON bytes — the same recipe #6/#7 use for their identities."""
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def content_sha256(payload: Any) -> str:
    return sha256_bytes(canonical_json(payload))


def read_table(
    path: Path | str,
    *,
    numeric_columns: Sequence[str] = (),
    string_columns: Sequence[str] = (),
) -> "Any":
    """Read a record CSV without pandas' dtype inference on string columns.

    The `arc` env's pandas (2.3.3) segfaults when inferring a column of 64-hex
    strings (hit by the checkpoint sha256 column); reading every field as text and
    coercing the *declared* numeric columns is both safe and deterministic.
    """
    import pandas as pd

    frame = pd.read_csv(path, dtype=str)
    for column in numeric_columns:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    for column in string_columns:
        if column in frame.columns:
            frame[column] = frame[column].astype(str)
    return frame


def load_json(path: Path | str) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_text_atomic(path: Path | str, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def write_json(path: Path | str, payload: Any) -> None:
    write_text_atomic(path, json.dumps(payload, indent=2, default=str) + "\n")


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def peak_rss_mb() -> float:
    import resource

    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def mem_available_gib() -> float:
    """MemAvailable from /proc/meminfo, in GiB (0.0 if unavailable)."""
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024.0 * 1024.0)
    except (OSError, ValueError):  # pragma: no cover - non-Linux
        pass
    return 0.0


# --------------------------------------------------------------------------- #
# Frozen-input identity — byte-compatible with #6/#7's recorded `data_identity`
# --------------------------------------------------------------------------- #
def frozen_identity(data_dir: Path | str) -> dict[str, Any]:
    """Identity of the frozen input: manifest hash, locked-config hash, artifacts.

    Identical to ``06.6_cm_random/execution/common.py:frozen_identity`` (same files,
    same normalisation, same ordering); the digest this produces is what #6/#7 wrote
    into every selected-checkpoint row.
    """
    data_dir = Path(data_dir)
    manifest_file = data_dir / "manifest.json"
    locked_file = data_dir / "locked_config.json"
    if not manifest_file.exists():
        raise EvaluationError(f"frozen input manifest is missing: {manifest_file}")
    if not locked_file.exists():
        raise EvaluationError(f"frozen locked config is missing: {locked_file}")
    modules = pipeline_modules()
    artifacts = modules["manifest"].normalize_artifact_manifest(load_json(manifest_file))
    return {
        "dir": str(data_dir),
        "manifest_sha256": sha256_file(manifest_file),
        "locked_config_sha256": sha256_file(locked_file),
        "artifacts": dict(sorted(artifacts.items())),
    }


def data_identity_sha256(frozen: Mapping[str, Any]) -> str:
    """One digest for "which frozen bytes were these models trained on" (#6/#7 recipe)."""
    return content_sha256(
        {
            "manifest_sha256": frozen["manifest_sha256"],
            "locked_config_sha256": frozen["locked_config_sha256"],
            "artifacts": frozen["artifacts"],
        }
    )


def format_lr_tag(lr: float) -> str:
    """``1e-3`` / ``3e-4`` / ``1e-4`` style tag (matches #6/#7 unit naming)."""
    exponent = f"{float(lr):.0e}"
    mantissa, power = exponent.split("e")
    return f"{mantissa}e{int(power)}"


def copy_json(source: Path, target: Path) -> None:
    shutil.copyfile(source, target)
