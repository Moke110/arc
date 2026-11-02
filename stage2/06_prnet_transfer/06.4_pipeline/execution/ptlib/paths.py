"""Repository-root and layout anchors for the 06.4 pipeline.

ADR 0030 §4: never compute the repo root by counting ``Path.parents`` or
``../..`` levels.  Walk up from a known file until a directory contains
``.git``; fail loudly when the walk reaches the filesystem root.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    "RepoRootError",
    "find_repo_root",
    "EXECUTION",
    "PIPELINE",
    "TRANSFER_STEP",
    "ROOT",
    "DEFAULT_FROZEN_INPUT",
    "DEFAULT_OUTPUT",
]


class RepoRootError(RuntimeError):
    """No ``.git`` directory exists at or above the start path."""


def find_repo_root(start: Path | str | None = None) -> Path:
    """Return the closest ancestor of *start* that contains ``.git``."""
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RepoRootError(f"no .git directory at or above {path}")


#: ``<repo>/stage2/06_prnet_transfer/06.4_pipeline/execution/``
EXECUTION = Path(__file__).resolve().parent.parent
#: ``<repo>/stage2/06_prnet_transfer/06.4_pipeline/``
PIPELINE = EXECUTION.parent
#: ``<repo>/stage2/06_prnet_transfer/``
TRANSFER_STEP = PIPELINE.parent
#: repository root (the directory containing ``.git``)
ROOT = find_repo_root()

#: default location of the frozen input built by ticket #3
DEFAULT_FROZEN_INPUT = TRANSFER_STEP / "06.3_input" / "output"
#: default compute output root (the run tree of #5–#11 lives beside the input)
DEFAULT_OUTPUT = TRANSFER_STEP / "06.4_pipeline" / "output"


def resolve(path: Path | str | None, root: Path = ROOT) -> Path | None:
    """Resolve *path* against the repository root; ``None`` stays ``None``."""
    if path is None:
        return None
    candidate = Path(os.path.expanduser(str(path)))
    return candidate if candidate.is_absolute() else (root / candidate)
