#!/usr/bin/env python
"""Resolve the repository root by walking up to the directory holding `.git`.

ADR 0030 §4 and the ticket-07 contract: no Stage 2 code computes the repository
root by counting `parents[N]` or `../..` segments.  A restructure that adds or
removes one level silently breaks such arithmetic; this module makes the lookup
layout-independent and *fails loudly* when no repository is found instead of
returning a plausible-looking wrong path.

The layout anchors below come from this file's own location
(`<substep>/execution/script/repo_root.py`), so callers never count levels for
those either:

    SCRIPT     <substep>/execution/script   (benchmark code)
    EXECUTION  <substep>/execution          (configs/, models/, script/, tests/)
    BENCH      <substep>                    (e.g. 03.1_cold_drug_18k; output/)
    ROOT       repository root (the directory containing `.git`)
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["RepoRootError", "find_repo_root", "ROOT", "SCRIPT", "EXECUTION", "BENCH"]


class RepoRootError(RuntimeError):
    """No `.git` directory exists at or above the start path."""


def find_repo_root(start: Path | str | None = None) -> Path:
    """Return the closest ancestor of *start* that contains `.git`.

    Raises:
        RepoRootError: when the walk reaches the filesystem root without finding
            a `.git` directory.  Never returns a guessed path.
    """
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RepoRootError(f"no .git directory at or above {path}")


# Layout anchors, derived from this module's own location, not from the root.
SCRIPT = Path(__file__).resolve().parent
EXECUTION = SCRIPT.parent
BENCH = EXECUTION.parent

# The single repository-root lookup; importing this module asserts it succeeds.
ROOT = find_repo_root()
