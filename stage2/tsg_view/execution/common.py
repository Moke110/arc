"""Single implementations of the small helpers every tsg_view module needs.

Each of these used to exist as a near-copy in two to four modules, which is how a
serializer, a path form or a per-source seed list can slowly drift apart. This module is
imported (never copied) by ``views``, ``matrix``, ``learner``, ``report``, ``render`` and
``truth``; ``views`` re-exports the names its callers already import from it.

Nothing here binds the Step 01 API, reads a payload or imports the fitting stack: it is
safe to import from a test that only inspects the schedule or the released constants.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

#: The published source order (tenv first; the report's tenv pass may publish alone).
SOURCE_ORDER = ("tenv", "tahoe_bulk")
#: The per-source view schedule, in publication order. ``solo-deepest`` is tenv-only.
VIEW_ORDER = {
    "tenv": ("full", "merge", "solo", "full-pooled-label", "solo-deepest"),
    "tahoe_bulk": ("full", "merge", "solo"),
}
#: Every main fit runs at both training seeds; ``solo-deepest`` is seed 42 only.
SEEDS = (42, 43)
VIEW_SEEDS = {"solo-deepest": (42,)}

#: Released-population constants (reconnaissance notes-bases-views §0/§5) and the
#: member-count histogram they imply. The three consumers assert the *observed* counts
#: against these values, so a released-input change fails closed instead of silently
#: publishing a differently shaped view/truth product.
RELEASED = {
    "tenv": {
        "rows": 4478, "treated": 4206, "vehicle": 272, "conditions": 1377,
        "member_histogram": {"1": 0, "2": 0, "3": 1320, "4": 48, "5": 0, "6": 9},
    },
    "tahoe_bulk": {
        "rows": 66718, "treated": 65318, "vehicle": 1400, "conditions": 52777,
        "member_histogram": {"1": 43036, "2": 7941, "3": 1400, "4": 350, "5": 0,
                             "6": 0, "7": 0, "8": 0, "9": 0, "10": 0, "11": 0,
                             "12": 0, "13": 0, "14": 0, "15": 0, "16": 50},
    },
}


def views_for(source):
    if source not in VIEW_ORDER:
        raise ValueError(f"unknown source: {source}")
    return VIEW_ORDER[source]


def seeds_for(view):
    return VIEW_SEEDS.get(view, SEEDS)


def released_conditions():
    return {source: value["conditions"] for source, value in RELEASED.items()}


def released_histogram(source):
    """Member-count histogram padded with the zero buckets the release carries."""
    return dict(RELEASED[source]["member_histogram"])


def replicated_conditions():
    """Conditions with at least two members (the merge row space)."""
    return {source: value["conditions"] - value["member_histogram"]["1"]
            for source, value in RELEASED.items()}


def repo_root(explicit=None):
    """The checkout that owns the pinned inputs and payloads.

    Never hard-code an absolute path: the repository is found by walking up from this
    module until a directory carrying ``.git`` is reached.
    """
    if explicit is not None:
        root = Path(explicit).resolve()
        if not (root / ".git").exists():
            raise ValueError(f"not a repository root: {root}")
        return root
    for path in Path(__file__).resolve().parents:
        if (path / ".git").exists():
            return path
    raise ValueError("cannot locate the repository root; pass --repo-root")


def relative(path, root):
    """Repo-relative posix path, else the bare file name (never an absolute path)."""
    try:
        return Path(path).absolute().relative_to(Path(root).absolute()).as_posix()
    except (ValueError, OSError):
        return Path(path).name


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dumps(value):
    """The canonical JSON text (no trailing newline); persisted files add one."""
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)


def write_json(path, value):
    """The repository's canonical deterministic JSON form.

    ``allow_nan=False`` matters: a non-finite score is a defect, and the literal ``NaN``
    is not JSON. Every caller writes finite products, so the bytes are unchanged.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps(value) + "\n")
    return path


def read_json(path):
    return json.loads(Path(path).read_text())
