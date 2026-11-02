#!/usr/bin/env python
"""Build the ticket-#7 queue file: 5 folds x 3 seeds x 3 LRs = 45 Tahoe-init units.

    python execution/build_queue.py [--data-dir <06.3_input/output>] \
        [--checkpoint-root <06.5 formal runs>] [--out-dir <06.7 output>]

The written ``units.json`` is the frozen definition of the matrix: it is consumed
by ``run_matrix.py`` (which hands it to the ticket-#4 ``run_queue.py``) and its
sha256 is recorded in the substep record.  Every seed's checkpoint is verified with
:mod:`warmstart` **before** the queue is written — the builder fails closed if a
checkpoint is missing, belongs to another seed, has a different gene axis or
structure, or would load only partially; a queue is never written from an
unverified initialisation.

Unit id: ``cm_f<fold>_s<seed>_lr<label>_tahoe``; args are exactly
``--fold --seed --lr --init pretrained --init-checkpoint <that seed's model.pt>
--device cuda [--hash-cache]``.  ``--out`` and ``--data-dir`` are owned by the queue.

``--profile`` writes a single-unit queue (``--profile-unit``, ``--profile-epochs``)
used to validate the wrapper on a bounded budget against the first available
formal checkpoint; that queue is marked ``complete: false`` and is never the matrix.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402
import matrix as matrix_mod  # noqa: E402
import warmstart  # noqa: E402

from ptlib import queue as queue_mod  # noqa: E402

QUEUE_SCHEMA = "cm_tahoe.queue/v1"


def _relpath(path: Path) -> str:
    try:
        return os.path.relpath(path, common.REPO_ROOT)
    except ValueError:  # pragma: no cover - different drive, not our layout
        return str(path)


def default_profile_unit(cands: list[matrix_mod.Candidate], seed: int) -> str:
    """Fold 0, highest LR, first seed — the mirror of #6's fastest candidate path."""
    fold0 = [candidate for candidate in cands if candidate.fold == 0]
    max_lr = max(candidate.lr for candidate in fold0)
    chosen = [candidate for candidate in fold0 if candidate.lr == max_lr]
    chosen.sort(key=lambda candidate: candidate.seed_index)
    if seed in {candidate.seed for candidate in chosen}:
        return next(candidate.unit_id for candidate in chosen if candidate.seed == seed)
    return chosen[0].unit_id


def build_queue(
    *,
    data_dir: Path,
    out_dir: Path,
    checkpoint_root: Path,
    queue_path: Path,
    hash_cache: str | None = None,
    seeds: list[int] | None = None,
    profile_unit: str | None = None,
    profile_epochs: int = 2,
    profile: bool = False,
) -> dict[str, Any]:
    """Verify the warm starts, then write + validate the queue; return its payload."""
    locked, _notes = common.read_locked_config(data_dir)
    spec = matrix_mod.matrix_spec(locked)
    structure = locked["structure"]
    cands = matrix_mod.candidates(spec)
    seeds = list(spec["seeds"]) if seeds is None else [int(seed) for seed in seeds]
    unknown = [seed for seed in seeds if seed not in spec["seeds"]]
    if unknown:
        raise matrix_mod.MatrixError(f"seeds {unknown} are not the frozen seeds {spec['seeds']}")

    selected: list[matrix_mod.Candidate]
    if profile:
        unit_id = profile_unit or default_profile_unit(cands, seeds[0])
        selected = [candidate for candidate in cands if candidate.unit_id == unit_id]
        if not selected:
            raise matrix_mod.MatrixError(f"--profile-unit {unit_id!r} is not a matrix unit")
        seeds = [selected[0].seed]
        if int(profile_epochs) < 1:
            raise matrix_mod.MatrixError("--profile-epochs must be positive")
    else:
        selected = [candidate for candidate in cands if candidate.seed in set(seeds)]

    # Verify the warm start(s) first: no usable checkpoint -> no queue file.
    verified = warmstart.verify_seeds(seeds, checkpoint_root, structure)
    checkpoints = verified["checkpoints"]

    units: list[dict[str, Any]] = []
    for candidate in selected:
        record = checkpoints[str(candidate.seed)]
        checkpoint = Path(record["path"])
        args = matrix_mod.train_args(candidate, checkpoint, hash_cache=hash_cache)
        if profile and int(profile_epochs) != 0:
            # bound the budget without breaking the locked min_epochs invariant
            args.extend(["--epochs", str(int(profile_epochs)), "--min-epochs", "1"])
        units.append(
            {
                "unit_id": candidate.unit_id,
                "stage": "finetune",
                "fold": candidate.fold,
                "seed": candidate.seed,
                "seed_index": candidate.seed_index,
                "lr": candidate.lr,
                "lr_label": candidate.lr_label,
                "order": candidate.order,
                "init": matrix_mod.INIT_MODE,
                "pipeline_init": matrix_mod.TRAIN_INIT_FLAG,
                "init_checkpoint": str(checkpoint),
                "init_checkpoint_worktree_relpath": _relpath(checkpoint),
                "init_checkpoint_sha256": record["sha256"],
                "init_checkpoint_config_sha256": record["checkpoint_config_sha256"],
                "init_checkpoint_epoch": record["checkpoint_epoch"],
                "args": args,
            }
        )

    payload = {
        "schema": QUEUE_SCHEMA,
        "ticket": 7,
        "issue": "Moke110/arc#7",
        "init": matrix_mod.INIT_MODE,
        "pipeline_init": matrix_mod.TRAIN_INIT_FLAG,
        "complete": bool(not profile and len(units) == spec["n_candidates"]),
        "profile": bool(profile),
        "data_dir": str(data_dir),
        "checkpoint_root": str(checkpoint_root),
        "hash_cache": hash_cache,
        "matrix": {
            "n_folds": spec["n_folds"],
            "folds": spec["folds"],
            "seeds": spec["seeds"],
            "lrs": spec["lrs"],
            "lr_labels": spec["lr_labels"],
            "n_candidates": spec["n_candidates"],
        },
        "warmstart": verified,
        "selection_rule": {
            "rule": "per fold, choose the LR maximising the mean over the three seeds of the "
                    "best validation direction AUPRC@200 (k=200)",
            "tie_break": "lower LR, then candidate order (ascending LR)",
            "test_usage": "never; the test split is scored only by ticket #9",
        },
        "units": units,
    }
    if profile:
        payload["profile_plan"] = {
            "unit_id": units[0]["unit_id"],
            "epochs": int(profile_epochs),
            "note": "bounded wrapper profile, not the formal matrix",
        }
    common.write_json_atomic(queue_path, payload)

    # Prove the definition with the pipeline's own parser/validator.
    loaded = queue_mod.load_queue(queue_path)
    if len(loaded) != len(units):
        raise SystemExit(f"queue has {len(loaded)} units, builder wrote {len(units)}")
    ids = [unit.unit_id for unit in loaded]
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate unit ids in the queue")
    for unit, expected in zip(loaded, units):
        if unit.unit_id != expected["unit_id"]:
            raise SystemExit(f"unit order/id mismatch: {unit.unit_id} != {expected['unit_id']}")
        base_args = expected["args"][: expected["args"].index("--device") + 2]
        if unit.args[: len(base_args)] != base_args:
            raise SystemExit(f"unit {unit.unit_id} args drifted: {unit.args}")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the 45-unit Tahoe-init CM candidate queue from the frozen locked config.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="frozen 06.3 input directory (default: local, else integration)")
    parser.add_argument("--checkpoint-root", type=Path, default=None,
                        help="ticket #5 formal runs root (default: local, else integration)")
    parser.add_argument("--out-dir", type=Path, default=common.DEFAULT_OUTPUT,
                        help="this substep's output directory (gitignored)")
    parser.add_argument("--queue", type=Path, default=None,
                        help="queue JSON path (default <out-dir>/units.json)")
    parser.add_argument("--hash-cache", type=Path, default=None,
                        help="frozen-artifact hash cache (default <out-dir>/hash_cache.json); "
                             "pass 'none' to disable")
    parser.add_argument("--seeds", default=None,
                        help="comma-separated frozen-seed subset (profile only)")
    parser.add_argument("--profile", action="store_true",
                        help="write a single-unit bounded profile queue instead of the matrix")
    parser.add_argument("--profile-unit", default=None,
                        help="profile unit id (default: fold 0, highest LR, first seed)")
    parser.add_argument("--profile-epochs", type=int, default=2)
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    data_dir = common.resolve_data_dir(ns.data_dir)
    checkpoint_root = common.resolve_checkpoint_root(ns.checkpoint_root)
    out_dir: Path = ns.out_dir
    if ns.profile:
        queue_path = ns.queue or (out_dir / "profile" / "profile_queue.json")
    else:
        queue_path = ns.queue or (out_dir / "units.json")
    if ns.hash_cache is None:
        hash_cache: str | None = str(out_dir / "hash_cache.json")
    elif str(ns.hash_cache).lower() in ("none", "off", ""):
        hash_cache = None
    else:
        hash_cache = str(ns.hash_cache)
    seeds = None
    if ns.seeds:
        seeds = [int(part) for part in str(ns.seeds).replace(";", ",").split(",") if part.strip()]
    try:
        payload = build_queue(
            data_dir=data_dir,
            out_dir=out_dir,
            checkpoint_root=checkpoint_root,
            queue_path=queue_path,
            hash_cache=hash_cache,
            seeds=seeds,
            profile=bool(ns.profile),
            profile_unit=ns.profile_unit,
            profile_epochs=int(ns.profile_epochs),
        )
    except (matrix_mod.MatrixError, warmstart.WarmStartError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({
        "queue": str(queue_path),
        "queue_sha256": common.sha256_file(queue_path),
        "complete": payload["complete"],
        "profile": payload["profile"],
        "data_dir": payload["data_dir"],
        "checkpoint_root": payload["checkpoint_root"],
        "hash_cache": payload["hash_cache"],
        "n_units": len(payload["units"]),
        "folds": payload["matrix"]["folds"],
        "seeds": payload["matrix"]["seeds"],
        "lrs": payload["matrix"]["lrs"],
        "warmstart": {
            "provenance_sha256": payload["warmstart"]["provenance_sha256"],
            "seeds": payload["warmstart"]["seeds"],
            "checkpoint_sha256": {
                seed: record["sha256"]
                for seed, record in payload["warmstart"]["checkpoints"].items()
            },
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
