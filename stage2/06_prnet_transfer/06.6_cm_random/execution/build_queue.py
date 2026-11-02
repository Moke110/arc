#!/usr/bin/env python
"""Build the ticket-#6 queue file: 5 folds x 3 seeds x 3 LRs = 45 random-init units.

    python execution/build_queue.py [--data-dir <06.3_input/output>] [--out-dir <06.6 output>]

The written ``units.json`` is the frozen definition of the matrix: it is consumed
by ``run_matrix.py`` (which hands it to the ticket-#4 ``run_queue.py``) and its
sha256 is recorded in the substep record.  The unit definition is validated with
the pipeline's **own** ``ptlib.queue.load_queue`` (not a second parser), so the
queue rejects anything this file would get wrong.

Unit id: ``cm_f<fold>_s<seed>_lr<label>_random``; args are exactly
``--fold --seed --lr --init random --device cuda [--hash-cache]``.  ``--out`` and
``--data-dir`` are owned by the queue.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402
import matrix as matrix_mod  # noqa: E402

from ptlib import queue as queue_mod  # noqa: E402

QUEUE_SCHEMA = "cm_random.queue/v1"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the 45-unit random-init CM candidate queue from the frozen locked config.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=common.DEFAULT_DATA_DIR,
                        help="frozen 06.3 input directory")
    parser.add_argument("--out-dir", type=Path, default=common.DEFAULT_OUTPUT,
                        help="this substep's output directory (gitignored)")
    parser.add_argument("--queue", type=Path, default=None,
                        help="queue JSON path (default <out-dir>/units.json)")
    parser.add_argument("--hash-cache", type=Path, default=None,
                        help="frozen-artifact hash cache (default <out-dir>/hash_cache.json); "
                             "pass 'none' to disable")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    data_dir: Path = ns.data_dir
    out_dir: Path = ns.out_dir
    queue_path: Path = ns.queue or (out_dir / "units.json")
    if ns.hash_cache is None:
        hash_cache: str | None = str(out_dir / "hash_cache.json")
    elif str(ns.hash_cache).lower() in ("none", "off", ""):
        hash_cache = None
    else:
        hash_cache = str(ns.hash_cache)

    locked, _notes = common.read_locked_config(data_dir)
    spec = matrix_mod.matrix_spec(locked)
    cands = matrix_mod.candidates(spec)

    units = [
        {
            "unit_id": candidate.unit_id,
            "stage": "finetune",
            "fold": candidate.fold,
            "seed": candidate.seed,
            "lr": candidate.lr,
            "lr_label": candidate.lr_label,
            "init": matrix_mod.INIT_MODE,
            "args": matrix_mod.train_args(candidate, hash_cache=hash_cache),
        }
        for candidate in cands
    ]
    payload = {
        "schema": QUEUE_SCHEMA,
        "ticket": 6,
        "issue": "Moke110/arc#6",
        "init": matrix_mod.INIT_MODE,
        "data_dir": str(data_dir),
        "hash_cache": hash_cache,
        "matrix": {
            "n_folds": spec["n_folds"],
            "folds": spec["folds"],
            "seeds": spec["seeds"],
            "lrs": spec["lrs"],
            "lr_labels": spec["lr_labels"],
            "n_candidates": spec["n_candidates"],
        },
        "selection_rule": {
            "rule": "per fold, choose the LR maximising the mean over the three seeds of the "
                    "best validation direction AUPRC@200 (k=200)",
            "tie_break": "lower LR, then candidate order (ascending LR)",
            "test_usage": "never; the test split is scored only by ticket #9",
        },
        "units": units,
    }
    common.write_json_atomic(queue_path, payload)

    # Prove the definition with the pipeline's own parser/validator.
    loaded = queue_mod.load_queue(queue_path)
    if len(loaded) != spec["n_candidates"]:
        raise SystemExit(f"queue has {len(loaded)} units, matrix says {spec['n_candidates']}")
    ids = [unit.unit_id for unit in loaded]
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate unit ids in the queue")
    for candidate, unit in zip(cands, loaded):
        if unit.unit_id != candidate.unit_id:
            raise SystemExit(f"unit order/id mismatch: {unit.unit_id} != {candidate.unit_id}")
        expected = [
            "--fold", str(candidate.fold),
            "--seed", str(candidate.seed),
            "--lr", repr(float(candidate.lr)),
            "--init", matrix_mod.INIT_MODE,
            "--device", "cuda",
        ]
        if unit.args[: len(expected)] != expected:
            raise SystemExit(f"unit {unit.unit_id} args drifted: {unit.args}")

    print(json.dumps({
        "queue": str(queue_path),
        "queue_sha256": common.sha256_file(queue_path),
        "data_dir": str(data_dir),
        "hash_cache": hash_cache,
        "n_candidates": spec["n_candidates"],
        "folds": spec["folds"],
        "seeds": spec["seeds"],
        "lrs": spec["lrs"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
