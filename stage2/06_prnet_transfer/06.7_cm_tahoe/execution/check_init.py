#!/usr/bin/env python
"""Verify ticket #5's formal pretrain checkpoints before any ticket-#7 training.

    python execution/check_init.py \
        --data-dir <06.3_input/output> \
        --checkpoint-root <06.5_tahoe_pretrain/output/formal/runs> \
        --out <06.7 output>/init/tahoe_init.json

For every frozen seed (or the ``--seeds`` subset used by the profile) this reads the
checkpoint, validates the payload contract, the 18,164-gene axis and the structure,
runs the pipeline's strict :func:`ptlib.checkpoint.load_provenance`, applies the
anti-random-init canary, and writes one provenance record per seed.  The command
exits non-zero (and writes nothing that a queue could consume) unless every
requested seed verifies — a seed mismatch is never downgraded to a warning.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402
import matrix as matrix_mod  # noqa: E402
import warmstart  # noqa: E402

SCHEMA = "cm_tahoe.init_check/v1"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _parse_seeds(raw: str | None, spec: dict[str, Any]) -> list[int]:
    if raw is None or str(raw).strip() == "":
        return list(spec["seeds"])
    seeds = [int(part) for part in str(raw).replace(";", ",").split(",") if part.strip()]
    if not seeds:
        raise SystemExit("--seeds was given but parsed to an empty list")
    unknown = [seed for seed in seeds if seed not in spec["seeds"]]
    if unknown:
        raise SystemExit(f"--seeds {unknown} are not the frozen seeds {spec['seeds']}")
    return seeds


def run(*, data_dir: Path, checkpoint_root: Path, seeds: list[int]) -> dict[str, Any]:
    locked, notes = common.read_locked_config(data_dir)
    spec = matrix_mod.matrix_spec(locked)
    structure = locked["structure"]
    verified = warmstart.verify_seeds(seeds, checkpoint_root, structure)
    return {
        "schema": SCHEMA,
        "created_utc": now_utc(),
        "command": sys.argv,
        "git": common.pt_manifest.git_info(common.REPO_ROOT),
        "data_dir": str(data_dir),
        "locked_config_sha256": common.sha256_file(Path(data_dir) / "locked_config.json"),
        "locked_config_notes": notes,
        "structure": structure,
        "matrix": {
            "n_folds": spec["n_folds"],
            "folds": spec["folds"],
            "seeds": spec["seeds"],
            "lrs": spec["lrs"],
            "lr_labels": spec["lr_labels"],
            "n_candidates": spec["n_candidates"],
        },
        "requested_seeds": seeds,
        "all_frozen_seeds_verified": sorted(seeds) == sorted(spec["seeds"]),
        "warmstart": verified,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify the ticket-#5 pretrain checkpoints for the ticket-#7 warm start.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="frozen 06.3 input directory (default: local, else integration)")
    parser.add_argument("--checkpoint-root", type=Path, default=None,
                        help="06.5 formal runs root (default: local, else integration)")
    parser.add_argument("--out", type=Path, default=None,
                        help="provenance JSON (default <output>/init/tahoe_init.json)")
    parser.add_argument("--seeds", default=None,
                        help="comma-separated subset of the frozen seeds (profile only)")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    data_dir = common.resolve_data_dir(ns.data_dir)
    checkpoint_root = common.resolve_checkpoint_root(ns.checkpoint_root)
    locked, _notes = common.read_locked_config(data_dir)
    spec = matrix_mod.matrix_spec(locked)
    seeds = _parse_seeds(ns.seeds, spec)
    try:
        payload = run(data_dir=data_dir, checkpoint_root=checkpoint_root, seeds=seeds)
    except (warmstart.WarmStartError, matrix_mod.MatrixError, FileNotFoundError) as exc:
        print(json.dumps({"schema": SCHEMA, "all_verified": False, "error": str(exc)}, indent=2))
        return 1
    out = ns.out or (common.DEFAULT_OUTPUT / "init" / "tahoe_init.json")
    common.write_json_atomic(out, payload)
    if not ns.quiet:
        summary = {
            "init_check": str(out),
            "all_verified": bool(payload["warmstart"]["all_verified"]),
            "requested_seeds": seeds,
            "provenance_sha256": payload["warmstart"]["provenance_sha256"],
            "checkpoints": {
                seed: {
                    "path": record["path"],
                    "sha256": record["sha256"],
                    "checkpoint_config_sha256": record["checkpoint_config_sha256"],
                    "checkpoint_epoch": record["checkpoint_epoch"],
                    "gene_axis": record["gene_axis"],
                    "parameter_count": record["parameter_count"],
                    "n_loaded_tensors": record["load_provenance"]["n_loaded_tensors"],
                    "missing_keys": record["load_provenance"]["missing_keys"],
                    "unexpected_keys": record["load_provenance"]["unexpected_keys"],
                    "weights_identical": record["load_provenance"]["weights_identical"],
                    "matches_fresh_random_init": record["anti_random_init"][
                        "matches_fresh_random_init"
                    ],
                }
                for seed, record in payload["warmstart"]["checkpoints"].items()
            },
        }
        print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
