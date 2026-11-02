#!/usr/bin/env python
"""Build and verify the ticket-#9 evaluation unit manifest.

    python execution/build_units.py            # the 30 units (fails closed until #7 lands)
    python execution/build_units.py --allow-incomplete   # harness dry run, random branch
    python execution/build_units.py --check    # re-verify an existing manifest

Writes ``output/record/units.json`` (the tracked manifest) and
``output/record/units_check.json`` (the verification summary). It reads only training-side
records and checkpoint bytes — it never reads treated test expression.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (  # noqa: E402
    DEFAULT_BASELINES_DIR,
    DEFAULT_DATA_DIR,
    DEFAULT_OUT,
    DEFAULT_RANDOM_RECORD,
    DEFAULT_TAHOE_RECORD,
    EvaluationError,
    write_json,
)
from units import build_unit_manifest, load_unit_manifest  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify + write the 30-unit evaluation manifest (ticket #9).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--random-record", type=Path, default=DEFAULT_RANDOM_RECORD)
    parser.add_argument("--tahoe-record", type=Path, default=DEFAULT_TAHOE_RECORD)
    parser.add_argument("--baselines-dir", type=Path, default=DEFAULT_BASELINES_DIR)
    parser.add_argument(
        "--out", type=Path, default=DEFAULT_OUT / "record" / "units.json",
        help="the tracked unit manifest",
    )
    parser.add_argument(
        "--check-json", type=Path, default=None,
        help="verification summary (default: <out dir>/units_check.json)",
    )
    parser.add_argument(
        "--allow-incomplete", action="store_true",
        help="build a partial manifest when a branch is not frozen yet (harness dry run)",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="re-verify the existing --out manifest instead of rebuilding it",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    check_path = ns.check_json or (ns.out.parent / "units_check.json")
    if ns.check:
        payload = load_unit_manifest(ns.out, require_complete=False)
        gates = payload.get("gates", {})
        summary = {
            "schema": payload.get("schema"),
            "complete": payload.get("complete"),
            "n_units": len(payload.get("units", [])),
            "data_identity": payload.get("data_identity"),
            "gates": {key: value for key, value in gates.items() if isinstance(value, bool)},
            "failed_gates": [
                key
                for key, value in gates.items()
                if isinstance(value, bool) and not value
            ],
        }
        write_json(check_path, summary)
        print(json.dumps(summary, indent=2))
        return 0 if not summary["failed_gates"] else 1

    try:
        payload = build_unit_manifest(
            data_dir=ns.data_dir,
            random_record=ns.random_record,
            tahoe_record=ns.tahoe_record if Path(ns.tahoe_record).exists() else ns.tahoe_record,
            baselines_dir=ns.baselines_dir,
            allow_incomplete=ns.allow_incomplete,
        )
    except EvaluationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if not payload["complete"] and not ns.allow_incomplete:
        print(
            "error: unit manifest is incomplete; pass --allow-incomplete for a dry run",
            file=sys.stderr,
        )
        return 2
    write_json(ns.out, payload)
    gates = payload.get("gates", {})
    summary = {
        "schema": payload["schema"],
        "complete": payload["complete"],
        "n_units": len(payload["units"]),
        "inits": {
            init: sum(1 for unit in payload["units"] if unit["init_name"] == init)
            for init in ("random", "tahoe")
        },
        "folds": sorted({int(unit["fold"]) for unit in payload["units"]}),
        "seeds": payload["requirements"]["seeds"],
        "selected_lr_by_fold_init": gates.get("lr_per_fold_init"),
        "data_identity": payload["data_identity"],
        "checkpoint_bytes_total": sum(
            int(unit["checkpoint"]["bytes"]) for unit in payload["units"]
        ),
        "gates": {key: value for key, value in gates.items() if isinstance(value, bool)},
        "failed_gates": [
            key for key, value in gates.items() if isinstance(value, bool) and not value
        ],
        "missing_inits": payload["missing_inits"],
        "manifest": str(ns.out),
    }
    write_json(check_path, summary)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
