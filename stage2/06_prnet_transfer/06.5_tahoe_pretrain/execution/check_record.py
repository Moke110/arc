#!/usr/bin/env python
"""Check a 06.5 pretrain record: structure, budget, digests, curve, resources.

    python check_record.py --record output/formal/tahoe_pretrain_record.json
    python check_record.py --record output/profile/profile_record.json --expect-mode profile
    python check_record.py --record output/formal/tahoe_pretrain_record.json --json

The checker is strict by default: it re-reads the frozen locked config, the
frozen input manifest, every run's own ``manifest.json`` / ``fit_history.json``
and (when present) the checkpoint files, and it refuses a formal record unless
it holds exactly the three frozen seeds, trained with the locked budget.
``--expect-mode profile`` accepts a short profile record but never lets it stand
in for the formal one.  Exit status 0 = pass, 1 = at least one error, 2 = the
record could not be read.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

EXECUTION = Path(__file__).resolve().parent
sys.path.insert(0, str(EXECUTION))

import record  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate a 06.5 Tahoe-pretrain record against the frozen contract.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--record", type=Path, required=True, help="record JSON to check")
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="frozen 06.3 input (default: the data_dir recorded in the record)")
    parser.add_argument("--expect-mode", choices=("profile", "formal", "any"), default="any",
                        help="require the record's mode (any = use the record's own mode)")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="accept a partial formal sweep (a per-seed run before the full sweep)")
    parser.add_argument("--no-frozen-verify", action="store_true",
                        help="do not re-read the frozen input manifest (offline record inspection)")
    parser.add_argument("--json", action="store_true", help="emit the full report as JSON")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = record.read_record(args.record)
    except record.RecordError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    identity = payload.get("input") or {}
    data_dir = Path(args.data_dir) if args.data_dir else Path(str(identity.get("data_dir", "")))
    locked = None
    if data_dir and (data_dir / "locked_config.json").exists():
        try:
            raw_locked = record.read_json(data_dir / "locked_config.json")
            locked = record.normalize_locked_config(raw_locked)
        except record.RecordError as exc:
            print(f"FAIL: cannot read the locked config: {exc}", file=sys.stderr)
            return 2
    else:
        print("WARNING: frozen locked config not found; skipping the locked-budget comparison",
              file=sys.stderr)
    expect_mode = None if args.expect_mode == "any" else args.expect_mode
    errors, warnings = record.check_record(
        payload,
        locked=locked,
        expect_mode=expect_mode,
        verify_frozen=not args.no_frozen_verify,
        expect_complete=not args.allow_incomplete,
    )
    summary = payload.get("summary") or {}
    if args.json:
        print(json.dumps({"record": str(args.record), "errors": errors, "warnings": warnings,
                          "summary": summary}, indent=2, default=str))
        return 0 if not errors else 1
    print(f"record: {args.record}")
    print(f"  mode={payload.get('mode')} complete={payload.get('complete')} "
          f"units={len(payload.get('units') or [])}")
    print(f"  seeds planned={payload.get('seeds_planned')} "
          f"completed={payload.get('seeds_completed')} failed={payload.get('seeds_failed')}")
    for unit in payload.get("units") or []:
        resources = unit.get("resources") or {}
        timings = unit.get("timings") or {}
        print(
            f"  - {unit.get('unit_id')}: {unit.get('status')} "
            f"epochs={unit.get('epochs_run')} stop={unit.get('stop_reason')} "
            f"best_valid_AUPRC@200={unit.get('best_valid_auprc')} "
            f"wall={timings.get('wall_s_pipeline')}s "
            f"host_peak={resources.get('host_peak_rss_mb_pipeline')}MiB "
            f"gpu_peak={resources.get('gpu_peak_mib_process') or resources.get('gpu_peak_mib_device')}MiB"
        )
    for warning in warnings:
        print(f"  warning: {warning}")
    if errors:
        print(f"FAIL: {len(errors)} error(s)")
        for error in errors:
            print(f"  - {error}")
        return 1
    print("PASS: the record is internally consistent and matches the frozen contract")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
