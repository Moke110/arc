#!/usr/bin/env python
"""Thin CLI: resume-safe parallel queue of training units (#5/#6/#7 driver).

    python execution/run_queue.py --queue units.json --out-root <out> \\
        --data-dir <frozen input> --workers 2

Queue format (JSON)::

    {"units": [{"unit_id": "...", "stage": "pretrain|finetune",
                "args": ["--seed", "...", "--lr", "...", ...]}]}

Completed units are skipped on a re-run (status CSV + run manifest + checkpoint
sha256 checked); ``--force`` re-runs everything; ``--dry-run`` prints the plan.
Each unit runs as its own ``train.py`` subprocess with ``OMP_NUM_THREADS`` set,
so the CM candidate cap (<=4 processes) and the pretrain memory rule are obeyed
by ``--workers``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptlib import paths  # noqa: E402
from ptlib import queue as queuelib  # noqa: E402
from ptlib import runtime  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Resume-safe parallel queue for PRnet formal-transfer training units.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--queue", type=Path, required=True, help="queue JSON file")
    parser.add_argument("--out-root", type=Path, required=True,
                        help="parent of one run directory per unit")
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT)
    parser.add_argument("--status-csv", type=Path, default=None,
                        help="default <out-root>/queue_status.csv")
    parser.add_argument("--workers", type=int, default=1,
                        help=f"1..{queuelib.MAX_WORKERS} concurrent training processes")
    parser.add_argument("--omp-threads", type=int, default=2,
                        help="OMP_NUM_THREADS per worker (brief: 2 for CM candidates)")
    parser.add_argument("--force", action="store_true", help="re-run even completed units")
    parser.add_argument("--dry-run", action="store_true", help="print the plan and exit")
    parser.add_argument("--python", default=sys.executable,
                        help="interpreter for the unit subprocesses")
    parser.add_argument("--table-cache", type=Path, default=None,
                        help="derived CSV table cache (pinned env has no parquet engine)")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ns = build_parser().parse_args(argv)
    runtime.setup_logging(None)
    options = queuelib.QueueOptions(
        queue=ns.queue,
        out_root=ns.out_root,
        data_dir=ns.data_dir,
        status_csv=ns.status_csv,
        workers=ns.workers,
        force=ns.force,
        dry_run=ns.dry_run,
        omp_threads=ns.omp_threads,
        table_cache=ns.table_cache,
        python=ns.python,
        train_script=Path(__file__).resolve().parent / "train.py",
    )
    summary = queuelib.run_queue(options)
    print(json.dumps(summary, indent=2, default=str))
    return 1 if summary.get("failed") else 0


if __name__ == "__main__":
    raise SystemExit(main())
