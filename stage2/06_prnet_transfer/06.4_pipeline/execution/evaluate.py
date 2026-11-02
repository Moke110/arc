#!/usr/bin/env python
"""Thin CLI: score a prediction directory with the frozen direction-AUPRC metric.

    python execution/evaluate.py --stage finetune --fold 0 --role test \\
        --pred-dir <pred-dir> --out <eval-dir>

Reads the truth Δ from the frozen input (CM: per-plate matched DMSO mean;
Tahoe: matched ``(plate, cell_line)`` control mean), aggregates replicates within
a condition, then computes direction AUPRC@{100,200,500} and Δ RMSE / Pearson /
Spearman.  Writes ``metrics.csv``, ``per_condition.csv``, ``aggregate.json``,
``deltas.npz`` and ``manifest.json``.  This is the only component that reads
validation/test treated expression.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptlib import evaluate as evallib  # noqa: E402
from ptlib import paths  # noqa: E402
from ptlib import runtime  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PRnet formal-transfer evaluator (truth Δ from matched controls).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--stage", choices=("pretrain", "finetune"), required=True)
    parser.add_argument("--pred-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="evaluation output directory")
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT)
    parser.add_argument("--fold", type=int, default=None)
    parser.add_argument("--role", choices=("train", "valid", "test"), default="test")
    parser.add_argument("--hash-cache", type=Path, default=None)
    parser.add_argument("--table-cache", type=Path, default=None,
                        help="derived CSV table cache (pinned env has no parquet engine)")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ns = build_parser().parse_args(argv)
    options = evallib.EvaluateOptions(
        stage=ns.stage,
        pred_dir=ns.pred_dir,
        out_dir=ns.out,
        data_dir=ns.data_dir,
        role=ns.role,
        fold=ns.fold,
        hash_cache=ns.hash_cache,
        table_cache=ns.table_cache,
        command=["evaluate.py", *argv],
    )
    runtime.setup_logging(options.out_dir / "run.log")
    result = evallib.run_evaluation(options)
    print(json.dumps(result["aggregate"], indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
