#!/usr/bin/env python
"""Thin CLI: deterministic, target-free prediction from a saved checkpoint.

    python execution/predict.py --stage finetune --fold 0 --role test \\
        --checkpoint <run-dir>/checkpoint/model.pt --out <pred-dir>

Writes ``pred_expr.npy`` (predicted log1p expression), ``delta.npy``
(``prediction - input baseline``; c* for CM, matched control for Tahoe),
``meta.csv`` and ``manifest.json``.  It never reads treated expression, so it
cannot leak a held-out label; ``evaluate.py`` is the only truth reader.
Repeated inference is run twice and required to be byte-identical.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse  # noqa: E402

from ptlib import paths  # noqa: E402
from ptlib import predict as predictlib  # noqa: E402
from ptlib import runtime  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PRnet formal-transfer prediction (target-free, deterministic).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--stage", choices=("pretrain", "finetune"), required=True)
    parser.add_argument("--checkpoint", type=Path, required=True,
                        help="checkpoint file or a run directory containing checkpoint/model.pt")
    parser.add_argument("--out", type=Path, required=True, help="prediction output directory")
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT)
    parser.add_argument("--fold", type=int, default=None)
    parser.add_argument("--role", choices=("train", "valid", "test"), default="test")
    parser.add_argument("--device", default="cuda", choices=("cpu", "cuda", "cuda:0"))
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=2,
                        help="run inference N times and require byte-identical output")
    parser.add_argument("--hash-cache", type=Path, default=None)
    parser.add_argument("--table-cache", type=Path, default=None,
                        help="derived CSV table cache (pinned env has no parquet engine)")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ns = build_parser().parse_args(argv)
    checkpoint = ns.checkpoint
    if checkpoint.is_dir():
        checkpoint = checkpoint / "checkpoint" / "model.pt"
    options = predictlib.PredictOptions(
        stage=ns.stage,
        checkpoint=checkpoint,
        out_dir=ns.out,
        data_dir=ns.data_dir,
        role=ns.role,
        fold=ns.fold,
        device=ns.device,
        batch_size=ns.batch_size,
        repeats=ns.repeats,
        hash_cache=ns.hash_cache,
        table_cache=ns.table_cache,
        command=["predict.py", *argv],
    )
    runtime.setup_logging(options.out_dir / "run.log")
    result = predictlib.run_predict(options)
    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
