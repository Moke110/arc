#!/usr/bin/env python
"""Thin CLI: train one PRnet formal-transfer unit (ticket #4 pipeline).

    # pretrain (Tahoe), matched controls, AUPRC@200 selection
    python execution/train.py --stage pretrain --out <run-dir> --seed 20261021

    # finetune one CM fold, fixed consensus c* at validation
    python execution/train.py --stage finetune --fold 0 --out <run-dir> \\
        --init random --lr 3e-4

    # strict warm start from a pretrain checkpoint
    python execution/train.py --stage finetune --fold 0 --out <run-dir> \\
        --init pretrained --init-checkpoint <pretrain-run>/checkpoint/model.pt

    # contract check only (no training)
    python execution/train.py --stage finetune --fold 0 --out /tmp/ignored --check

The library lives in ``ptlib/train.py``; this file only parses, logs and exits.
Tickets #5/#6/#7 invoke it (directly or through ``run_queue.py``); #10/#11 only
pass ``--structure`` / ``--train-drugs``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptlib import runtime  # noqa: E402
from ptlib import train as trainlib  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    options = trainlib.options_from_args(argv)
    runtime.setup_logging(options.out_dir / "run.log")
    if options.check:
        report = trainlib.check_run(options)
        print(json.dumps(report, indent=2, default=str))
        return 0
    result = trainlib.run_training(options)
    print(json.dumps(result.summary(), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
