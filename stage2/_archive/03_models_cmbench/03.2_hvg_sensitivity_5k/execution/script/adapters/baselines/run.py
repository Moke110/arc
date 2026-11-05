#!/usr/bin/env python
"""Subprocess entry point for the mandatory Stage 2 baselines.

The `--baseline` flag selects one of:
`delta_mean`, `delta_mean_dose`, `ridge_fp`, `nearest_drug`.
The implementation lives in `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/baselines.py`; this file
only implements the ticket-05 CLI and checkpoint handshake.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from repo_root import EXECUTION  # noqa: E402

from baselines import fit_baseline, predict_baseline  # noqa: E402


def _fit(args: argparse.Namespace) -> int:
    config = json.loads(Path(args.config).read_text())
    checkpoint = fit_baseline(
        args.baseline,
        Path(args.train),
        Path(args.valid),
        config,
        Path(args.run_dir),
    )
    print(f"{args.baseline} fit: {checkpoint}")
    return 0


def _predict(args: argparse.Namespace) -> int:
    pred = predict_baseline(
        args.baseline,
        Path(args.cases),
        Path(args.checkpoint),
        Path(args.out),
    )
    print(f"{args.baseline} predict: {pred.shape}")
    return 0


def main(argv: list[str] | None = None) -> int:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--baseline", required=True, choices=["delta_mean", "delta_mean_dose", "ridge_fp", "nearest_drug"])

    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    fit = sub.add_parser("fit", parents=[common])
    fit.add_argument("--train", required=True)
    fit.add_argument("--valid", required=True)
    fit.add_argument("--config", required=True)
    fit.add_argument("--run-dir", required=True)
    fit.add_argument("--controls")
    fit.set_defaults(func=_fit)
    predict = sub.add_parser("predict", parents=[common])
    predict.add_argument("--cases", required=True)
    predict.add_argument("--checkpoint", required=True)
    predict.add_argument("--out", required=True)
    predict.add_argument("--controls")
    predict.set_defaults(func=_predict)
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
