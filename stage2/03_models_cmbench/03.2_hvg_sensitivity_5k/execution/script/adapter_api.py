"""Shared adapter scaffold for the Stage 2 benchmark (ticket 05/12).

Adapters are still isolated subprocesses; this module only standardises the
two-subcommand CLI and the checkpoint handshake so tickets 13–15 do not each
re-implement argument parsing.  A real adapter subclasses `Adapter` and
implements `fit` / `predict`; the runner invokes `run.py` in the adapter's own
environment.
"""

from __future__ import annotations

import argparse
import json
import sys
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np

# Model code is vendored read-only under `pert_model/`; never leave __pycache__.
sys.dont_write_bytecode = True


class Adapter(ABC):
    """Minimal fit/predict seam from ticket 05."""

    name: str = "adapter"

    @abstractmethod
    def fit(
        self,
        train: Path,
        valid: Path,
        config: dict,
        run_dir: Path,
        controls: Path | None = None,
    ) -> Path:
        """Train on `train`, early-stop on `valid`, and return a checkpoint path."""

    @abstractmethod
    def predict(
        self,
        cases: Path,
        checkpoint: Path,
        run_dir: Path,
        controls: Path | None = None,
    ) -> np.ndarray:
        """Return a per-well Δ of shape `(n_cases, 5000)` in row order."""


def run_cli(adapter: Adapter, argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=getattr(adapter, "name", "adapter"))
    sub = parser.add_subparsers(dest="command", required=True)

    fit = sub.add_parser("fit")
    fit.add_argument("--train", required=True)
    fit.add_argument("--valid", required=True)
    fit.add_argument("--config", required=True)
    fit.add_argument("--run-dir", required=True)
    fit.add_argument("--controls")

    predict = sub.add_parser("predict")
    predict.add_argument("--cases", required=True)
    predict.add_argument("--checkpoint", required=True)
    predict.add_argument("--out", required=True)
    predict.add_argument("--controls")

    args = parser.parse_args(argv)
    if args.command == "fit":
        config = json.loads(Path(args.config).read_text())
        run_dir = Path(args.run_dir)
        checkpoint = adapter.fit(
            Path(args.train),
            Path(args.valid),
            config,
            run_dir,
            controls=Path(args.controls) if args.controls else None,
        )
        (run_dir / "checkpoint_path.txt").write_text(
            str(Path(checkpoint).resolve()) + "\n", encoding="utf-8"
        )
        print(f"{adapter.name} fit: {checkpoint}")
        return 0

    run_dir = Path(args.out).resolve().parent
    pred = adapter.predict(
        Path(args.cases),
        Path(args.checkpoint),
        run_dir,
        controls=Path(args.controls) if args.controls else None,
    )
    pred = np.asarray(pred, dtype=np.float32)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.save(args.out, pred)
    print(f"{adapter.name} predict: {pred.shape}")
    return 0
