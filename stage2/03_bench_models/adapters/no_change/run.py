#!/usr/bin/env python
"""CPU smoke adapter: Δ = 0 for every gene and every well.

This is the adapter the runner uses to prove the full fit → predict → score
path before any model environment is touched.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import sys

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH))

from adapter_api import Adapter, run_cli  # noqa: E402


class NoChangeAdapter(Adapter):
    name = "no_change"

    def fit(
        self,
        train: Path,
        valid: Path,
        config: dict,
        run_dir: Path,
        controls: Path | None = None,
    ) -> Path:
        checkpoint = run_dir / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        meta = {
            "model": self.name,
            "prediction": "zeros",
            "n_genes": 18164,
            "seed": config.get("seed"),
        }
        (checkpoint / "meta.json").write_text(
            json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return checkpoint

    def predict(
        self,
        cases: Path,
        checkpoint: Path,
        run_dir: Path,
        controls: Path | None = None,
    ) -> np.ndarray:
        with np.load(cases, allow_pickle=False) as npz:
            n = int(len(npz["well_id"]))
        return np.zeros((n, 18164), dtype=np.float32)


if __name__ == "__main__":
    raise SystemExit(run_cli(NoChangeAdapter()))
