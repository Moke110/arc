#!/usr/bin/env python
"""Independent acceptance checks for the ticket-16 benchmark run.

Reads only the per-unit artifacts and the generated leaderboard; does not
import any model.  Fails closed on a missing/NaN prediction, a wrong group
count, a non-completed unit, or a leaderboard row that is not backed by five
completed folds.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from repo_root import BENCH, EXECUTION  # noqa: E402

OUTPUT = BENCH / "output" / "units"
RESULTS = BENCH / "output" / "results"
REGISTRY = json.loads((EXECUTION / "configs/models.json").read_text())
N_GENES = int(json.loads((BENCH / "output/input/manifest.json").read_text())["n_genes"])
N_TEST_GROUPS = 276


def _checkpoint_path(run_dir: Path) -> Path:
    """Resolve a unit checkpoint the way `runner._checkpoint_path` does.

    The recorded `checkpoint_path.txt` may hold a stale absolute path from an
    earlier tree layout; the runner already falls back to the run's default
    `checkpoint/` directory, and this validator must use the same rule.
    """
    marker = run_dir / "checkpoint_path.txt"
    if marker.exists():
        recorded = Path(marker.read_text().strip())
        if recorded.exists():
            return recorded
    return run_dir / "checkpoint"


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> int:
    status_path = OUTPUT / "run_status.csv"
    if not status_path.exists():
        print(f"missing {status_path}", file=sys.stderr)
        return 2
    status = pd.read_csv(status_path)
    expected = len(REGISTRY) * 5
    problems: list[str] = []
    if len(status) != expected:
        problems.append(f"run_status has {len(status)} rows, expected {expected}")
    bad = status[status["status"] != "completed"]
    if len(bad):
        problems.append(f"non-completed units: {bad[['model', 'fold', 'status']].to_dict('records')}")

    for _, unit in status.iterrows():
        model, fold = str(unit["model"]), int(unit["fold"])
        run_dir = OUTPUT / model / f"fold_{fold}"
        pred = np.load(run_dir / "pred_test.npy", mmap_mode="r")
        if pred.shape[1] != N_GENES:
            problems.append(f"{model}/{fold}: pred_test width {pred.shape[1]}")
        if not np.isfinite(np.asarray(pred, dtype=np.float32)).all():
            problems.append(f"{model}/{fold}: pred_test non-finite")
        if pred.shape[0] != int(unit["n_test"]):
            problems.append(f"{model}/{fold}: pred_test rows {pred.shape[0]} != {unit['n_test']}")
        del pred
        pred = np.load(run_dir / "pred_valid.npy", mmap_mode="r")
        if pred.shape[1] != N_GENES or pred.shape[0] != int(unit["n_valid"]):
            problems.append(f"{model}/{fold}: pred_valid shape {pred.shape}")
        del pred
        groups = _json(run_dir / "group_ids.json")["test"]
        if len(set(groups)) != N_TEST_GROUPS:
            problems.append(f"{model}/{fold}: {len(set(groups))} test groups, expected {N_TEST_GROUPS}")
        metrics = _json(run_dir / "metrics.json")
        if not np.isfinite(float(metrics["auprc"]["200"])):
            problems.append(f"{model}/{fold}: non-finite AUPRC")
        if not (run_dir / "checkpoint_path.txt").exists():
            problems.append(f"{model}/{fold}: missing checkpoint_path.txt")
        elif not _checkpoint_path(run_dir).exists():
            problems.append(f"{model}/{fold}: checkpoint path does not exist")

    leaderboard = pd.read_csv(RESULTS / "leaderboard.csv")
    if len(leaderboard) != len(REGISTRY):
        problems.append(f"leaderboard has {len(leaderboard)} rows, expected {len(REGISTRY)}")
    if (leaderboard["n_folds"] != 5).any() or not leaderboard["winner_eligible"].all():
        problems.append("leaderboard has a row without five completed folds")
    if not np.isfinite(leaderboard["headline_auprc"]).all():
        problems.append("leaderboard has a non-finite headline")

    print(f"units checked: {len(status)} / {expected}")
    print(f"leaderboard rows: {len(leaderboard)}")
    if problems:
        print("\nFAILED:")
        for problem in problems:
            print(" -", problem)
        return 1
    print("all acceptance checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
