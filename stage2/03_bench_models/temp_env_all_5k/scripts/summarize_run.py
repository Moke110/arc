#!/usr/bin/env python
"""Ticket-16 run cards: aggregate wall-clock, peak GPU, and final epoch/step.

Reads `output/run_status.csv` plus each `output/<model>/fold_<k>/` unit and
writes `results/run_cards.csv`.  Peak GPU comes from (in order) the adapter's
`fit_history.json`, a `peak GPU` line in `run.log`, or the external
`output/gpu_watch.csv` sample window.  No model internals are imported.
"""

from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cases import BenchData  # noqa: E402

ARCHIVE = HERE.parent
OUTPUT = ARCHIVE / "output"
RESULTS = ARCHIVE / "results"
DATA = ARCHIVE / "data"

COLUMNS = [
    "model", "fold", "status", "wall_clock_s", "fit_wall_s",
    "predict_valid_s", "predict_test_s", "peak_rss_mb", "peak_gpu_mb",
    "final_epoch_or_step", "best_epoch_or_step", "headline_auprc",
    "n_train", "n_valid", "n_test", "error",
]


def _json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def _peak_from_log(path: Path) -> float:
    if not path.exists():
        return float("nan")
    text = path.read_text(errors="ignore")
    values = [float(m) for m in re.findall(r"peak GPU\s+([0-9.]+)\s*MB", text)]
    return max(values) if values else float("nan")


def _epochs(model: str, history: dict) -> tuple:
    """Return `(final_epoch_or_step, best_epoch_or_step)` for one model."""

    if model in {"prnet", "chemcpa"}:
        rows = history.get("history", [])
        if not rows:
            return (np.nan, np.nan)
        valid = [float(r.get("valid_loss", np.nan)) for r in rows]
        best = int(np.nanargmin(valid)) if np.isfinite(valid).any() else -1
        return (int(rows[-1].get("epoch", -1)) + 1, int(rows[best].get("epoch", -1)) + 1)
    if model == "cpa":
        text = history.get("epoch_history_csv")
        if not text:
            return (np.nan, np.nan)
        frame = pd.read_csv(pd.io.common.StringIO(text))
        if "epoch" not in frame:
            return (np.nan, np.nan)
        final = int(frame["epoch"].max()) + 1
        valid = frame[frame.get("mode") == "valid"] if "mode" in frame else frame
        if len(valid) and "recon_loss" in valid:
            best = int(valid.loc[valid["recon_loss"].idxmin(), "epoch"]) + 1
        else:
            best = np.nan
        return (final, best)
    if model == "biolord":
        rows = history.get("history", [])
        best = history.get("best", {}).get("epoch", np.nan)
        final = int(rows[-1]["epoch"]) if rows else np.nan
        return (final, best)
    if model == "perturbnet":
        rows = history.get("cinn", [])
        return (
            int(rows[-1]["epoch"]) + 1 if rows else np.nan,
            int(history.get("best_epoch", -1)) + 1,
        )
    if model == "squidiff":
        return (history.get("steps_run", np.nan), history.get("best", {}).get("step", np.nan))
    if model == "cmonge":
        rows = history.get("history", [])
        final = int(rows[-1]["iter"]) if rows else np.nan
        if rows:
            totals = [float(r.get("total", np.nan)) for r in rows]
            best = int(rows[int(np.nanargmin(totals))]["iter"]) if np.isfinite(totals).any() else np.nan
        else:
            best = np.nan
        return (final, best)
    return (np.nan, np.nan)


def _gpu_window(started: str, finished: str, samples: pd.DataFrame) -> float:
    if samples.empty or not started or not finished:
        return float("nan")
    try:
        start = datetime.fromisoformat(started).timestamp()
        stop = datetime.fromisoformat(finished).timestamp()
    except Exception:
        return float("nan")
    window = samples[(samples["timestamp_s"] >= start - 5) & (samples["timestamp_s"] <= stop + 5)]
    return float(window["total_used_mb"].max()) if len(window) else float("nan")


def main() -> int:
    status_path = OUTPUT / "run_status.csv"
    if not status_path.exists():
        print(f"missing {status_path}; run the queue first", file=sys.stderr)
        return 2
    status = pd.read_csv(status_path)
    samples = pd.read_csv(OUTPUT / "gpu_watch.csv") if (OUTPUT / "gpu_watch.csv").exists() else pd.DataFrame()
    data = BenchData(DATA)

    rows = []
    for _, unit in status.iterrows():
        model, fold = str(unit["model"]), int(unit["fold"])
        run_dir = OUTPUT / model / f"fold_{fold}"
        history = _json(run_dir / "fit_history.json")
        detail = _json(run_dir / "status.json")
        metrics = _json(run_dir / "metrics.json")
        final, best = _epochs(model, history)
        # Prefer the external sampler: it is the same measurement (total device
        # memory, including the ~2.5 GB desktop baseline) for every unit.  The
        # adapters' own `peak_gpu_mb` is a process-allocation figure and is only
        # a fallback for units shorter than the 30 s sampling interval.
        gpu = _gpu_window(detail.get("started_at"), detail.get("finished_at"), samples)
        if not np.isfinite(float(gpu)):
            gpu = history.get("peak_gpu_mb", np.nan)
        if not np.isfinite(float(gpu)):
            gpu = _peak_from_log(run_dir / "run.log")
        rows.append(
            {
                "model": model,
                "fold": fold,
                "status": unit.get("status"),
                "wall_clock_s": unit.get("wall_clock_s"),
                "fit_wall_s": detail.get("fit", {}).get("wall_s"),
                "predict_valid_s": detail.get("predict_valid", {}).get("wall_s"),
                "predict_test_s": detail.get("predict_test", {}).get("wall_s"),
                "peak_rss_mb": unit.get("peak_rss_mb"),
                "peak_gpu_mb": gpu,
                "final_epoch_or_step": final,
                "best_epoch_or_step": best,
                "headline_auprc": metrics.get("auprc", {}).get("200", np.nan),
                "n_train": unit.get("n_train"),
                "n_valid": unit.get("n_valid"),
                "n_test": unit.get("n_test"),
                "error": (str(unit.get("error"))[:300] if pd.notna(unit.get("error")) else ""),
            }
        )

    frame = pd.DataFrame(rows, columns=COLUMNS)
    frame = frame.sort_values(["model", "fold"], kind="mergesort").reset_index(drop=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    frame.to_csv(RESULTS / "run_cards.csv", index=False, lineterminator="\n")
    print(frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
