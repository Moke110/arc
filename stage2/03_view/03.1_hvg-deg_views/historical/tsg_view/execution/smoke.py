#!/usr/bin/env python3
"""Bounded smoke fit for one view: train, replay, score and replay again.

Each phase is a separate process (``learner.py``) so the exact fresh-process replay rule
is exercised rather than asserted. The driver writes a small receipt outside the run
directory so the run manifest stays complete.
"""
import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from views import SOURCES, VIEWS, repo_root, write_json  # noqa: E402

SCHEMA = "tsg-view-smoke-v1"


def default_source(root):
    return Path(root) / "pert_model/TranSiGen/source"


def smoke(*, view, source, freeze, output, rows, epochs, device="cpu", fold=0, seed=42,
          repo_root_path=None, sealed=None, source_dir=None, record=None):
    if view not in VIEWS or source not in SOURCES:
        raise ValueError(f"unknown smoke identity: {view}/{source}")
    if not 2 <= int(rows) <= 64:
        raise ValueError("smoke rows must be in 2..64")
    if not 1 <= int(epochs) <= 3:
        raise ValueError("smoke epochs must be in 1..3")
    if device not in ("cuda", "cpu"):
        raise ValueError("device must be cuda or cpu")
    root = repo_root(repo_root_path)
    native = Path(source_dir) if source_dir else default_source(root)
    if not (native / "src/model.py").is_file():
        raise ValueError(f"the TranSiGen source is unavailable: {native}")
    freeze = Path(freeze)
    if not (freeze / "manifest.json").is_file():
        raise ValueError(f"the view freeze is unavailable: {freeze}")
    run = Path(output) / source / view / f"fold-{fold}" / f"seed-{seed}"
    if run.exists() and any(run.iterdir()):
        raise ValueError(f"the smoke run directory is not empty: {run}")
    run.mkdir(parents=True, exist_ok=True)
    worker = HERE / "learner.py"
    common = [sys.executable, str(worker), str(freeze), str(native), str(run)]
    if sealed is not None:
        common += ["--sealed", str(sealed)]
    if repo_root_path is not None:
        common += ["--repo-root", str(repo_root_path)]
    phases = []
    started = time.monotonic()
    for phase, extra in (
            ("train", ["--rows", str(int(rows)), "--epochs", str(int(epochs)),
                       "--device", device, "--fold", str(int(fold)),
                       "--seed", str(int(seed))]),
            ("replay", []),
            ("test", []),
            ("replay", [])):
        phase_started = time.monotonic()
        completed = subprocess.run([*common, phase, *extra], check=True)
        phases.append({"phase": phase,
                       "elapsed_seconds": round(time.monotonic() - phase_started, 3),
                       "returncode": completed.returncode})
    terminal = json.loads((run / "terminal.json").read_text())
    receipt = {
        "schema": SCHEMA,
        "view": view,
        "source": source,
        "device": device,
        "seed": int(seed),
        "fold": int(fold),
        "rows": int(rows),
        "epochs": int(epochs),
        "status": terminal["status"],
        "best_epoch": terminal["best_epoch"],
        "epochs_run": terminal["epochs"],
        "validation_mse": terminal["validation_mse"],
        "training_elapsed_seconds": terminal["elapsed_seconds"],
        "training_peak_rss_kib": terminal["peak_rss_kib"],
        "driver_elapsed_seconds": round(time.monotonic() - started, 3),
        "children_peak_rss_kib": int(resource.getrusage(
            resource.RUSAGE_CHILDREN).ru_maxrss),
        "phases": phases,
        "run": str(run),
        "interpreter": sys.executable,
        "replay": "train -> replay -> test -> replay, one process per phase",
    }
    if record is not None:
        write_json(record, receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--view", required=True, choices=VIEWS)
    parser.add_argument("--source", required=True, choices=SOURCES)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "checkpoints/smoke")
    parser.add_argument("--rows", type=int, choices=range(2, 65), metavar="2..64",
                        default=32)
    parser.add_argument("--epochs", type=int, choices=range(1, 4), metavar="1..3", default=2)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--fold", type=int, default=0, choices=range(3))
    parser.add_argument("--seed", type=int, default=42, choices=(42, 43))
    parser.add_argument("--sealed", type=Path)
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--source-path", type=Path)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    receipt = smoke(view=args.view, source=args.source, freeze=args.freeze,
                    output=args.output, rows=args.rows, epochs=args.epochs,
                    device=args.device, fold=args.fold, seed=args.seed,
                    repo_root_path=args.repo_root, sealed=args.sealed,
                    source_dir=args.source_path, record=args.record)
    sys.stdout.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
