#!/usr/bin/env python
"""Phase-B real-data smoke: the ticket #4 acceptance evidence.

Runs the whole path on the **real frozen input** at small scale, on the pinned GPU
env, as separate CLI subprocesses (so each run records its own manifest):

1. ``--stage pretrain`` 2 epochs on a frozen Tahoe subset (random init), with
   per-epoch timing and peak RSS;
2. ``--stage finetune`` 2 epochs on CM fold 0, **random init**;
3. ``--stage finetune`` 2 epochs on CM fold 0, **warm start** from (1);
4. predict + evaluate for both CM initialisations (valid and test) and for the
   pretrain checkpoint (valid);
5. regenerate model cards with the evaluation numbers attached;
6. write ``output/smoke/smoke_record.json`` + ``smoke_summary.csv``.

Everything printed here is a **pipeline number**, not a science result: 2 epochs,
small row subsets, one fold.

    make gpu-smoke                      # from execution/
    python smoke_real.py --data-dir <06.3_input/output> [--build-cache]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptlib import paths  # noqa: E402

EXECUTION = Path(__file__).resolve().parent
SMOKE_OUT = paths.PIPELINE / "output" / "smoke"


def _run(command: list[str], log_path: Path, env: dict[str, str]) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.flush()
        code = subprocess.run(command, env=env, stdout=handle, stderr=subprocess.STDOUT).returncode
    elapsed = time.perf_counter() - started
    print(f"  [{'ok' if code == 0 else 'FAIL'}] {Path(command[1]).stem} "
          f"{' '.join(command[2:5])} ({elapsed:.1f}s) -> {log_path.name}")
    return code


def _manifest(run_dir: Path) -> dict:
    path = run_dir / "manifest.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def _history(run_dir: Path) -> dict:
    path = run_dir / "fit_history.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def _aggregate(eval_dir: Path) -> dict:
    path = eval_dir / "aggregate.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase-B real-data smoke for the 06.4 pipeline (pipeline evidence, not science).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT)
    parser.add_argument("--out-root", type=Path, default=SMOKE_OUT)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--cm-lr", type=float, default=3e-4)
    parser.add_argument("--pretrain-rows", type=int, default=3000)
    parser.add_argument("--valid-rows", type=int, default=1000)
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--tahoe-cache", type=Path, default=None)
    parser.add_argument("--build-cache", action="store_true",
                        help="build the derived 4.9 GB log1p cache first (recommended for #5)")
    parser.add_argument("--concurrency-probe", action="store_true",
                        help="run two 1-epoch pretrains concurrently and record both peak RSS")
    parser.add_argument("--warm-rerun", action="store_true",
                        help="run one extra 1-epoch pretrain to measure a fully warm-cache epoch")
    args = parser.parse_args(argv)

    data_dir = Path(args.data_dir)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    hash_cache = out_root / "hash_cache.json"
    cache = Path(args.tahoe_cache) if args.tahoe_cache else (
        paths.PIPELINE / "output" / "cache" / "tahoe_log1p_cp10k_v1.npy"
    )
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["OMP_NUM_THREADS"] = "4"
    train = str(EXECUTION / "train.py")
    predict = str(EXECUTION / "predict.py")
    evaluate = str(EXECUTION / "evaluate.py")
    card = str(EXECUTION / "make_model_card.py")
    common = ["--data-dir", str(data_dir), "--hash-cache", str(hash_cache), "--device", args.device]
    record: dict = {
        "label": "PIPELINE SMOKE — not a science result",
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "data_dir": str(data_dir),
        "device": args.device,
        "seed": int(args.seed),
        "fold": int(args.fold),
        "epochs": int(args.epochs),
        "pretrain_rows": int(args.pretrain_rows),
        "valid_rows": int(args.valid_rows),
        "steps": {},
        "errors": [],
    }
    started_all = time.perf_counter()

    def must(code: int, what: str) -> bool:
        if code != 0:
            record["errors"].append(what)
            print(f"smoke FAILED at: {what}")
            return False
        return True

    # 0) optional derived-cache build + verification (one-off for #5)
    if args.build_cache:
        print("building derived Tahoe log1p cache (one-off) …")
        code = _run(
            [args.python, train, "--stage", "pretrain", "--out", str(out_root / "_cache_build_check"),
             "--check", "--cache-mode", "build", "--tahoe-cache", str(cache), *common],
            out_root / "logs" / "cache_build.log", env,
        )
        if not must(code, "cache build"):
            return 1
        record["steps"]["cache_build"] = {"path": str(cache), "bytes": cache.stat().st_size}
        sidecar = json.loads(Path(str(cache) + ".json").read_text(encoding="utf-8"))
        record["steps"]["cache_build"].update(
            {
                "layer_sha256": sidecar.get("layer_sha256"),
                "sample_max_abs_diff": sidecar.get("sample_max_abs_diff"),
                "build_seconds": sidecar.get("build_seconds"),
            }
        )

    # 1) Tahoe pretrain, 2 epochs, subset
    pretrain_dir = out_root / f"pretrain_seed{args.seed}"
    code = _run(
        [args.python, train, "--stage", "pretrain", "--out", str(pretrain_dir),
         "--seed", str(args.seed), "--smoke", "--epochs", str(args.epochs),
         "--max-train-rows", str(args.pretrain_rows),
         "--max-valid-rows", str(args.valid_rows),
         "--tahoe-cache", str(cache), "--cache-mode", "auto", *common],
        out_root / "logs" / "pretrain.log", env,
    )
    if not must(code, "pretrain"):
        return 1
    pretrain_manifest = _manifest(pretrain_dir)
    pretrain_history = _history(pretrain_dir).get("history", [])
    record["steps"]["pretrain"] = {
        "run_dir": str(pretrain_dir),
        "checkpoint": (pretrain_manifest.get("checkpoint") or {}).get("path"),
        "checkpoint_sha256": (pretrain_manifest.get("checkpoint") or {}).get("sha256"),
        "peak_rss_mb": pretrain_manifest.get("peak_rss_mb"),
        "wall_s": (pretrain_manifest.get("timings") or {}).get("wall_s"),
        "per_epoch_seconds": [entry.get("seconds") for entry in pretrain_history],
        "per_epoch_valid_seconds": [entry.get("valid_seconds") for entry in pretrain_history],
        "per_epoch_valid_split": {
            "gather_s": [entry.get("valid_gather_seconds") for entry in pretrain_history],
            "forward_s": [entry.get("valid_forward_seconds") for entry in pretrain_history],
            "auprc_s": [entry.get("valid_metrics_seconds") for entry in pretrain_history],
        },
        "best_valid_auprc": pretrain_manifest.get("best_valid_auprc"),
        "tahoe_cache": pretrain_manifest.get("tahoe_cache"),
    }

    # 1b) optional warm-cache epoch + concurrency probe
    if args.warm_rerun:
        warm_dir = out_root / "pretrain_warm_1epoch"
        code = _run(
            [args.python, train, "--stage", "pretrain", "--out", str(warm_dir),
             "--seed", str(args.seed + 1), "--epochs", "1", "--min-epochs", "1",
             "--max-train-rows", str(args.pretrain_rows),
             "--max-valid-rows", str(args.valid_rows),
             "--tahoe-cache", str(cache), "--cache-mode", "require", *common],
            out_root / "logs" / "pretrain_warm.log", env,
        )
        if must(code, "warm-cache pretrain"):
            warm_history = _history(warm_dir).get("history", [])
            record["steps"]["pretrain_warm"] = {
                "run_dir": str(warm_dir),
                "peak_rss_mb": _manifest(warm_dir).get("peak_rss_mb"),
                "per_epoch_seconds": [entry.get("seconds") for entry in warm_history],
                "per_epoch_valid_seconds": [entry.get("valid_seconds") for entry in warm_history],
            }

    if args.concurrency_probe:
        probes = []
        procs = []
        for offset in (0, 1):
            probe_dir = out_root / f"concurrency_probe_{offset}"
            command = [
                args.python, train, "--stage", "pretrain", "--out", str(probe_dir),
                "--seed", str(args.seed + 10 + offset), "--epochs", "1",
                "--min-epochs", "1", "--max-train-rows", str(args.pretrain_rows),
                "--max-valid-rows", str(args.valid_rows),
                "--tahoe-cache", str(cache), "--cache-mode", "auto", *common,
            ]
            log_path = out_root / "logs" / f"concurrency_probe_{offset}.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handle = open(log_path, "w", encoding="utf-8")
            handle.write("$ " + " ".join(command) + "\n")
            handle.flush()
            procs.append((subprocess.Popen(command, env=env, stdout=handle, stderr=subprocess.STDOUT),
                          probe_dir, handle))
        probe_started = time.perf_counter()
        for proc, probe_dir, handle in procs:
            code = proc.wait()
            handle.close()
            probes.append({
                "run_dir": str(probe_dir),
                "exit_code": int(code),
                "peak_rss_mb": _manifest(probe_dir).get("peak_rss_mb"),
            })
        record["steps"]["concurrency_probe"] = {
            "concurrent": 2,
            "wall_s": time.perf_counter() - probe_started,
            "runs": probes,
        }

    # 2/3) CM finetune fold, random and warm-started
    init_specs = {
        f"cm_fold{args.fold}_random": ["--init", "random"],
        f"cm_fold{args.fold}_tahoe": [
            "--init", "pretrained",
            "--init-checkpoint", str(pretrain_dir / "checkpoint" / "model.pt"),
        ],
    }
    for name, init_args in init_specs.items():
        run_dir = out_root / name
        code = _run(
            [args.python, train, "--stage", "finetune", "--fold", str(args.fold),
             "--out", str(run_dir), "--seed", str(args.seed), "--smoke",
             "--epochs", str(args.epochs), "--lr", str(args.cm_lr), *init_args, *common],
            out_root / "logs" / f"{name}.log", env,
        )
        if not must(code, f"finetune {name}"):
            return 1
        manifest = _manifest(run_dir)
        history = _history(run_dir).get("history", [])
        record["steps"][name] = {
            "tag": name,
            "init": manifest.get("init"),
            "checkpoint": manifest.get("checkpoint"),
            "peak_rss_mb": manifest.get("peak_rss_mb"),
            "wall_s": (manifest.get("timings") or {}).get("wall_s"),
            "per_epoch_seconds": [entry.get("seconds") for entry in history],
            "per_epoch_valid_seconds": [entry.get("valid_seconds") for entry in history],
            "per_epoch_valid_split": {
                "gather_s": [entry.get("valid_gather_seconds") for entry in history],
                "forward_s": [entry.get("valid_forward_seconds") for entry in history],
                "auprc_s": [entry.get("valid_metrics_seconds") for entry in history],
            },
            "best_valid_auprc": manifest.get("best_valid_auprc"),
            "rows": (manifest.get("data_counts") or {}),
        }

        # 4) predictions + evaluation for valid and test
        for role in ("valid", "test"):
            pred_dir = run_dir / "predictions" / role
            code = _run(
                [args.python, predict, "--stage", "finetune", "--fold", str(args.fold),
                 "--role", role, "--checkpoint", str(run_dir / "checkpoint" / "model.pt"),
                 "--out", str(pred_dir), *common],
                out_root / "logs" / f"{name}_predict_{role}.log", env,
            )
            if not must(code, f"predict {name} {role}"):
                return 1
            eval_dir = run_dir / "evaluation" / role
            code = _run(
                [args.python, evaluate, "--stage", "finetune", "--fold", str(args.fold),
                 "--role", role, "--pred-dir", str(pred_dir), "--out", str(eval_dir), *common],
                out_root / "logs" / f"{name}_eval_{role}.log", env,
            )
            if not must(code, f"evaluate {name} {role}"):
                return 1
            aggregate = _aggregate(eval_dir)
            record["steps"][name][f"eval_{role}"] = {
                "pred_dir": str(pred_dir),
                "eval_dir": str(eval_dir),
                "n_conditions": aggregate.get("n_conditions"),
                "n_rows": aggregate.get("n_rows"),
                "headline_k": aggregate.get("headline_k"),
                "headline_macro_auprc": (
                    (aggregate.get("by_k") or {}).get(str(aggregate.get("headline_k")), {})
                    .get("macro", {}).get("mean")
                ),
                "regression": aggregate.get("regression"),
            }
            if role == "test":
                code = _run(
                    [args.python, card, "--run-dir", str(run_dir), "--eval-dir", str(eval_dir)],
                    out_root / "logs" / f"{name}_card.log", env,
                )
                if not must(code, f"model card {name}"):
                    return 1

    # pretrain prediction/evaluation on the frozen valid split (plate13/14)
    pred_dir = pretrain_dir / "predictions" / "valid"
    code = _run(
        [args.python, predict, "--stage", "pretrain", "--role", "valid",
         "--checkpoint", str(pretrain_dir / "checkpoint" / "model.pt"),
         "--out", str(pred_dir), *common],
        out_root / "logs" / "pretrain_predict_valid.log", env,
    )
    if must(code, "pretrain predict valid"):
        eval_dir = pretrain_dir / "evaluation" / "valid"
        code = _run(
            [args.python, evaluate, "--stage", "pretrain", "--role", "valid",
             "--pred-dir", str(pred_dir), "--out", str(eval_dir), *common],
            out_root / "logs" / "pretrain_eval_valid.log", env,
        )
        if must(code, "pretrain evaluate valid"):
            aggregate = _aggregate(eval_dir)
            record["steps"]["pretrain"]["eval_valid"] = {
                "pred_dir": str(pred_dir),
                "eval_dir": str(eval_dir),
                "n_conditions": aggregate.get("n_conditions"),
                "n_rows": aggregate.get("n_rows"),
                "headline_k": aggregate.get("headline_k"),
                "headline_macro_auprc": (
                    (aggregate.get("by_k") or {}).get(str(aggregate.get("headline_k")), {})
                    .get("macro", {}).get("mean")
                ),
                "regression": aggregate.get("regression"),
            }
            _run(
                [args.python, card, "--run-dir", str(pretrain_dir), "--eval-dir", str(eval_dir)],
                out_root / "logs" / "pretrain_card.log", env,
            )

    record["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    record["total_wall_s"] = time.perf_counter() - started_all
    record["git"] = (pretrain_manifest.get("git") or {})
    record["data_artifacts"] = ((pretrain_manifest.get("data") or {}).get("artifacts") or {})

    # summary table (tracked record; force-added in Phase B)
    rows = [
        ["step", "run_dir", "peak_rss_mb", "wall_s", "per_epoch_s", "per_epoch_valid_s",
         "headline_auprc", "checkpoint_sha256"],
    ]

    def add_row(tag: str, step: dict, auprc: float | None = None) -> None:
        rows.append([
            tag, step.get("run_dir", ""), step.get("peak_rss_mb", ""), step.get("wall_s", ""),
            ";".join(str(x) for x in (step.get("per_epoch_seconds") or []) if x is not None),
            ";".join(str(x) for x in (step.get("per_epoch_valid_seconds") or []) if x is not None),
            "" if auprc is None else auprc,
            ((step.get("checkpoint") or {}).get("sha256") or "")[:16],
        ])

    pretrain_step = record["steps"].get("pretrain", {})
    add_row(
        "pretrain", pretrain_step,
        (pretrain_step.get("eval_valid") or {}).get("headline_macro_auprc"),
    )
    for name in init_specs:
        step = record["steps"].get(name, {})
        add_row(name, step, (step.get("eval_test") or {}).get("headline_macro_auprc"))

    import csv

    out_root.mkdir(parents=True, exist_ok=True)
    with open(out_root / "smoke_summary.csv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(rows)
    (out_root / "smoke_record.json").write_text(
        json.dumps(record, indent=2, sort_keys=False, default=str) + "\n", encoding="utf-8"
    )
    print("\nSMOKE SUMMARY (pipeline evidence only — not a science result)")
    for row in rows[1:]:
        print("  ", row[0], "| AUPRC", row[6], "| peak_rss_mb", row[2],
              "| epochs(s)", row[4], "| valid(s)", row[5])
    print(f"\nwrote {out_root / 'smoke_record.json'} and smoke_summary.csv "
          f"({record['total_wall_s']:.1f}s)")
    return 1 if record["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
