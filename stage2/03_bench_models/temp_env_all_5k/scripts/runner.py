#!/usr/bin/env python
"""Fold driver for the Stage 2 benchmark (ticket 05/10/12).

Responsibilities
-----------------
- Build the shared per-fold case bundles once (train/valid with expression,
  test metadata only) and verify the leakage gates.
- Invoke each adapter through its declared conda/venv/python environment; no
  benchmark module imports model internals.
- Run fit -> predict(valid) -> predict(test), enforce the shape/finiteness
  contract, average per-well predictions within `group_id`, and write the
  per-fold artifacts.
- Resume only units whose `status.json` is `completed` and whose artifacts are
  present; never fabricate a prediction for a failed/timed-out unit.
- Call the evaluator to write the leaderboard.

Baselines are registered in `configs/models.json`.  Real model entries are
added by tickets 13–15 without changing this driver.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import signal
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
ROOT = HERE.parents[3]
DATA = ARCHIVE / "data"
REGISTRY_PATH = ARCHIVE / "configs/models.json"

sys.path.insert(0, str(HERE))

from cases import BenchData, sha256_file  # noqa: E402
import evaluate  # noqa: E402


class RunnerError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Registry and environment bridge


def load_registry(path: Path = REGISTRY_PATH) -> dict:
    if not path.exists():
        raise RunnerError(f"missing model registry: {path}")
    registry = json.loads(path.read_text())
    if not isinstance(registry, dict):
        raise RunnerError("models.json must be an object keyed by model name")
    return registry


def python_command(env: dict) -> list[str]:
    kind = env.get("kind", "python")
    if kind == "python":
        return [sys.executable]
    if kind == "venv":
        path = Path(env["path"])
        if not path.is_absolute():
            path = ROOT / path
        return [str(path / "bin" / "python")]
    if kind == "conda":
        conda = os.environ.get("CONDA_EXE") or shutil.which("conda") or "conda"
        return [conda, "run", "--no-capture-output", "-n", str(env["name"]), "python"]
    raise RunnerError(f"unknown env kind: {kind}")


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _run_subprocess(
    command: list[str],
    log_path: Path,
    timeout_s: float,
    env_extra: dict[str, str],
) -> dict:
    """Run one adapter command, capture output, and kill the process group on timeout."""

    env = os.environ.copy()
    env.update(env_extra)
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("MKL_NUM_THREADS", "4")
    env.setdefault("OPENBLAS_NUM_THREADS", "4")
    env.setdefault("NUMEXPR_NUM_THREADS", "4")

    started = time.time()
    timed_out = False
    stdout = ""
    stderr = ""
    returncode: int | None = None
    try:
        process = subprocess.Popen(
            command,
            cwd=str(ROOT),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                stdout, stderr = process.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                except ProcessLookupError:
                    pass
                stdout, stderr = process.communicate()
        returncode = process.returncode
    except Exception as exc:  # noqa: BLE001 - surfaced in status.json
        stderr = f"{stderr}\nrunner exception: {exc!r}\n"
        returncode = -1
    wall = time.time() - started

    with open(log_path, "a", encoding="utf-8") as log:
        log.write(f"\n=== {_iso_now()} ===\n$ {' '.join(command)}\n")
        log.write(stdout or "")
        log.write(stderr or "")
        log.write(f"\nexit={returncode} timed_out={timed_out} wall_s={wall:.3f}\n")

    return {
        "returncode": returncode,
        "timed_out": timed_out,
        "stdout": stdout,
        "stderr": stderr,
        "wall_s": wall,
    }


# ---------------------------------------------------------------------------
# Case bundles


def ensure_cases(
    data: BenchData, fold: int, output_dir: Path, force: bool = False
) -> tuple[dict[str, Path], dict]:
    case_dir = output_dir / "_cases" / f"fold_{fold}"
    paths = data.write_cases(fold, case_dir, force=force)
    info = data.verify_cases(fold, paths)
    manifest_path = case_dir / "cases.json"
    if force or not manifest_path.exists():
        _write_json(manifest_path, data.cases_manifest(fold, paths))
    return paths, info


# ---------------------------------------------------------------------------
# Unit execution


REQUIRED_ARTIFACTS = (
    "config.json",
    "run.log",
    "checkpoint_path.txt",
    "pred_valid.npy",
    "pred_test.npy",
    "metrics.json",
    "valid_metrics.json",
    "well_ids.json",
    "group_ids.json",
    "status.json",
)


def unit_complete(run_dir: Path) -> bool:
    status_path = run_dir / "status.json"
    if not status_path.exists():
        return False
    try:
        status = json.loads(status_path.read_text())
    except Exception:
        return False
    if status.get("status") != "completed":
        return False
    if not all((run_dir / name).exists() for name in REQUIRED_ARTIFACTS):
        return False
    checkpoint = run_dir / "checkpoint_path.txt"
    try:
        checkpoint_path = Path(checkpoint.read_text().strip())
    except Exception:
        return False
    return checkpoint_path.exists()


def _load_case_n(path: Path) -> int:
    with np.load(path, allow_pickle=False) as npz:
        return int(len(npz["well_id"]))


def _validate_prediction(path: Path, n_expected: int, n_genes: int) -> np.ndarray:
    pred = np.load(path, allow_pickle=False)
    if pred.ndim != 2 or pred.shape[1] != n_genes:
        raise RunnerError(f"prediction shape {pred.shape}, expected (n, {n_genes})")
    if pred.shape[0] != n_expected:
        raise RunnerError(f"prediction rows {pred.shape[0]} != cases {n_expected}")
    pred = np.asarray(pred, dtype=np.float32)
    if not np.isfinite(pred).all():
        raise RunnerError("prediction contains NaN/Inf")
    np.save(path, pred)
    return pred


def _checkpoint_path(run_dir: Path) -> Path:
    marker = run_dir / "checkpoint_path.txt"
    if marker.exists():
        path = Path(marker.read_text().strip())
        if path.exists():
            return path
    default = run_dir / "checkpoint"
    if default.exists():
        return default
    raise RunnerError("adapter did not leave a checkpoint")


def _group_arrays(
    pred: np.ndarray, group_ids: np.ndarray, data: BenchData
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    unique, group_pred, counts = evaluate.group_mean_predictions(pred, group_ids)
    missing = [g for g in unique if g not in data.group_row]
    if missing:
        raise RunnerError(f"prediction references unknown groups: {missing[:3]}")
    true = np.stack([data.group_delta[data.group_row[g]] for g in unique]).astype(np.float32)
    return unique, group_pred.astype(np.float32), true, counts


def run_unit(
    data: BenchData,
    model: str,
    fold: int,
    entry: dict,
    case_paths: dict[str, Path],
    case_info: dict,
    output_dir: Path,
    seed: int,
    timeout_s: float,
    dry_run: bool,
    force: bool,
) -> dict:
    run_dir = output_dir / model / f"fold_{fold}"
    run_dir.mkdir(parents=True, exist_ok=True)
    status = {
        "model": model,
        "fold": fold,
        "status": "running",
        "started_at": _iso_now(),
        "seed": seed,
        "adapter": str(entry["adapter"]),
        "env": entry.get("env", {"kind": "python"}),
        "git_commit": _git_commit(),
        "n_train": case_info["n_train"],
        "n_valid": case_info["n_valid"],
        "n_test": case_info["n_test"],
        "n_train_groups": case_info["n_train_groups"],
        "n_valid_groups": case_info["n_valid_groups"],
        "n_test_groups": case_info["n_test_groups"],
        "dry_run": dry_run,
        "peak_gpu_mb": 0,
        "peak_rss_mb": 0.0,
        "fit": {},
        "predict_valid": {},
        "predict_test": {},
        "error": None,
    }

    if not force and unit_complete(run_dir):
        existing = json.loads((run_dir / "status.json").read_text())
        return existing

    adapter = (HERE / entry["adapter"]).resolve()
    if not adapter.exists():
        raise RunnerError(f"missing adapter: {adapter}")
    env_command = python_command(entry.get("env", {"kind": "python"}))
    extra = [str(x) for x in entry.get("extra_args", [])]
    controls = str((DATA / "controls.npz").resolve())
    config_path = run_dir / "config.json"
    log_path = run_dir / "run.log"
    if force and log_path.exists():
        log_path.unlink()

    config = {
        "model": model,
        "fold": fold,
        "seed": seed,
        "dry_run": dry_run,
        "adapter": str(adapter.relative_to(HERE)),
        "env": entry.get("env", {"kind": "python"}),
        "extra_args": extra,
        "baseline": bool(entry.get("baseline", False)),
        "baseline_params": entry.get("baseline_params", {}),
        "author_defaults": entry.get("author_defaults", {}),
        "data": {
            "controls": controls,
            "manifest_sha256": sha256_file(DATA / "manifest.json"),
        },
        "cases": {role: str(path.resolve()) for role, path in case_paths.items()},
        "budget": {"timeout_s": timeout_s},
    }
    _write_json(config_path, config)

    env_extra = {
        "OMP_NUM_THREADS": "4",
        "MKL_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "4",
        "NUMEXPR_NUM_THREADS": "4",
    }
    if dry_run:
        env_extra["CUDA_VISIBLE_DEVICES"] = ""
        env_extra["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

    def rss_mb() -> float:
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        return float(usage.ru_maxrss) / 1024.0

    try:
        fit_cmd = env_command + [
            str(adapter),
            "fit",
            "--train",
            str(case_paths["train"]),
            "--valid",
            str(case_paths["valid"]),
            "--config",
            str(config_path),
            "--run-dir",
            str(run_dir),
            "--controls",
            controls,
            *extra,
        ]
        fit = _run_subprocess(fit_cmd, log_path, timeout_s, env_extra)
        status["fit"] = {
            "returncode": fit["returncode"],
            "timed_out": fit["timed_out"],
            "wall_s": fit["wall_s"],
        }
        status["peak_rss_mb"] = max(status["peak_rss_mb"], rss_mb())
        if fit["timed_out"]:
            status["status"] = "timed_out"
            status["error"] = "fit timed out; predicting from the last checkpoint if present"
        elif fit["returncode"] != 0:
            status["status"] = "failed"
            status["error"] = (fit["stderr"] or "")[-2000:]
        else:
            status["status"] = "completed"

        if status["status"] == "failed":
            raise RunnerError(status["error"] or "fit failed")

        # A timed-out fit may still leave a usable checkpoint.  Crash failures
        # never fabricate a prediction.
        checkpoint = _checkpoint_path(run_dir)
        status["checkpoint"] = str(checkpoint)
        for role in ("valid", "test"):
            n_cases = _load_case_n(case_paths[role])
            pred_path = run_dir / f"pred_{role}.npy"
            predict_cmd = env_command + [
                str(adapter),
                "predict",
                "--cases",
                str(case_paths[role]),
                "--checkpoint",
                str(checkpoint),
                "--out",
                str(pred_path),
                "--controls",
                controls,
                *extra,
            ]
            result = _run_subprocess(predict_cmd, log_path, timeout_s, env_extra)
            status[f"predict_{role}"] = {
                "returncode": result["returncode"],
                "timed_out": result["timed_out"],
                "wall_s": result["wall_s"],
            }
            status["peak_rss_mb"] = max(status["peak_rss_mb"], rss_mb())
            if result["timed_out"]:
                raise RunnerError(f"predict {role} timed out")
            if result["returncode"] != 0:
                raise RunnerError(f"predict {role} failed: {(result['stderr'] or '')[-1000:]}")
            _validate_prediction(pred_path, n_cases, int(data.manifest["n_genes"]))
        # Keep `timed_out` if the fit exceeded the budget; otherwise completed.
        if status["status"] != "timed_out":
            status["status"] = "completed"
    except Exception as exc:  # noqa: BLE001 - recorded, never fabricated
        status["status"] = status["status"] if status["status"] in {"timed_out", "failed"} else "failed"
        status["error"] = f"{exc!r}\n{traceback.format_exc()[-2000:]}"
        status["finished_at"] = _iso_now()
        status["wall_clock_s"] = (
            datetime.fromisoformat(status["finished_at"]) - datetime.fromisoformat(status["started_at"])
        ).total_seconds()
        _write_json(run_dir / "status.json", status)
        return status

    # Group average, truth, and per-fold metrics.
    with np.load(case_paths["valid"], allow_pickle=False) as valid_npz:
        valid_wells = [str(x) for x in valid_npz["well_id"]]
        valid_groups = np.asarray(valid_npz["group_id"]).astype(str)
    with np.load(case_paths["test"], allow_pickle=False) as test_npz:
        test_wells = [str(x) for x in test_npz["well_id"]]
        test_groups = np.asarray(test_npz["group_id"]).astype(str)
    _write_json(run_dir / "well_ids.json", {"valid": valid_wells, "test": test_wells})
    _write_json(
        run_dir / "group_ids.json",
        {"valid": valid_groups.tolist(), "test": test_groups.tolist()},
    )

    pred_valid = np.load(run_dir / "pred_valid.npy", allow_pickle=False)
    pred_test = np.load(run_dir / "pred_test.npy", allow_pickle=False)
    valid_unique, valid_pred_groups, valid_true, _ = _group_arrays(pred_valid, valid_groups, data)
    test_unique, test_pred_groups, test_true, _ = _group_arrays(pred_test, test_groups, data)

    np.save(run_dir / "group_pred_valid.npy", valid_pred_groups)
    np.save(run_dir / "group_true_valid.npy", valid_true)
    np.save(run_dir / "group_pred_test.npy", test_pred_groups)
    np.save(run_dir / "group_true_test.npy", test_true)
    _write_json(
        run_dir / "group_order.json",
        {"valid": valid_unique.tolist(), "test": test_unique.tolist()},
    )

    valid_summary = evaluate.summarize_predictions(valid_true, pred_valid, valid_groups)
    test_summary = evaluate.summarize_predictions(test_true, pred_test, test_groups)
    valid_summary.update({"model": model, "fold": fold, "role": "valid"})
    test_summary.update({"model": model, "fold": fold, "role": "test"})
    _write_json(run_dir / "valid_metrics.json", valid_summary)
    _write_json(run_dir / "metrics.json", test_summary)

    status["finished_at"] = _iso_now()
    status["wall_clock_s"] = (
        datetime.fromisoformat(status["finished_at"]) - datetime.fromisoformat(status["started_at"])
    ).total_seconds()
    status["n_groups_valid"] = int(len(valid_unique))
    status["n_groups_test"] = int(len(test_unique))
    status["artifacts"] = {
        name: {"sha256": sha256_file(run_dir / name), "bytes": (run_dir / name).stat().st_size}
        for name in (
            "config.json",
            "run.log",
            "pred_valid.npy",
            "pred_test.npy",
            "metrics.json",
            "valid_metrics.json",
            "well_ids.json",
            "group_ids.json",
        )
        if (run_dir / name).exists()
    }
    _write_json(run_dir / "status.json", status)
    return status


# ---------------------------------------------------------------------------
# Driver


def _select_models(registry: dict, requested: list[str] | None, dry_run: bool) -> list[str]:
    if requested:
        unknown = sorted(set(requested) - set(registry))
        if unknown:
            raise RunnerError(f"unknown model(s): {unknown}")
        return requested
    if dry_run:
        return sorted(name for name, entry in registry.items() if entry.get("baseline"))
    return sorted(registry)


def _select_folds(requested: list[int] | None, dry_run: bool) -> list[int]:
    if requested:
        folds = sorted(set(requested))
    elif dry_run:
        folds = [0]
    else:
        folds = [0, 1, 2, 3, 4]
    bad = [f for f in folds if f not in range(5)]
    if bad:
        raise RunnerError(f"fold(s) outside 0..4: {bad}")
    return folds


def run(args: argparse.Namespace) -> dict:
    registry = load_registry()
    models = _select_models(registry, args.model, args.dry_run)
    folds = _select_folds(args.fold, args.dry_run)

    output_dir = Path(args.output_dir) if args.output_dir else (
        ARCHIVE / "output" / "dry_run" if args.dry_run else ARCHIVE / "output"
    )
    results_dir = Path(args.results_dir) if args.results_dir else (
        ARCHIVE / "results" / "dry_run" if args.dry_run else ARCHIVE / "results"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    timeout_s = float(args.timeout) if args.timeout else (600.0 if args.dry_run else 21_600.0)
    data = BenchData(DATA)
    case_paths: dict[int, dict[str, Path]] = {}
    case_infos: dict[int, dict] = {}
    for fold in folds:
        paths, info = ensure_cases(data, fold, output_dir, force=args.force)
        case_paths[fold] = paths
        case_infos[fold] = info

    statuses: list[dict] = []
    for model in models:
        entry = registry[model]
        for fold in folds:
            status = run_unit(
                data,
                model,
                fold,
                entry,
                case_paths[fold],
                case_infos[fold],
                output_dir,
                args.seed,
                timeout_s,
                args.dry_run,
                args.force,
            )
            statuses.append(status)
            print(
                f"[{status.get('status'):9s}] {model} fold {fold} "
                f"({status.get('wall_clock_s', 0):.1f}s)",
                flush=True,
            )
            if status.get("status") == "failed":
                print(f"  error: {str(status.get('error'))[:300]}", flush=True)
            elif status.get("status") == "timed_out":
                print("  timed out; prediction not fabricated", flush=True)

    status_columns = [
        "model", "fold", "status", "wall_clock_s", "n_train", "n_valid", "n_test",
        "fit_returncode", "predict_valid_returncode", "predict_test_returncode",
        "peak_rss_mb", "error",
    ]
    rows = []
    for s in statuses:
        rows.append(
            {
                "model": s.get("model"),
                "fold": s.get("fold"),
                "status": s.get("status"),
                "wall_clock_s": s.get("wall_clock_s"),
                "n_train": s.get("n_train"),
                "n_valid": s.get("n_valid"),
                "n_test": s.get("n_test"),
                "fit_returncode": s.get("fit", {}).get("returncode"),
                "predict_valid_returncode": s.get("predict_valid", {}).get("returncode"),
                "predict_test_returncode": s.get("predict_test", {}).get("returncode"),
                "peak_rss_mb": s.get("peak_rss_mb"),
                "error": (str(s.get("error"))[:500] if s.get("error") else ""),
            }
        )
    run_status = pd.DataFrame(rows, columns=status_columns)
    run_status = run_status.sort_values(["model", "fold"], kind="mergesort").reset_index(drop=True)
    run_status.to_csv(output_dir / "run_status.csv", index=False, lineterminator="\n")

    evaluation = evaluate.evaluate_run(
        output_dir,
        results_dir,
        models=models,
        folds=folds,
        registry=registry,
        data_dir=DATA,
    )
    return {
        "models": models,
        "folds": folds,
        "statuses": statuses,
        "output_dir": output_dir,
        "results_dir": results_dir,
        "evaluation": evaluation,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="CPU baselines, one fold")
    parser.add_argument("--model", action="append", default=None, help="model name (repeatable)")
    parser.add_argument("--fold", action="append", type=int, default=None, help="fold 0..4 (repeatable)")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--timeout", type=float, default=None, help="seconds per fit/predict call")
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--force", action="store_true", help="rerun completed units")
    args = parser.parse_args(argv)

    result = run(args)
    failed = [s for s in result["statuses"] if s.get("status") not in {"completed"}]
    print(
        f"\noutput: {result['output_dir']}\nresults: {result['results_dir']}\n"
        f"units: {len(result['statuses'])}, non-completed: {len(failed)}"
    )
    if len(result["evaluation"]["leaderboard"]):
        print(
            result["evaluation"]["leaderboard"][
                ["rank", "model", "headline_auprc", "n_folds", "winner_eligible"]
            ].to_string(index=False)
        )
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
