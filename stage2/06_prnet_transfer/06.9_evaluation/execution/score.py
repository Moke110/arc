#!/usr/bin/env python
"""Score the 30 selected checkpoints on the held-out split — one unit at a time.

For each unit the driver runs the ticket-#4 CLIs, never a copy of them:

    predict.py  --stage finetune --fold F --role {valid,test} --checkpoint <sha-verified>
    evaluate.py --stage finetune --fold F --role {valid,test} --pred-dir <pred>

``predict.py`` is target-free by construction (it cannot read treated expression), so the
prediction it writes is safely reproducible before any truth is read. ``evaluate.py`` is the
one evaluator that reads truth. The driver adds: resume semantics, per-step wall time and
peak RSS, whole-GPU peak sampling, verification of every produced manifest against the
verified checkpoint/unit identity, and the **test gate** (the ``test`` phase refuses to run
without ``--confirm-frozen-selection``, because the 30 units must be selected before the
first test row is scored).

Outputs per unit under ``output/units/<unit_id>/``::

    pred/    pred_expr.npy, delta.npy, meta.csv, manifest.json
    eval/    metrics.csv, per_condition.csv, aggregate.json, deltas.npz, manifest.json
    record.json

Run-level evidence: ``output/record/score_runs.csv`` and ``output/record/score_summary.json``.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Mapping

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402
from common import (  # noqa: E402
    DEFAULT_DATA_DIR,
    DEFAULT_GPU_PY,
    DEFAULT_OUT,
    DEFAULT_TABLE_CACHE,
    EvaluationError,
    load_json,
    mem_available_gib,
    now_utc,
    peak_rss_mb,
    sha256_file,
    write_json,
    write_text_atomic,
)
from units import fold_reference_for, load_unit_manifest  # noqa: E402

SCHEMA = "stage2/06.9_evaluation/score/v1"
PIPELINE_EXECUTION = Path(__file__).resolve().parents[2] / "06.4_pipeline" / "execution"
PREDICT_CLI = PIPELINE_EXECUTION / "predict.py"
EVALUATE_CLI = PIPELINE_EXECUTION / "evaluate.py"

RUN_COLUMNS = (
    "unit_id", "phase", "status", "fold", "seed", "init", "lr", "started_utc",
    "predict_exit", "predict_wall_s", "predict_peak_rss_mb",
    "evaluate_exit", "evaluate_wall_s", "evaluate_peak_rss_mb",
    "gpu_peak_mb", "eval_manifest_sha256", "error",
)


class ScoreError(EvaluationError):
    """A unit cannot be scored safely."""


# --------------------------------------------------------------------------- #
# Small process helpers
# --------------------------------------------------------------------------- #
def _run_child(command: list[str], cwd: Path) -> dict[str, Any]:
    """Run one child; return exit code, wall time and the child's peak RSS (KiB)."""
    env = dict(os.environ)
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    env.setdefault("OMP_NUM_THREADS", "2")
    env.setdefault("MKL_NUM_THREADS", "2")
    env.setdefault("OPENBLAS_NUM_THREADS", "2")
    started = time.perf_counter()
    started_utc = now_utc()
    process = subprocess.Popen(command, cwd=str(cwd), env=env)
    try:
        _pid, status, usage = os.wait4(process.pid, 0)
        exit_code = os.waitstatus_to_exitcode(status)
        process.returncode = exit_code
    except BaseException:  # pragma: no cover - interrupted run: never leave a zombie
        process.kill()
        process.wait()
        raise
    return {
        "command": list(command),
        "exit_code": int(exit_code),
        "wall_s": time.perf_counter() - started,
        "peak_rss_mb": float(usage.ru_maxrss) / 1024.0,
        "started_utc": started_utc,
        "finished_utc": now_utc(),
    }


class GpuSampler(threading.Thread):
    """Sample whole-GPU memory use while a unit runs (no per-process attribution expected)."""

    def __init__(self, interval_s: float = 2.0) -> None:
        super().__init__(daemon=True)
        self.interval_s = float(interval_s)
        self._stop = threading.Event()
        self.samples: list[int] = []
        self.available = True

    def run(self) -> None:  # pragma: no cover - hardware dependent
        while not self._stop.is_set():
            try:
                done = subprocess.run(
                    [
                        "nvidia-smi",
                        "--query-gpu=memory.used",
                        "--format=csv,noheader,nounits",
                    ],
                    capture_output=True,
                    text=True,
                    timeout=5,
                    check=False,
                )
                if done.returncode == 0 and done.stdout.strip():
                    values = [
                        int(line.strip())
                        for line in done.stdout.splitlines()
                        if line.strip().isdigit()
                    ]
                    if values:
                        self.samples.append(max(values))
                else:
                    self.available = False
                    return
            except (OSError, subprocess.SubprocessError, ValueError):
                self.available = False
                return
            self._stop.wait(self.interval_s)

    def stop(self) -> dict[str, Any]:
        self._stop.set()
        self.join(timeout=5)
        return {
            "sampled": bool(self.samples),
            "peak_mb": max(self.samples) if self.samples else None,
            "baseline_mb": self.samples[0] if self.samples else None,
            "n_samples": len(self.samples),
            "interval_s": self.interval_s,
        }


# --------------------------------------------------------------------------- #
# Pre-flight and reference checks
# --------------------------------------------------------------------------- #
def _preflight(min_available_gib: float) -> dict[str, Any]:
    available = mem_available_gib()
    record: dict[str, Any] = {
        "mem_available_gib": available,
        "min_available_gib": float(min_available_gib),
        "gpu_mem_used_mb": None,
    }
    try:
        done = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10, check=False,
        )
        if done.returncode == 0 and done.stdout.strip():
            record["gpu_mem_used_mb"] = max(
                int(line.strip()) for line in done.stdout.splitlines() if line.strip()
            )
    except (OSError, subprocess.SubprocessError, ValueError):
        pass
    if available and available < float(min_available_gib):
        raise ScoreError(
            f"only {available:.1f} GiB RAM available (need >= {min_available_gib:.1f}); "
            "abort before loading the frozen input"
        )
    return record


def _verify_unit_outputs(
    unit: Mapping[str, Any],
    phase: str,
    pred_dir: Path,
    eval_dir: Path,
    reference: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Cross-check the produced prediction/evaluation manifests against the unit identity."""
    pred_manifest = load_json(pred_dir / "manifest.json")
    eval_manifest = load_json(eval_dir / "manifest.json")
    aggregate = load_json(eval_dir / "aggregate.json")

    # evaluate.py records the scored checkpoint under predictions.checkpoint_sha256
    evaluated_checkpoint = (eval_manifest.get("predictions") or {}).get("checkpoint_sha256")
    if evaluated_checkpoint is None:
        evaluated_checkpoint = (eval_manifest.get("checkpoint") or {}).get("sha256")
    if str(evaluated_checkpoint) != str(unit["checkpoint"]["sha256"]):
        raise ScoreError(
            f"{unit['unit_id']}: evaluation manifest checkpoint sha256 "
            f"{evaluated_checkpoint} != unit {unit['checkpoint']['sha256']}"
        )
    pred_checkpoint = (pred_manifest.get("checkpoint") or {}).get("sha256")
    if str(pred_checkpoint) != str(unit["checkpoint"]["sha256"]):
        raise ScoreError(f"{unit['unit_id']}: prediction manifest checkpoint disagrees")
    if str(eval_manifest.get("stage")) != "finetune":
        raise ScoreError(f"{unit['unit_id']}: evaluation stage is not finetune")
    if str(eval_manifest.get("role")) != phase:
        raise ScoreError(f"{unit['unit_id']}: evaluation role {eval_manifest.get('role')} != {phase}")
    if int(eval_manifest.get("fold", -1)) != int(unit["fold"]):
        raise ScoreError(f"{unit['unit_id']}: evaluation fold disagrees")
    # predict.py records the raw-array digest; evaluate.py hashes the .npy file, so the
    # cross-check is the file digest the evaluator actually scored.
    predicted_sha = ((pred_manifest.get("predictions") or {}).get("pred_expr_sha256"))
    evaluated_sha = ((eval_manifest.get("predictions") or {}).get("pred_expr_sha256"))
    pred_expr_file_sha = sha256_file(pred_dir / "pred_expr.npy")
    if pred_expr_file_sha != evaluated_sha:
        raise ScoreError(
            f"{unit['unit_id']}: evaluator scored a different prediction array "
            f"(file {pred_expr_file_sha[:12]}…, manifest {evaluated_sha})"
        )
    if not bool(pred_manifest.get("repeated_inference_identical")):
        raise ScoreError(f"{unit['unit_id']}: repeated inference was not byte-identical")

    n_rows = int(aggregate.get("n_rows", -1))
    n_conditions = int(aggregate.get("n_conditions", -1))
    headline_k = int(aggregate.get("headline_k", -1))
    if headline_k != 200:
        raise ScoreError(f"{unit['unit_id']}: headline k is {headline_k}, expected 200")
    macro = (((aggregate.get("by_k") or {}).get("200")) or {}).get("macro") or {}
    if not isinstance(macro.get("mean"), (int, float)) or macro["mean"] != macro["mean"]:
        raise ScoreError(f"{unit['unit_id']}: no finite headline AUPRC@200")
    checks: dict[str, Any] = {
        "n_rows": n_rows,
        "n_conditions": n_conditions,
        "headline_k": headline_k,
        "checkpoint_sha256": str(unit["checkpoint"]["sha256"]),
        "pred_expr_sha256": predicted_sha,
        "pred_expr_file_sha256": pred_expr_file_sha,
        "eval_manifest_sha256": sha256_file(eval_dir / "manifest.json"),
        "pred_manifest_sha256": sha256_file(pred_dir / "manifest.json"),
        "repeated_inference_identical": bool(pred_manifest.get("repeated_inference_identical")),
    }
    if reference:
        if phase == "test":
            if n_rows != int(reference["n_test_rows"]):
                raise ScoreError(
                    f"{unit['unit_id']}: scored {n_rows} test rows, #8's identical population "
                    f"has {reference['n_test_rows']}"
                )
            if n_conditions != int(reference["n_test_conditions"]):
                raise ScoreError(
                    f"{unit['unit_id']}: scored {n_conditions} test conditions, #8 records "
                    f"{reference['n_test_conditions']}"
                )
        elif phase == "valid":
            if n_rows != int(reference["n_valid_rows"]):
                raise ScoreError(
                    f"{unit['unit_id']}: scored {n_rows} valid rows, #8 records "
                    f"{reference['n_valid_rows']}"
                )
        checks["fold_reference"] = dict(reference)
    return checks


# --------------------------------------------------------------------------- #
# One unit
# --------------------------------------------------------------------------- #
def score_unit(
    unit: Mapping[str, Any],
    *,
    phase: str,
    units_payload: Mapping[str, Any],
    out: Path,
    data_dir: Path,
    table_cache: Path | None,
    hash_cache: Path | None,
    gpu_py: str,
    force: bool,
) -> dict[str, Any]:
    unit_id = str(unit["unit_id"])
    unit_dir = out / "units" / unit_id
    pred_dir = unit_dir / "pred"
    eval_dir = unit_dir / "eval"
    record_path = unit_dir / "record.json"
    reference = fold_reference_for(units_payload, int(unit["fold"]))
    base = {
        "schema": SCHEMA,
        "kind": "evaluation_unit",
        "unit_id": unit_id,
        "phase": phase,
        "fold": int(unit["fold"]),
        "seed": int(unit["seed"]),
        "init": unit["init_name"],
        "lr": float(unit["lr"]),
        "checkpoint_sha256": str(unit["checkpoint"]["sha256"]),
        "checkpoint_path": str(unit["checkpoint"]["path"]),
        "config_sha256": str(unit["config_sha256"]),
        "data_identity": str(unit["data_identity"]),
        "created_utc": now_utc(),
    }

    if record_path.exists() and not force:
        previous = load_json(record_path)
        if (
            previous.get("status") == "completed"
            and previous.get("phase") == phase
            and str(previous.get("checkpoint_sha256")) == str(unit["checkpoint"]["sha256"])
        ):
            stored_eval = (previous.get("verification") or {}).get("eval_manifest_sha256")
            if (
                (eval_dir / "manifest.json").exists()
                and (pred_dir / "manifest.json").exists()
                and stored_eval == sha256_file(eval_dir / "manifest.json")
            ):
                previous["status"] = "skipped_complete"
                previous["skipped_utc"] = now_utc()
                return previous

    if force:
        for directory in (pred_dir, eval_dir):
            if directory.exists():
                for child in sorted(directory.rglob("*"), reverse=True):
                    if child.is_file():
                        child.unlink()
                    elif child.is_dir():
                        child.rmdir()
                directory.rmdir()
    pred_dir.mkdir(parents=True, exist_ok=True)
    eval_dir.mkdir(parents=True, exist_ok=True)

    common_data = ["--data-dir", str(data_dir)]
    if hash_cache is not None:
        common_data += ["--hash-cache", str(hash_cache)]
    if table_cache is not None:
        common_data += ["--table-cache", str(table_cache)]
    predict_command = [
        gpu_py, str(PREDICT_CLI), "--stage", "finetune",
        "--fold", str(unit["fold"]), "--role", phase,
        "--checkpoint", str(unit["checkpoint"]["path"]),
        "--out", str(pred_dir), *common_data,
        "--device", "cuda", "--repeats", "2",
    ]
    evaluate_command = [
        gpu_py, str(EVALUATE_CLI), "--stage", "finetune",
        "--fold", str(unit["fold"]), "--role", phase,
        "--pred-dir", str(pred_dir), "--out", str(eval_dir), *common_data,
    ]

    record = dict(base)
    record["status"] = "running"
    record["steps"] = {}
    sampler = GpuSampler(interval_s=2.0)
    sampler.start()
    try:
        record["steps"]["predict"] = _run_child(predict_command, unit_dir)
        if record["steps"]["predict"]["exit_code"] != 0:
            raise ScoreError(
                f"{unit_id}: predict.py exited {record['steps']['predict']['exit_code']}"
            )
        record["steps"]["evaluate"] = _run_child(evaluate_command, unit_dir)
        if record["steps"]["evaluate"]["exit_code"] != 0:
            raise ScoreError(
                f"{unit_id}: evaluate.py exited {record['steps']['evaluate']['exit_code']}"
            )
        verification = _verify_unit_outputs(
            unit, phase, pred_dir, eval_dir, reference
        )
        record["verification"] = verification
        record["status"] = "completed"
        record["finished_utc"] = now_utc()
        record["wall_s"] = sum(
            step["wall_s"] for step in record["steps"].values()
        )
        record["host_peak_rss_mb"] = max(
            step["peak_rss_mb"] for step in record["steps"].values()
        )
    except EvaluationError as exc:
        record["status"] = "failed"
        record["error"] = str(exc)
        record["finished_utc"] = now_utc()
        record["wall_s"] = sum(
            step.get("wall_s", 0.0) for step in record["steps"].values()
        )
        record["host_peak_rss_mb"] = max(
            (step.get("peak_rss_mb", 0.0) for step in record["steps"].values()),
            default=0.0,
        )
    finally:
        record["gpu"] = sampler.stop()
        write_json(record_path, record)
    return record


# --------------------------------------------------------------------------- #
# Queue + record
# --------------------------------------------------------------------------- #
def _run_row(record: Mapping[str, Any]) -> dict[str, str]:
    """Flatten one unit record to the ``score_runs.csv`` row."""
    steps = record.get("steps") or {}
    predict = steps.get("predict") or {}
    evaluate = steps.get("evaluate") or {}
    return {
        "unit_id": str(record.get("unit_id", "")),
        "phase": str(record.get("phase", "")),
        "status": str(record.get("status", "")),
        "fold": str(record.get("fold", "")),
        "seed": str(record.get("seed", "")),
        "init": str(record.get("init", "")),
        "lr": str(record.get("lr", "")),
        "started_utc": str(
            predict.get("started_utc") or record.get("created_utc") or ""
        ),
        "predict_exit": str(predict.get("exit_code", "")),
        "predict_wall_s": f"{predict.get('wall_s', ''):.3f}" if predict else "",
        "predict_peak_rss_mb": f"{predict.get('peak_rss_mb', ''):.1f}" if predict else "",
        "evaluate_exit": str(evaluate.get("exit_code", "")),
        "evaluate_wall_s": f"{evaluate.get('wall_s', ''):.3f}" if evaluate else "",
        "evaluate_peak_rss_mb": (
            f"{evaluate.get('peak_rss_mb', ''):.1f}" if evaluate else ""
        ),
        "gpu_peak_mb": str((record.get("gpu") or {}).get("peak_mb", "")),
        "eval_manifest_sha256": str(
            (record.get("verification") or {}).get("eval_manifest_sha256", "")
        ),
        "error": str(record.get("error") or "").replace(",", ";"),
    }


def _read_runs_rows(out: Path) -> list[dict[str, str]]:
    path = out / "record" / "score_runs.csv"
    if not path.exists():
        return []
    try:
        frame = common.read_table(path)
    except (OSError, ValueError):  # pragma: no cover - corrupt scratch only
        return []
    return [dict(row) for row in frame.to_dict("records")]


def _row_started_utc(out: Path, unit_id: str) -> str:
    path = out / "units" / str(unit_id) / "record.json"
    if not path.exists():
        return ""
    try:
        record = load_json(path)
    except (OSError, ValueError):  # pragma: no cover
        return ""
    steps = record.get("steps") or {}
    return str(
        (steps.get("predict") or {}).get("started_utc")
        or record.get("created_utc")
        or ""
    )


def _write_runs_csv(out: Path, records: list[Mapping[str, Any]]) -> list[dict[str, str]]:
    """Write the run table cumulatively: a subset re-run must not drop earlier units."""
    rows_by_unit: dict[str, dict[str, str]] = {}
    for row in _read_runs_rows(out):
        unit_id = str(row.get("unit_id", ""))
        if not unit_id:
            continue
        if not row.get("started_utc"):
            row["started_utc"] = _row_started_utc(out, unit_id)
        rows_by_unit[unit_id] = row
    for record in records:
        row = _run_row(record)
        rows_by_unit[row["unit_id"]] = row
    rows = list(rows_by_unit.values())
    header = ",".join(RUN_COLUMNS)
    lines = [header]
    for row in rows:
        lines.append(",".join(str(row.get(column, "")) for column in RUN_COLUMNS))
    write_text_atomic(out / "record" / "score_runs.csv", "\n".join(lines) + "\n")
    return rows


def _summary_from_rows(
    rows: list[Mapping[str, Any]],
    *,
    phase: str,
    preflight: Mapping[str, Any],
    units_manifest: Path,
    units_manifest_sha256: str,
    selection_complete: bool,
) -> dict[str, Any]:
    """Cumulative run summary over every unit in ``score_runs.csv`` (all invocations)."""

    def number(value: Any) -> float:
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0.0

    completed = [r for r in rows if r.get("status") == "completed"]
    skipped = [r for r in rows if r.get("status") == "skipped_complete"]
    failed = [r for r in rows if r.get("status") == "failed"]
    wall = sum(
        number(r.get("predict_wall_s")) + number(r.get("evaluate_wall_s")) for r in rows
    )
    host_peak = max(
        (
            max(number(r.get("predict_peak_rss_mb")), number(r.get("evaluate_peak_rss_mb")))
            for r in rows
        ),
        default=0.0,
    )
    gpu_peak = max((number(r.get("gpu_peak_mb")) for r in rows), default=0.0)
    started = min(
        (str(r.get("started_utc")) for r in rows if r.get("started_utc")), default=""
    )
    return {
        "schema": SCHEMA + "/summary",
        "kind": "score_summary",
        "phase": phase,
        "created_utc": now_utc(),
        "started_utc": started,
        "units_planned": len(rows),
        "units_completed": len(completed),
        "units_skipped": len(skipped),
        "units_failed": len(failed),
        "failures": [
            {
                "unit_id": r.get("unit_id"),
                "status": r.get("status"),
                "error": r.get("error"),
            }
            for r in failed
        ],
        "wall_s": wall,
        "host_peak_rss_mb": host_peak,
        "gpu_peak_mb": int(gpu_peak) if gpu_peak else 0,
        "preflight": dict(preflight),
        "units_manifest": str(units_manifest),
        "units_manifest_sha256": units_manifest_sha256,
        "selection_complete": bool(selection_complete),
        "test_role_scored": phase == "test",
        "peak_rss_mb_process": peak_rss_mb(),
        "note": (
            "cumulative over every unit in score_runs.csv (a subset --only-unit re-run "
            "does not drop earlier units); wall_s is the sum of per-unit predict+evaluate "
            "wall time"
        ),
    }


def run(
    *,
    units_path: Path,
    phase: str,
    out: Path,
    data_dir: Path,
    table_cache: Path | None,
    hash_cache: Path | None,
    gpu_py: str,
    force: bool,
    only_units: list[str] | None,
    limit: int | None,
    plan: bool,
    min_available_gib: float,
) -> int:
    out = Path(out).resolve()
    data_dir = Path(data_dir).resolve()
    units_path = Path(units_path).resolve()
    try:
        units_payload = load_unit_manifest(units_path, require_complete=(phase == "test"))
    except EvaluationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    units = list(units_payload["units"])
    if only_units:
        wanted = set(only_units)
        units = [unit for unit in units if unit["unit_id"] in wanted]
        missing = wanted - {unit["unit_id"] for unit in units}
        if missing:
            print(f"error: unknown unit ids: {sorted(missing)}", file=sys.stderr)
            return 2
    if limit is not None:
        units = units[: int(limit)]
    if not units:
        print("error: no units selected", file=sys.stderr)
        return 2
    print(
        f"{phase}-phase scoring: {len(units)} unit(s), folds "
        f"{sorted({int(u['fold']) for u in units})}, inits "
        f"{sorted({u['init_name'] for u in units})}"
    )

    out.mkdir(parents=True, exist_ok=True)
    runs: list[dict[str, Any]] = []
    started = time.perf_counter()
    started_utc = now_utc()
    preflight: dict[str, Any] = {}
    try:
        preflight = _preflight(min_available_gib)
    except EvaluationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 3
    print(f"pre-flight: {json.dumps(preflight)}")
    if plan:
        for unit in units:
            reference = fold_reference_for(units_payload, int(unit["fold"])) or {}
            expected = (
                reference.get("n_test_rows")
                if phase == "test"
                else reference.get("n_valid_rows")
            )
            print(
                f"  {unit['unit_id']}: fold {unit['fold']} seed {unit['seed']} "
                f"init {unit['init_name']} lr {unit['lr']} rows~{expected}"
            )
        return 0

    for index, unit in enumerate(units, start=1):
        print(f"[{index}/{len(units)}] {unit['unit_id']} ...", flush=True)
        record = score_unit(
            unit,
            phase=phase,
            units_payload=units_payload,
            out=out,
            data_dir=data_dir,
            table_cache=table_cache,
            hash_cache=hash_cache,
            gpu_py=gpu_py,
            force=force,
        )
        runs.append(record)
        _write_runs_csv(out, runs)
        print(
            f"    {record.get('status')} in "
            f"{float(record.get('wall_s') or 0.0):.1f} s "
            f"(host peak {(record.get('host_peak_rss_mb') or 0.0):.0f} MiB, "
            f"GPU peak {(record.get('gpu') or {}).get('peak_mb')} MiB)",
            flush=True,
        )

    merged_rows = _write_runs_csv(out, runs)
    summary = _summary_from_rows(
        merged_rows,
        phase=phase,
        preflight=preflight,
        units_manifest=units_path,
        units_manifest_sha256=sha256_file(units_path),
        selection_complete=bool(units_payload.get("complete")),
    )
    write_json(out / "record" / "score_summary.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "failures"}, indent=2))
    if summary["units_failed"]:
        print(json.dumps({"failures": summary["failures"]}, indent=2), file=sys.stderr)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Score the selected PRnet transfer units (ticket #9).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--units", type=Path, default=DEFAULT_OUT / "record" / "units.json"
    )
    parser.add_argument("--phase", choices=("valid", "test"), required=True)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--table-cache", type=Path, default=DEFAULT_TABLE_CACHE)
    parser.add_argument(
        "--hash-cache", type=Path, default=DEFAULT_OUT / "hash_cache.json",
        help="frozen-artifact hash memo (first run hashes 16.4 GB once)",
    )
    parser.add_argument("--gpu-py", default=DEFAULT_GPU_PY)
    parser.add_argument("--force", action="store_true", help="re-score even completed units")
    parser.add_argument(
        "--only-unit", action="append", default=None, metavar="UNIT_ID",
        help="score a subset (harness smoke); repeatable",
    )
    parser.add_argument("--limit", type=int, default=None, help="score at most N units")
    parser.add_argument("--plan", action="store_true", help="print the plan, score nothing")
    parser.add_argument(
        "--min-available-gib", type=float, default=6.0,
        help="abort below this MemAvailable",
    )
    parser.add_argument(
        "--confirm-frozen-selection", action="store_true",
        help="required for --phase test: asserts #6/#7 selections are final and frozen",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    # Resolve every path once: a relative --out would otherwise be re-interpreted by a
    # child whose cwd is the (relative) unit directory, nesting the outputs.
    ns.units = Path(ns.units).resolve()
    ns.out = Path(ns.out).resolve()
    ns.data_dir = Path(ns.data_dir).resolve()
    if ns.table_cache is not None:
        ns.table_cache = Path(ns.table_cache).resolve()
    if ns.hash_cache is not None:
        ns.hash_cache = Path(ns.hash_cache).resolve()
    if ns.phase == "test" and not ns.confirm_frozen_selection:
        print(
            "refusing --phase test: the 30 units must be selected and frozen before the "
            "first test row is scored. Re-run with --confirm-frozen-selection only after "
            "the Lead confirms both #6 and #7 matrices are final and #5's sweep is done.",
            file=sys.stderr,
        )
        return 4
    if ns.phase == "test":
        print(
            "=== TEST PHASE ===\n"
            "Frozen selection confirmed by the caller. This is the single touch of the "
            "test split for the first-round study; no result here may feed back into LR, "
            "epoch, seed or capacity choice.",
            flush=True,
        )
    return run(
        units_path=ns.units,
        phase=ns.phase,
        out=ns.out,
        data_dir=ns.data_dir,
        table_cache=ns.table_cache if ns.table_cache and Path(ns.table_cache).exists() else None,
        hash_cache=ns.hash_cache,
        gpu_py=ns.gpu_py,
        force=ns.force,
        only_units=ns.only_unit,
        limit=ns.limit,
        plan=ns.plan,
        min_available_gib=ns.min_available_gib,
    )


if __name__ == "__main__":
    raise SystemExit(main())
