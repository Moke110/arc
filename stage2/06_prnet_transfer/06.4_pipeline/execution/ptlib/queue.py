"""Resume-safe parallel run queue for training units (#5/#6/#7 driver).

A queue is a JSON file of units::

    {"units": [
       {"unit_id": "tahoe_seed20261021",
        "stage": "pretrain",
        "args": ["--seed", "20261021", "--lr", "0.001"]},
       {"unit_id": "fold0_seed20261021_lr3e-4_tahoe",
        "stage": "finetune",
        "args": ["--fold", "0", "--seed", "20261021", "--lr", "3e-4",
                 "--init", "pretrained", "--init-checkpoint", "/path/model.pt"]}
    ]}

Each unit runs ``train.py`` as a subprocess with its own run directory
``<out-root>/<unit_id>``.  The queue is **resume-safe**: a unit whose run
directory holds a completed manifest with a matching config fingerprint (and a
checkpoint whose sha256 matches) is skipped; a unit whose manifest is missing,
failed, tampered or built with different arguments is re-run.  A status CSV is
rewritten after every completion with exit codes and timings.

Parallelism: ``--workers N`` (1–4; the brief caps CM candidate training at 4
concurrent processes with ``OMP_NUM_THREADS=2`` each).  The pretrain memory rule
in the README bounds concurrent *pretrain* seeds to 2 on this host.
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from . import config as config_mod
from . import manifest as manifest_mod
from . import paths

__all__ = [
    "QueueError",
    "QueueUnit",
    "QueueOptions",
    "load_queue",
    "run_queue",
    "unit_run_dir",
    "expected_config_sha256",
]

log = logging.getLogger("ptlib.queue")

MAX_WORKERS = 4
STATUS_COLUMNS = (
    "unit_id", "stage", "status", "exit_code", "started_utc", "finished_utc",
    "wall_s", "run_dir", "error",
)


class QueueError(RuntimeError):
    """The queue file or a unit definition is invalid."""


@dataclass
class QueueUnit:
    unit_id: str
    stage: str
    args: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.unit_id or any(ch in self.unit_id for ch in ("/", "\\", "..")):
            raise QueueError(f"unsafe unit_id {self.unit_id!r}")
        if self.stage not in ("pretrain", "finetune"):
            raise QueueError(f"unit {self.unit_id}: stage must be pretrain or finetune")
        if self.stage == "finetune" and "--fold" not in self.args:
            raise QueueError(f"unit {self.unit_id}: finetune units need --fold")
        if "--out" in self.args or "--data-dir" in self.args:
            raise QueueError(
                f"unit {self.unit_id}: do not pass --out/--data-dir; the queue owns them"
            )
        if "--stage" in self.args:
            raise QueueError(f"unit {self.unit_id}: do not pass --stage; the unit declares it")


@dataclass
class QueueOptions:
    queue: Path
    out_root: Path
    data_dir: Path = paths.DEFAULT_FROZEN_INPUT
    status_csv: Path | None = None
    workers: int = 1
    force: bool = False
    dry_run: bool = False
    omp_threads: int = 2
    table_cache: Path | None = None
    python: str = sys.executable
    train_script: Path = field(default_factory=lambda: paths.EXECUTION / "train.py")


def load_queue(path: Path | str) -> list[QueueUnit]:
    path = Path(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise QueueError(f"queue file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise QueueError(f"queue file is not valid JSON: {path} ({exc})") from exc
    units_raw = payload.get("units") if isinstance(payload, Mapping) else payload
    if not isinstance(units_raw, list) or not units_raw:
        raise QueueError(f"queue {path} must hold a non-empty 'units' list")
    units: list[QueueUnit] = []
    seen: set[str] = set()
    for record in units_raw:
        if not isinstance(record, Mapping):
            raise QueueError("each queue unit must be an object")
        unit_id = str(record.get("unit_id") or record.get("id") or "")
        stage = str(record.get("stage") or "")
        raw_args = record.get("args", record.get("argv", []))
        if isinstance(raw_args, Mapping):
            args: list[str] = []
            for key, value in raw_args.items():
                flag = key if key.startswith("--") else "--" + str(key).replace("_", "-")
                if value is True:
                    args.append(flag)
                elif value in (False, None):
                    continue
                else:
                    args.extend([flag, str(value)])
        elif isinstance(raw_args, list):
            args = [str(x) for x in raw_args]
        else:
            raise QueueError(f"unit {unit_id}: args must be a list or an object")
        unit = QueueUnit(unit_id=unit_id, stage=stage, args=args)
        if unit.unit_id in seen:
            raise QueueError(f"duplicate unit_id {unit.unit_id!r}")
        seen.add(unit.unit_id)
        units.append(unit)
    return units


def unit_run_dir(out_root: Path | str, unit: QueueUnit) -> Path:
    return Path(out_root) / unit.unit_id


def _locked_config(data_dir: Path) -> dict[str, Any]:
    cfg, _notes = config_mod.read_locked_config(Path(data_dir) / "locked_config.json")
    return cfg


def expected_config_sha256(unit: QueueUnit, data_dir: Path | str) -> str:
    """Resolve the unit's run config and return its fingerprint (no hashing)."""
    from . import train as train_mod

    locked = _locked_config(Path(data_dir))
    out_dir = Path(".")  # placeholder; not part of the fingerprint
    argv = [
        "--stage", unit.stage,
        "--out", str(out_dir),
        *unit.args,
        "--data-dir", str(data_dir),
    ]
    options = train_mod.options_from_args(argv)
    cfg, _notes = train_mod.resolve_run_config(locked, options)
    return config_mod.config_fingerprint(cfg)


def _unit_command(
    unit: QueueUnit,
    data_dir: Path,
    run_dir: Path,
    python: str,
    script: Path,
    table_cache: Path | None = None,
) -> list[str]:
    command = [
        python, str(script),
        "--stage", unit.stage,
        "--out", str(run_dir),
        "--data-dir", str(data_dir),
        *unit.args,
    ]
    if table_cache is not None:
        command.extend(["--table-cache", str(table_cache)])
    return command


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _write_status(path: Path, rows: list[dict[str, Any]]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in STATUS_COLUMNS})
    tmp.replace(path)


def _complete(unit: QueueUnit, run_dir: Path, expected: str) -> tuple[bool, str]:
    return manifest_mod.unit_complete(
        run_dir, expected_config_sha256=expected, require_checkpoint=True
    )


def run_queue(options: QueueOptions) -> dict[str, Any]:
    """Execute every not-yet-complete unit; return a summary dict."""
    options.out_root = Path(options.out_root)
    if not 1 <= int(options.workers) <= MAX_WORKERS:
        raise SystemExit(
            f"--workers must be in 1..{MAX_WORKERS} (the brief caps CM candidate training "
            "at 4 concurrent processes)"
        )
    units = load_queue(options.queue)
    data_dir = Path(options.data_dir)
    if not (data_dir / "locked_config.json").exists():
        raise QueueError(f"--data-dir has no locked_config.json: {data_dir}")
    expected = {unit.unit_id: expected_config_sha256(unit, data_dir) for unit in units}
    run_dirs = {unit.unit_id: unit_run_dir(options.out_root, unit) for unit in units}
    status_path = Path(options.status_csv) if options.status_csv else options.out_root / "queue_status.csv"
    logs_dir = options.out_root / "queue_logs"

    status: dict[str, dict[str, Any]] = {}
    for unit in units:
        complete, reason = (False, "not started")
        if not options.force:
            complete, reason = _complete(unit, run_dirs[unit.unit_id], expected[unit.unit_id])
            if complete:
                status[unit.unit_id] = {
                    "unit_id": unit.unit_id, "stage": unit.stage,
                    "status": "skipped_complete", "exit_code": 0,
                    "started_utc": "", "finished_utc": "", "wall_s": 0.0,
                    "run_dir": str(run_dirs[unit.unit_id]), "error": reason,
                }
        if options.dry_run:
            log.info(
                "PLAN %-40s %-8s %s (%s)",
                unit.unit_id, unit.stage,
                "would skip (complete)" if complete else "would run",
                reason,
            )
    if options.dry_run:
        planned = [u for u in units if u.unit_id not in status]
        log.info(
            "dry run: %d unit(s); %d complete, %d to run, workers=%d",
            len(units), len(status), len(planned), options.workers,
        )
        for unit in planned:
            log.info("  %s: %s", unit.unit_id, _unit_command(
                unit, data_dir, run_dirs[unit.unit_id], options.python, options.train_script,
                options.table_cache))
        return {"units": len(units), "skipped": len(status), "planned": len(planned)}

    _write_status(status_path, list(status.values()))
    pending = [unit for unit in units if unit.unit_id not in status]
    log.info(
        "queue: %d unit(s), %d complete, %d pending, workers=%d",
        len(units), len(status), len(pending), options.workers,
    )
    lock = threading.Lock()

    def _launch(unit: QueueUnit) -> None:
        run_dir = run_dirs[unit.unit_id]
        run_dir.mkdir(parents=True, exist_ok=True)
        logs_dir.mkdir(parents=True, exist_ok=True)
        log_file = logs_dir / f"{unit.unit_id}.log"
        command = _unit_command(
            unit, data_dir, run_dir, options.python, options.train_script,
            options.table_cache,
        )
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(options.omp_threads)
        env["MKL_NUM_THREADS"] = str(options.omp_threads)
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        started_utc = _now()
        started = time.perf_counter()
        record: dict[str, Any] = {
            "unit_id": unit.unit_id, "stage": unit.stage, "status": "failed",
            "exit_code": -1, "started_utc": started_utc, "finished_utc": "",
            "wall_s": 0.0, "run_dir": str(run_dir), "error": "",
        }
        try:
            if not options.force:
                complete, reason = _complete(
                    unit, run_dir, expected[unit.unit_id]
                )
                if complete:
                    record.update(status="skipped_complete", exit_code=0, error=reason)
                    return
            with open(log_file, "w", encoding="utf-8") as handle:
                handle.write("$ " + " ".join(command) + "\n")
                handle.flush()
                proc = subprocess.run(command, env=env, stdout=handle, stderr=subprocess.STDOUT)
            record["exit_code"] = int(proc.returncode)
            record["status"] = "completed" if proc.returncode == 0 else "failed"
            if proc.returncode != 0:
                record["error"] = f"exit code {proc.returncode}; see {log_file}"
        except OSError as exc:  # pragma: no cover - defensive
            record["error"] = f"{type(exc).__name__}: {exc}"
        finally:
            record["finished_utc"] = _now()
            record["wall_s"] = round(time.perf_counter() - started, 3)
            with lock:
                status[unit.unit_id] = dict(record)
                _write_status(status_path, list(status.values()))
            log.info(
                "unit %-40s %-8s exit=%s wall=%.1fs",
                unit.unit_id, record["status"], record["exit_code"], record["wall_s"],
            )

    if pending:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        with ThreadPoolExecutor(max_workers=int(options.workers)) as pool:
            futures = {pool.submit(_launch, unit): unit for unit in pending}
            for future in as_completed(futures):
                future.result()

    rows = list(status.values())
    failed = [row for row in rows if row["status"] == "failed"]
    summary = {
        "units": len(units),
        "completed": sum(1 for row in rows if row["status"] == "completed"),
        "skipped_complete": sum(1 for row in rows if row["status"] == "skipped_complete"),
        "failed": len(failed),
        "status_csv": str(status_path),
        "failed_units": [row["unit_id"] for row in failed],
    }
    log.info(
        "queue done: %d completed, %d skipped (already complete), %d failed; status %s",
        summary["completed"], summary["skipped_complete"], summary["failed"], status_path,
    )
    return summary
