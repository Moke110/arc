#!/usr/bin/env python
"""Ticket-#7 queue wrapper: run the 45-unit Tahoe-init matrix, or a bounded profile.

    # the matrix (after the Lead-run formal pretrain sweep has finished)
    python execution/run_matrix.py --scope "systemd-run --user --scope ..." \
        --workers 2 --omp-threads 2

    # bounded wrapper profile against the first available formal checkpoint
    python execution/run_matrix.py --profile --profile-epochs 2 --workers 2

What it adds around the ticket-#4 ``run_queue.py`` (which stays the only queue
implementation):

1. **Pre-flight** — refuse to start when ``MemAvailable`` is under ``--free-gib-min``
   (the ticket's 12 GiB rule) or when more than one other heavy GPU process is
   already running (the repo-wide "two heavy jobs" rule; a GPU process whose
   memory this WSL2 driver reports as ``[N/A]`` counts as heavy because it cannot
   be shown to be light); record both readings.
2. **Sampling** — 1 Hz host process-tree RSS, scope cgroup memory and GPU used
   memory, so the record carries measured peaks, not estimates.
3. **Run-level record** — ``<out>/record/matrix_run.json``: wall, peaks, scope,
   pre-flight readings, the verified warm-start provenance digest, queue summary
   and per-unit status, next to the queue's own ``queue_status.csv``.
4. **Profile mode** (``--profile``) — builds a one-unit queue with a bounded
   ``--epochs`` budget, runs it through the same queue path, records the raw
   warm-start/startup log lines, and proves the training path is identical to
   ticket #6's apart from the initialisation.  Written to ``<out>/profile/`` and
   never accepted as the formal matrix.
5. **Resume check** (``--verify-resume``) — re-run the queue over the finished
   matrix and assert all 45 units are skipped as complete; then copy one finished
   run, delete its manifest/checkpoint to simulate an interrupted unit, and prove
   a real re-run completes it (and a second pass then skips it).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import socket
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_queue  # noqa: E402
import common  # noqa: E402
import matrix as matrix_mod  # noqa: E402
import warmstart  # noqa: E402

from ptlib import config as pt_config  # noqa: E402
from ptlib import queue as pt_queue  # noqa: E402
from ptlib import train as pt_train  # noqa: E402

RUN_SCHEMA = "cm_tahoe.matrix_run/v1"
PROFILE_SCHEMA = "cm_tahoe.profile/v1"
RESUME_SCHEMA = "cm_tahoe.resume_check/v1"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# --------------------------------------------------------------------------- #
# host readings
# --------------------------------------------------------------------------- #
def meminfo() -> dict[str, int]:
    out: dict[str, int] = {}
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        parts = value.split()
        if parts and parts[0].isdigit():
            out[key] = int(parts[0])  # kB
    return out


def mem_available_gib() -> float:
    return meminfo().get("MemAvailable", 0) / 1024.0 / 1024.0


def cgroup_memory_current_bytes() -> int | None:
    try:
        lines = Path("/proc/self/cgroup").read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    rel = None
    for line in lines:
        if line.startswith("0::"):
            rel = line[3:].strip()
            break
    if not rel:
        return None
    path = Path("/sys/fs/cgroup") / rel.lstrip("/") / "memory.current"
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except (OSError, ValueError):
        return None


def gpu_reading() -> dict[str, Any]:
    try:
        query = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10, check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:  # pragma: no cover - host guard
        return {"available": False, "error": f"{type(exc).__name__}: {exc}"}
    if query.returncode != 0:
        return {"available": False, "error": query.stderr.strip()[:200]}
    rows = [row.strip() for row in query.stdout.splitlines() if row.strip()]
    if not rows:
        return {"available": False, "error": "nvidia-smi returned no GPU"}
    name, used, total = [part.strip() for part in rows[0].split(",")]
    return {
        "available": True,
        "name": name,
        "used_mib": float(used),
        "total_mib": float(total),
        "n_gpus": len(rows),
    }


def gpu_processes() -> list[dict[str, Any]]:
    """Compute/graphics processes holding GPU memory; ``[N/A]`` stays unattributable.

    This WSL2 driver reports ``[N/A]`` for ``used_memory`` (documented by ticket
    #5's profile), so a process with no attributable number is recorded as
    ``used_mib: None`` and treated as heavy by the pre-flight rule.
    """
    try:
        query = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=pid,used_memory",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10, check=False,
        )
    except (OSError, subprocess.SubprocessError):  # pragma: no cover - host guard
        return []
    rows: list[dict[str, Any]] = []
    if query.returncode != 0:
        return rows
    for line in query.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 2 or not parts[0].isdigit():
            continue
        try:
            used_mib: float | None = float(parts[1])
        except ValueError:
            used_mib = None
        rows.append({
            "pid": int(parts[0]),
            "used_mib": used_mib,
            "memory_attributable": used_mib is not None,
        })
    return rows


def process_tree_rss_mib(root_pid: int | None = None) -> float:
    """Sum of RSS over ``root_pid`` and all its descendants (``/proc`` walk)."""
    root_pid = int(os.getpid() if root_pid is None else root_pid)
    ppid: dict[int, int] = {}
    rss_kb: dict[int, int] = {}
    for entry in os.scandir("/proc"):
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
            fields = stat[stat.rindex(")") + 2:].split()
            ppid[pid] = int(fields[1])
            rss_kb[pid] = int(fields[21]) * (os.sysconf("SC_PAGE_SIZE") // 1024)
        except (OSError, ValueError, IndexError):
            continue
    children: dict[int, list[int]] = {}
    for pid, parent in ppid.items():
        children.setdefault(parent, []).append(pid)
    total = 0
    stack = [root_pid]
    seen: set[int] = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += rss_kb.get(pid, 0)
        stack.extend(children.get(pid, ()))
    return total / 1024.0


class Sampler(threading.Thread):
    """1 Hz peak sampler: process-tree RSS, scope cgroup memory, GPU memory."""

    def __init__(self, interval: float = 1.0, include_processes: bool = True) -> None:
        super().__init__(name="cm-tahoe-sampler", daemon=True)
        self.interval = float(interval)
        self.include_processes = include_processes
        self._stop_event = threading.Event()
        self.samples: list[dict[str, Any]] = []
        self.errors: list[str] = []

    def sample_once(self) -> dict[str, Any]:
        reading: dict[str, Any] = {
            "t": round(time.time(), 3),
            "rss_mib": round(process_tree_rss_mib(), 1),
            "cgroup_mib": None,
        }
        current = cgroup_memory_current_bytes()
        if current is not None:
            reading["cgroup_mib"] = round(current / 1024.0 / 1024.0, 1)
        if self.include_processes:
            gpu = gpu_reading()
            if gpu.get("available"):
                reading["gpu_mib"] = gpu["used_mib"]
        self.samples.append(reading)
        return reading

    def run(self) -> None:  # pragma: no cover - thread body
        try:
            self.sample_once()
        except Exception as exc:  # noqa: BLE001 - sampling must never kill the run
            self.errors.append(f"{type(exc).__name__}: {exc}")
        while not self._stop_event.wait(self.interval):
            try:
                self.sample_once()
            except Exception as exc:  # noqa: BLE001
                self.errors.append(f"{type(exc).__name__}: {exc}")

    def stop(self) -> None:
        self._stop_event.set()

    def peaks(self) -> dict[str, Any]:
        rss = [s["rss_mib"] for s in self.samples if s.get("rss_mib") is not None]
        cgroup = [s["cgroup_mib"] for s in self.samples if s.get("cgroup_mib") is not None]
        gpu = [s["gpu_mib"] for s in self.samples if s.get("gpu_mib") is not None]
        return {
            "n_samples": len(self.samples),
            "host_process_tree_rss_mib": max(rss) if rss else None,
            "scope_cgroup_memory_current_mib": max(cgroup) if cgroup else None,
            "gpu_used_mib": max(gpu) if gpu else None,
            "host_rss_start_mib": rss[0] if rss else None,
            "errors": self.errors[:10],
        }


# --------------------------------------------------------------------------- #
# queue invocation
# --------------------------------------------------------------------------- #
def run_queue_cli(
    *,
    python: str,
    queue_path: Path,
    out_root: Path,
    data_dir: Path,
    status_csv: Path,
    workers: int,
    omp_threads: int,
    table_cache: Path | None,
    dry_run: bool = False,
    log_path: Path | None = None,
) -> tuple[int, dict[str, Any], float]:
    """Invoke the ticket-#4 ``run_queue.py`` CLI; return (exit, record-json, wall)."""
    command = [
        python, str(common.PIPELINE_EXEC / "run_queue.py"),
        "--queue", str(queue_path),
        "--out-root", str(out_root),
        "--data-dir", str(data_dir),
        "--status-csv", str(status_csv),
        "--workers", str(int(workers)),
        "--omp-threads", str(int(omp_threads)),
        "--python", python,
    ]
    if table_cache is not None:
        command.extend(["--table-cache", str(table_cache)])
    if dry_run:
        command.append("--dry-run")
    started = time.perf_counter()
    if log_path is None:
        proc = subprocess.run(command, capture_output=True, text=True, check=False)
        stdout, stderr = proc.stdout, proc.stderr
    else:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w", encoding="utf-8") as handle:
            handle.write("$ " + " ".join(command) + "\n\n")
            handle.flush()
            proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
        stdout, stderr = "", ""
    wall = time.perf_counter() - started
    summary: dict[str, Any] = {}
    text = stdout
    if not text and log_path is not None and log_path.exists():
        text = log_path.read_text(encoding="utf-8", errors="replace")
    if text:
        start = text.rfind("\n{")
        blob = text[start + 1:] if start >= 0 else text[text.find("{"):]
        try:
            summary = json.loads(blob)
        except json.JSONDecodeError:
            summary = {}
    record = {
        "command": command,
        "exit_code": int(proc.returncode),
        "wall_s": round(wall, 3),
        "summary": summary,
        "stderr_tail": (stderr or "")[-2000:] if stderr else "",
    }
    return int(proc.returncode), record, wall


def _per_unit_rows(status_csv: Path) -> list[dict[str, Any]]:
    if not status_csv.exists():
        return []
    with open(status_csv, newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


# --------------------------------------------------------------------------- #
# pre-flight
# --------------------------------------------------------------------------- #
def preflight(ns: argparse.Namespace, *, n_units: int, queue_path: Path) -> tuple[dict[str, Any], str | None]:
    """Read the host, run the CUDA check, and return (record, abort-reason)."""
    available = mem_available_gib()
    gpu = gpu_reading()
    processes = gpu_processes()
    other_heavy = [
        row for row in processes
        if row["used_mib"] is None or row["used_mib"] >= float(ns.heavy_gpu_mib)
    ]
    record = {
        "created_utc": now_utc(),
        "host": socket.gethostname(),
        "nproc": os.cpu_count(),
        "python": sys.executable,
        "queue": str(queue_path),
        "queue_sha256": common.sha256_file(queue_path),
        "units": int(n_units),
        "workers": int(ns.workers),
        "omp_threads": int(ns.omp_threads),
        "scope": ns.scope,
        "mem_available_gib": round(available, 2),
        "gpu": gpu,
        "gpu_processes": processes,
        "gpu_memory_attributable": all(row["memory_attributable"] for row in processes),
        "other_heavy_gpu_processes": other_heavy,
        "heavy_gpu_threshold_mib": float(ns.heavy_gpu_mib),
        "max_other_heavy_gpu": int(ns.max_other_heavy_gpu),
        "free_gib_min": float(ns.free_gib_min),
        "cuda_check": None,
    }
    check = subprocess.run(
        [sys.executable, "-c", "import torch;print(int(torch.cuda.is_available()))"],
        capture_output=True, text=True, check=False,
    )
    record["cuda_check"] = {
        "returncode": int(check.returncode),
        "stdout": check.stdout.strip()[:200],
        "stderr": check.stderr.strip()[-500:],
    }
    if available < float(ns.free_gib_min):
        return record, f"MemAvailable {available:.2f} GiB below the {ns.free_gib_min} GiB floor"
    if len(other_heavy) > int(ns.max_other_heavy_gpu):
        return record, (
            f"{len(other_heavy)} other heavy GPU process(es) running "
            f"(> --max-other-heavy-gpu {ns.max_other_heavy_gpu})"
        )
    if record["cuda_check"]["stdout"] != "1":
        return record, f"queue interpreter {sys.executable} has no usable CUDA"
    return record, None


def _write_abort(out_dir: Path, record: dict[str, Any], reason: str, schema: str) -> None:
    common.write_json_atomic(
        out_dir / ("profile" if schema == PROFILE_SCHEMA else "record") / "matrix_run.json",
        {
            "schema": schema,
            "created_utc": now_utc(),
            "status": "aborted_preflight",
            "reason": reason,
            "preflight": record,
        },
    )


# --------------------------------------------------------------------------- #
# formal matrix
# --------------------------------------------------------------------------- #
def run_matrix(ns: argparse.Namespace) -> int:
    out_dir: Path = ns.out_dir
    queue_path: Path = ns.queue or (out_dir / "units.json")
    out_root: Path = out_dir / "runs"
    status_csv: Path = out_dir / "queue_status.csv"
    table_cache = None if ns.table_cache is None else Path(ns.table_cache)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not queue_path.exists():
        print(f"error: queue file not found: {queue_path} (run build_queue.py first)", file=sys.stderr)
        return 2
    queue_payload = common.load_json(queue_path)
    n_units = len(queue_payload.get("units") or [])
    if n_units != int(ns.expect_units):
        print(f"error: queue holds {n_units} units, expected {ns.expect_units}", file=sys.stderr)
        return 2
    if not queue_payload.get("complete", True):
        print("error: queue is marked partial; the formal matrix needs all 45 units", file=sys.stderr)
        return 2
    warm = queue_payload.get("warmstart") or {}
    if warm.get("all_verified") is not True:
        print("error: queue holds no verified warm-start provenance", file=sys.stderr)
        return 2

    record, abort_reason = preflight(ns, n_units=n_units, queue_path=queue_path)
    if abort_reason is not None:
        print(f"abort: {abort_reason}", file=sys.stderr)
        _write_abort(out_dir, record, abort_reason, RUN_SCHEMA)
        return 3
    print(json.dumps({"preflight": record}, indent=2)[:4000])
    sys.stdout.flush()

    sampler = Sampler(interval=float(ns.interval_s))
    sampler.start()
    started = time.perf_counter()
    try:
        code, queue_record, wall = run_queue_cli(
            python=sys.executable,
            queue_path=queue_path,
            out_root=out_root,
            data_dir=Path(ns.data_dir),
            status_csv=status_csv,
            workers=int(ns.workers),
            omp_threads=int(ns.omp_threads),
            table_cache=table_cache,
            dry_run=bool(ns.dry_run),
            log_path=(out_dir / "logs" / "run_queue.log") if not ns.dry_run else None,
        )
    finally:
        wall = time.perf_counter() - started
        sampler.stop()
        sampler.join(timeout=10)

    peaks = sampler.peaks()
    summary = queue_record.get("summary") or {}
    per_unit = _per_unit_rows(status_csv)
    payload = {
        "schema": RUN_SCHEMA,
        "created_utc": now_utc(),
        "status": "dry_run" if ns.dry_run else ("completed" if code == 0 else "failed"),
        "preflight": record,
        "warmstart": {
            "provenance_sha256": warm.get("provenance_sha256"),
            "seeds": warm.get("seeds"),
            "checkpoint_sha256": {
                seed: item.get("sha256") for seed, item in (warm.get("checkpoints") or {}).items()
            },
            "checkpoint_config_sha256": {
                seed: item.get("checkpoint_config_sha256")
                for seed, item in (warm.get("checkpoints") or {}).items()
            },
        },
        "run": {
            "wall_s": round(wall, 3),
            "exit_code": code,
            "queue_summary": summary,
            "queue_status_csv": str(status_csv),
            "queue_log": str(out_dir / "logs" / "run_queue.log"),
            "n_unit_rows": len(per_unit),
            "units_failed": [row.get("unit_id") for row in per_unit if row.get("status") == "failed"],
        },
        "peaks": peaks,
        "git": common.pt_manifest.git_info(common.REPO_ROOT),
        "data_dir": str(ns.data_dir),
        "checkpoint_root": queue_payload.get("checkpoint_root"),
    }
    path = out_dir / "record" / "matrix_run.json"
    common.write_json_atomic(path, payload)
    print(json.dumps({
        "matrix_run": str(path),
        "status": payload["status"],
        "wall_s": payload["run"]["wall_s"],
        "peaks": peaks,
        "queue_summary": summary,
    }, indent=2))
    if ns.dry_run:
        return 0
    return code


# --------------------------------------------------------------------------- #
# bounded profile (wrapper validation, never the matrix)
# --------------------------------------------------------------------------- #
def _arg_flags(args: list[str]) -> dict[str, list[str]]:
    """Parse ``--flag value`` pairs into ``{flag: [values]}`` (value may be missing)."""
    out: dict[str, list[str]] = {}
    index = 0
    while index < len(args):
        token = args[index]
        if token.startswith("--"):
            if index + 1 < len(args) and not args[index + 1].startswith("--"):
                out.setdefault(token, []).append(args[index + 1])
                index += 2
            else:
                out.setdefault(token, []).append("")
                index += 1
        else:
            index += 1
    return out


def expected_config(
    data_dir: Path, *, fold: int, seed: int, lr: float, init_checkpoint: Path
) -> dict[str, Any]:
    """Resolve the **formal** (no profile budget override) run config for one unit."""
    locked, _notes = common.read_locked_config(data_dir)
    argv = [
        "--stage", "finetune",
        "--out", "/nonexistent/placeholder",
        "--fold", str(int(fold)),
        "--seed", str(int(seed)),
        "--lr", repr(float(lr)),
        "--init", matrix_mod.TRAIN_INIT_FLAG,
        "--init-checkpoint", str(init_checkpoint),
        "--data-dir", str(data_dir),
    ]
    options = pt_train.options_from_args(argv)
    cfg, notes = pt_train.resolve_run_config(locked, options)
    return {"config": cfg, "notes": notes}


def random_reference(unit: dict[str, Any], random_runs_root: Path) -> dict[str, Any]:
    """Load ticket #6's run for the same ``(fold, seed, lr)`` and compare paths."""
    random_unit_id = (
        f"cm_f{int(unit['fold'])}_s{int(unit['seed'])}_lr{matrix_mod.lr_label(float(unit['lr']))}_random"
    )
    run_dir = random_runs_root / random_unit_id
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.exists():
        return {
            "available": False,
            "random_unit_id": random_unit_id,
            "random_runs_root": str(random_runs_root),
            "note": "ticket #6's run manifest is not present in this checkout",
        }
    manifest = common.load_json(manifest_path)
    return {
        "available": True,
        "random_unit_id": random_unit_id,
        "random_run_dir": str(run_dir),
        "random_manifest_sha256": common.sha256_file(manifest_path),
        "random_config": manifest.get("config"),
        "random_config_sha256": manifest.get("config_sha256"),
        "random_command": manifest.get("command"),
        "random_init": manifest.get("init"),
        "random_checkpoint_sha256": (manifest.get("checkpoint") or {}).get("sha256"),
    }


def compare_training_path(
    unit: dict[str, Any],
    *,
    data_dir: Path,
    random_runs_root: Path,
) -> dict[str, Any]:
    """Prove the Tahoe path is #6's path with only the initialisation changed."""
    reference = random_reference(unit, random_runs_root)
    expected = expected_config(
        data_dir,
        fold=int(unit["fold"]),
        seed=int(unit["seed"]),
        lr=float(unit["lr"]),
        init_checkpoint=Path(unit["init_checkpoint"]),
    )
    tahoe_config = expected["config"]
    tahoe_sha = pt_config.config_fingerprint(tahoe_config)
    result: dict[str, Any] = {
        "ticket6_reference": reference,
        "tahoe_expected_config": tahoe_config,
        "tahoe_expected_config_sha256": tahoe_sha,
        "tahoe_unit_args": list(unit["args"]),
        "expected_config_notes": expected["notes"],
    }
    if not reference["available"]:
        result.update({"comparable": False, "only_init_differs": None})
        return result
    random_config = reference["random_config"] or {}
    keys = set(tahoe_config) | set(random_config)
    differing = sorted(key for key in keys if tahoe_config.get(key) != random_config.get(key))
    init_pair = {key: {"random": random_config.get(key), "tahoe": tahoe_config.get(key)} for key in differing}
    # Unit args: every flag except the initialisation ones must be identical.
    # The initialisation flags are the intended difference; the paths and the
    # profile's bounded-budget overrides are not part of the formal training path.
    shared_flags = {"--init", "--init-checkpoint", "--out", "--data-dir",
                    "--hash-cache", "--table-cache", "--epochs", "--min-epochs"}
    tahoe_flags = {k: v for k, v in _arg_flags(unit["args"]).items() if k not in shared_flags}
    random_flags_all = _arg_flags(list(reference.get("random_command") or []))
    random_flags = {key: random_flags_all.get(key) for key in tahoe_flags}
    shared_args_identical = all(
        random_flags_all.get(flag) == values for flag, values in tahoe_flags.items()
    )
    only_init = differing == ["init"] and init_pair["init"] == {"random": "random", "tahoe": "pretrained"}
    result.update({
        "comparable": True,
        "differing_config_keys": differing,
        "differing_config_values": {
            key: {"random": random_config.get(key), "tahoe": tahoe_config.get(key)}
            for key in differing
        },
        "expected_init_only_difference": init_pair,
        "only_init_differs": bool(only_init),
        "shared_unit_args_identical": bool(shared_args_identical),
        "tahoe_shared_flags": tahoe_flags,
        "random_shared_flags": random_flags,
    })
    result["passed"] = bool(only_init and shared_args_identical)
    return result


def _raw_startup(log_path: Path, *, warm_lines: int = 5, head_lines: int = 25) -> dict[str, Any]:
    """The raw startup/axis-check output of the actual training subprocess."""
    if not log_path.exists():
        return {"log": str(log_path), "available": False}
    text = log_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    return {
        "log": str(log_path),
        "available": True,
        "head": lines[:head_lines],
        "warm_start_lines": [line for line in lines if "warm start" in line][:warm_lines],
        "config_lines": [line for line in lines if "resolved" in line.lower()][:5],
        "n_lines": len(lines),
    }


def verify_profile_candidate(
    unit: dict[str, Any],
    *,
    out_root: Path,
    expected_checkpoint: dict[str, Any],
    profile_epochs: int,
) -> tuple[dict[str, Any], list[str]]:
    run_dir = out_root / unit["unit_id"]
    manifest_path = run_dir / "manifest.json"
    problems: list[str] = []
    if not manifest_path.exists():
        return {"run_dir": str(run_dir), "manifest": None}, ["profile run produced no manifest"]
    manifest = common.load_json(manifest_path)
    init = manifest.get("init")
    problems.extend(warmstart.verify_manifest_init(init, expected_checkpoint))
    if manifest.get("status") != "completed":
        problems.append(f"profile run manifest status is {manifest.get('status')!r}")
    if int(manifest.get("epochs_run", -1)) != int(profile_epochs):
        problems.append(
            f"profile ran {manifest.get('epochs_run')} epochs, expected {profile_epochs}"
        )
    if int((manifest.get("config") or {}).get("recipe", {}).get("max_epochs", -1)) != int(profile_epochs):
        problems.append("profile run config.recipe.max_epochs differs from the bounded budget")
    if manifest.get("smoke") is not False:
        problems.append("profile run was marked smoke; expected a real bounded run")
    if (manifest.get("checkpoint") or {}).get("reload_identical") is not True:
        problems.append("profile run checkpoint did not strict-reload identically")
    candidate = {
        "run_dir": str(run_dir),
        "manifest": str(manifest_path),
        "manifest_sha256": common.sha256_file(manifest_path),
        "status": manifest.get("status"),
        "config_sha256": manifest.get("config_sha256"),
        "init": init,
        "best_valid_auprc": manifest.get("best_valid_auprc"),
        "best_epoch": manifest.get("best_epoch"),
        "epochs_run": manifest.get("epochs_run"),
        "stop_reason": manifest.get("stop_reason"),
        "checkpoint": manifest.get("checkpoint"),
        "peak_rss_mb": manifest.get("peak_rss_mb"),
        "timings": manifest.get("timings"),
    }
    return candidate, problems


def _stage_checkpoint(source: Path, staged: Path) -> dict[str, Any]:
    """Byte-identical snapshot of a live checkpoint, verified by sha256.

    A running formal sweep rewrites ``model.pt`` whenever validation improves, so a
    verification/read pair on the live path could straddle a replacement.  Copying
    reads one inode atomically (``save_checkpoint`` uses ``os.replace``), and the
    sha256 of source-before, source-after and the copy must all agree.
    """
    staged.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, 4):
        before = common.sha256_file(source)
        shutil.copy2(source, staged)
        after = common.sha256_file(source)
        staged_sha = common.sha256_file(staged)
        if before == after == staged_sha:
            return {
                "source_path": str(source),
                "staged_path": str(staged),
                "source_sha256_at_copy": before,
                "staged_sha256": staged_sha,
                "attempts": attempt,
                "byte_identical": True,
                "note": "staged snapshot of the live formal checkpoint; sha256-equal to the "
                        "source inode at copy time",
            }
    raise SystemExit(f"could not snapshot {source} consistently after 3 attempts")


def _run_axis_check(
    *, data_dir: Path, checkpoint_root: Path, seed: int, out_path: Path
) -> dict[str, Any]:
    """Run ``check_init.py`` as its own process and capture its raw output."""
    command = [
        sys.executable,
        str(Path(__file__).resolve().parent / "check_init.py"),
        "--data-dir", str(data_dir),
        "--checkpoint-root", str(checkpoint_root),
        "--seeds", str(int(seed)),
        "--out", str(out_path),
    ]
    proc = subprocess.run(command, capture_output=True, text=True, check=False)
    raw = proc.stdout
    if proc.stderr:
        raw += f"\n[stderr]\n{proc.stderr}"
    log_path = out_path.with_suffix(".log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("$ " + " ".join(command) + "\n\n" + raw, encoding="utf-8")
    return {
        "command": command,
        "exit_code": int(proc.returncode),
        "raw_output": raw,
        "log": str(log_path),
        "out": str(out_path),
    }


def run_profile(ns: argparse.Namespace) -> int:
    out_dir: Path = ns.out_dir
    profile_dir = out_dir / "profile"
    out_root: Path = out_dir / "runs"
    status_csv = profile_dir / "queue_status.csv"
    profile_dir.mkdir(parents=True, exist_ok=True)

    # Stage a stable snapshot of the first available formal checkpoint; the live
    # sweep may replace it while this profile runs.
    locked, _notes = common.read_locked_config(Path(ns.data_dir))
    spec = matrix_mod.matrix_spec(locked)
    cands = matrix_mod.candidates(spec)
    unit_id = ns.profile_unit or build_queue.default_profile_unit(cands, spec["seeds"][0])
    candidate = next((item for item in cands if item.unit_id == unit_id), None)
    if candidate is None:
        print(f"error: --profile-unit {unit_id!r} is not a matrix unit", file=sys.stderr)
        return 2
    source_checkpoint = matrix_mod.checkpoint_path(Path(ns.checkpoint_root), candidate.seed)
    if not source_checkpoint.is_file():
        print(f"error: no formal checkpoint at {source_checkpoint}", file=sys.stderr)
        return 2
    staged_root = profile_dir / "staged_checkpoints"
    staging = _stage_checkpoint(
        source_checkpoint, matrix_mod.checkpoint_path(staged_root, candidate.seed)
    )
    axis_check = _run_axis_check(
        data_dir=Path(ns.data_dir),
        checkpoint_root=staged_root,
        seed=candidate.seed,
        out_path=profile_dir / f"init_seed{candidate.seed}.json",
    )
    print(axis_check["raw_output"])
    if axis_check["exit_code"] != 0:
        print(f"error: axis/structure check failed for seed {candidate.seed}", file=sys.stderr)
        return 2

    queue_payload = build_queue.build_queue(
        data_dir=Path(ns.data_dir),
        out_dir=out_dir,
        checkpoint_root=staged_root,
        queue_path=profile_dir / "profile_queue.json",
        hash_cache=None if str(ns.hash_cache).lower() in ("none", "off", "") else str(ns.hash_cache),
        seeds=None,
        profile=True,
        profile_unit=ns.profile_unit,
        profile_epochs=int(ns.profile_epochs),
    )
    unit = queue_payload["units"][0]
    expected_checkpoint = queue_payload["warmstart"]["checkpoints"][str(unit["seed"])]
    queue_path = profile_dir / "profile_queue.json"

    record, abort_reason = preflight(ns, n_units=1, queue_path=queue_path)
    if abort_reason is not None:
        print(f"abort: {abort_reason}", file=sys.stderr)
        _write_abort(out_dir, record, abort_reason, PROFILE_SCHEMA)
        return 3
    print(json.dumps({"preflight": record}, indent=2)[:4000])
    print(json.dumps({"warmstart": {
        "seed": unit["seed"],
        "sha256": expected_checkpoint["sha256"],
        "checkpoint_config_sha256": expected_checkpoint["checkpoint_config_sha256"],
        "gene_axis": expected_checkpoint["gene_axis"],
        "load_provenance": {
            key: expected_checkpoint["load_provenance"][key]
            for key in ("strict_load", "weights_identical", "missing_keys", "unexpected_keys",
                        "shape_mismatches", "dtype_mismatches", "n_loaded_tensors",
                        "n_loaded_parameter_tensors", "n_loaded_buffer_tensors")
        },
        "anti_random_init": expected_checkpoint["anti_random_init"],
    }}, indent=2))
    sys.stdout.flush()

    sampler = Sampler(interval=float(ns.interval_s))
    sampler.start()
    started = time.perf_counter()
    try:
        code, queue_record, wall = run_queue_cli(
            python=sys.executable,
            queue_path=queue_path,
            out_root=out_root,
            data_dir=Path(ns.data_dir),
            status_csv=status_csv,
            workers=int(ns.workers),
            omp_threads=int(ns.omp_threads),
            table_cache=None if ns.table_cache is None else Path(ns.table_cache),
            log_path=profile_dir / "logs" / "run_queue.log",
        )
    finally:
        wall = time.perf_counter() - started
        sampler.stop()
        sampler.join(timeout=10)
    peaks = sampler.peaks()

    candidate, candidate_problems = verify_profile_candidate(
        unit,
        out_root=out_root,
        expected_checkpoint=expected_checkpoint,
        profile_epochs=int(ns.profile_epochs),
    )
    path_comparison = compare_training_path(
        unit,
        data_dir=Path(ns.data_dir),
        random_runs_root=Path(ns.random_runs_root),
    )
    startup = _raw_startup(out_root / "queue_logs" / f"{unit['unit_id']}.log")
    problems = list(candidate_problems)
    if code != 0:
        problems.append(f"queue exit code {code}")
    if not path_comparison.get("comparable"):
        problems.append("ticket #6 reference run not found; the training-path comparison could not run")
    elif not path_comparison.get("passed"):
        problems.append("training path differs from ticket #6 beyond the initialisation")
    passed = not problems
    payload = {
        "schema": PROFILE_SCHEMA,
        "created_utc": now_utc(),
        "mode": "profile",
        "command": [sys.executable, *sys.argv],
        "git": common.pt_manifest.git_info(common.REPO_ROOT),
        "data_dir": str(ns.data_dir),
        "checkpoint_root": str(ns.checkpoint_root),
        "random_runs_root": str(ns.random_runs_root),
        "bounded_budget": {
            "epochs": int(ns.profile_epochs),
            "note": "wrapper/profile validation only; never the formal 100-epoch locked budget",
        },
        "unit": {
            "unit_id": unit["unit_id"],
            "fold": unit["fold"],
            "seed": unit["seed"],
            "lr": unit["lr"],
            "lr_label": unit["lr_label"],
            "args": unit["args"],
        },
        "warmstart": {
            "provenance_sha256": queue_payload["warmstart"]["provenance_sha256"],
            "record": expected_checkpoint,
        },
        "preflight": record,
        "run": {
            "wall_s": round(wall, 3),
            "exit_code": code,
            "queue_summary": queue_record.get("summary"),
            "queue_log": str(profile_dir / "logs" / "run_queue.log"),
            "status_csv": str(status_csv),
        },
        "peaks": peaks,
        "candidate": candidate,
        "path_comparison": path_comparison,
        "profile_staging": {
            "checkpoint": staging,
            "axis_check": axis_check,
        },
        "raw_startup": startup,
        "notes": [
            "bounded wrapper profile only: 2 epochs on one unit, never the formal 100-epoch budget",
            "run against the first available #5 formal checkpoint while the seed-20261021 pretrain "
            "was still improving; the checkpoint was snapshotted byte-identically first, so the "
            "live sweep's next os.replace cannot tear the verification/load",
            "the GPU device peak includes the concurrent #5 pretrain; this WSL2 driver reports "
            "[N/A] for per-process GPU memory, so the unit's own GPU footprint is not claimed",
            "the ticket-#6 comparison uses the same (fold, seed, lr) random-init run manifest and "
            "proves the resolved config differs only in `init`",
        ],
        "problems": problems,
        "passed": passed,
    }
    path = profile_dir / "profile_record.json"
    common.write_json_atomic(path, payload)
    print(json.dumps({
        "profile_record": str(path),
        "passed": passed,
        "problems": problems,
        "wall_s": payload["run"]["wall_s"],
        "peaks": peaks,
        "candidate": {
            key: candidate.get(key)
            for key in ("status", "best_valid_auprc", "best_epoch", "epochs_run", "stop_reason")
        },
        "path_comparison": {
            key: path_comparison.get(key)
            for key in ("comparable", "only_init_differs", "differing_config_keys",
                        "shared_unit_args_identical", "passed")
        },
    }, indent=2))
    return 0 if passed else 1


# --------------------------------------------------------------------------- #
# resume proof
# --------------------------------------------------------------------------- #
def verify_resume(ns: argparse.Namespace) -> int:
    """Prove resume: 45/45 skipped when complete, and an interrupted unit re-runs."""
    out_dir: Path = ns.out_dir
    out_root: Path = out_dir / "runs"
    queue_path: Path = ns.queue or (out_dir / "units.json")
    probe_root: Path = out_dir / "_resume_probe"
    probe_root.mkdir(parents=True, exist_ok=True)
    queue_payload = common.load_json(queue_path)
    units = queue_payload["units"]
    report: dict[str, Any] = {
        "schema": RESUME_SCHEMA,
        "created_utc": now_utc(),
        "queue": str(queue_path),
        "warmstart_provenance_sha256": (queue_payload.get("warmstart") or {}).get("provenance_sha256"),
        "n_units": len(units),
        "checks": [],
    }

    def check(cid: str, ok: bool, detail: Any) -> None:
        report["checks"].append({"id": cid, "ok": bool(ok), "detail": detail})
        print(f"  [{'ok' if ok else 'FAIL'}] {cid}: {json.dumps(detail)[:600]}")

    code, record, _wall = run_queue_cli(
        python=sys.executable,
        queue_path=queue_path,
        out_root=out_root,
        data_dir=Path(ns.data_dir),
        status_csv=probe_root / "whole_matrix_status.csv",
        workers=1,
        omp_threads=2,
        table_cache=None if ns.table_cache is None else Path(ns.table_cache),
        dry_run=True,
    )
    summary = record.get("summary") or {}
    check(
        "finished_matrix_all_skipped",
        code == 0
        and int(summary.get("units", -1)) == len(units)
        and int(summary.get("skipped", -1)) == len(units)
        and int(summary.get("planned", -1)) == 0,
        {"exit_code": code, "summary": summary},
    )

    victim = None
    for unit in units:
        if (out_root / unit["unit_id"] / "manifest.json").exists():
            victim = unit
            break
    if victim is None:
        check("interrupted_unit_recoverable", False, "no completed unit to copy")
    else:
        victim_id = victim["unit_id"]
        scratch = probe_root / "runs" / victim_id
        if scratch.exists():
            shutil.rmtree(scratch)
        shutil.copytree(out_root / victim_id, scratch)
        original_manifest = common.load_json(out_root / victim_id / "manifest.json")
        for name in ("manifest.json", "fit_history.json", "summary_table.csv"):
            (scratch / name).unlink(missing_ok=True)
        shutil.rmtree(scratch / "checkpoint", ignore_errors=True)
        single_queue = probe_root / "single_unit.json"
        common.write_json_atomic(
            single_queue,
            {"schema": queue_payload.get("schema"), "units": [victim],
             "data_dir": queue_payload.get("data_dir")},
        )
        code, record, _wall = run_queue_cli(
            python=sys.executable,
            queue_path=single_queue,
            out_root=probe_root / "runs",
            data_dir=Path(ns.data_dir),
            status_csv=probe_root / "interrupted_status.csv",
            workers=1,
            omp_threads=2,
            table_cache=None if ns.table_cache is None else Path(ns.table_cache),
        )
        summary = record.get("summary") or {}
        rerun_manifest_path = scratch / "manifest.json"
        rerun_manifest = (
            common.load_json(rerun_manifest_path) if rerun_manifest_path.exists() else {}
        )
        check(
            "interrupted_unit_recoverable",
            code == 0 and int(summary.get("completed", 0)) == 1
            and rerun_manifest.get("config_sha256") == original_manifest.get("config_sha256"),
            {
                "unit_id": victim_id,
                "exit_code": code,
                "summary": summary,
                "config_sha256_original": original_manifest.get("config_sha256"),
                "config_sha256_rerun": rerun_manifest.get("config_sha256"),
                "rerun_wall_s": record.get("wall_s"),
            },
        )
        code2, record2, _wall2 = run_queue_cli(
            python=sys.executable,
            queue_path=single_queue,
            out_root=probe_root / "runs",
            data_dir=Path(ns.data_dir),
            status_csv=probe_root / "recovered_status.csv",
            workers=1,
            omp_threads=2,
            table_cache=None if ns.table_cache is None else Path(ns.table_cache),
            dry_run=True,
        )
        summary2 = record2.get("summary") or {}
        check(
            "recovered_unit_then_skipped",
            code2 == 0 and int(summary2.get("skipped", -1)) == 1
            and int(summary2.get("planned", -1)) == 0,
            {"exit_code": code2, "summary": summary2},
        )
        report["interrupted_unit"] = {
            "unit_id": victim_id,
            "original_manifest_sha256": common.sha256_file(out_root / victim_id / "manifest.json"),
            "rerun_manifest_sha256": (
                common.sha256_file(rerun_manifest_path) if rerun_manifest_path.exists() else None
            ),
            "rerun_config_sha256": rerun_manifest.get("config_sha256"),
            "original_config_sha256": original_manifest.get("config_sha256"),
        }

    report["passed"] = all(item["ok"] for item in report["checks"])
    path = out_dir / "record" / "resume_check.json"
    common.write_json_atomic(path, report)
    print(json.dumps({"resume_check": str(path), "passed": report["passed"],
                      "checks": report["checks"]}, indent=2)[:4000])
    return 0 if report["passed"] else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the ticket-#7 45-candidate Tahoe-init matrix through the ticket-#4 queue.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--checkpoint-root", type=Path, default=None,
                        help="ticket #5 formal runs root (default: local, else integration)")
    parser.add_argument("--out-dir", type=Path, default=common.DEFAULT_OUTPUT)
    parser.add_argument("--queue", type=Path, default=None)
    parser.add_argument("--table-cache", type=Path, default=None,
                        help="derived CSV table cache for the pinned env (default: none)")
    parser.add_argument("--random-runs-root", type=Path, default=None,
                        help="ticket #6 run tree for the profile training-path comparison")
    parser.add_argument("--hash-cache", type=Path, default=None,
                        help="frozen-artifact hash cache (profile mode)")
    parser.add_argument("--workers", type=int, default=2, help="concurrent training processes (1..4)")
    parser.add_argument("--omp-threads", type=int, default=2)
    parser.add_argument("--expect-units", type=int, default=45)
    parser.add_argument("--free-gib-min", type=float, default=12.0,
                        help="abort when MemAvailable is below this many GiB")
    parser.add_argument("--heavy-gpu-mib", type=float, default=1024.0,
                        help="a GPU process using at least this much memory counts as heavy")
    parser.add_argument("--max-other-heavy-gpu", type=int, default=1,
                        help="abort when more than this many other heavy GPU processes run "
                             "(one is allowed: the repo-wide two-heavy-jobs rule)")
    parser.add_argument("--interval-s", type=float, default=1.0, help="sampler period")
    parser.add_argument("--scope", default="", help="the systemd-run scope string, recorded verbatim")
    parser.add_argument("--dry-run", action="store_true", help="ask the queue for its plan only")
    parser.add_argument("--profile", action="store_true",
                        help="bounded one-unit wrapper profile instead of the matrix")
    parser.add_argument("--profile-unit", default=None,
                        help="profile unit id (default: fold 0, highest LR, first seed)")
    parser.add_argument("--profile-epochs", type=int, default=2,
                        help="bounded epoch budget for --profile (never the locked 100)")
    parser.add_argument("--verify-resume", action="store_true",
                        help="run the resume proof instead of the matrix")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    ns.data_dir = common.resolve_data_dir(ns.data_dir)
    ns.checkpoint_root = common.resolve_checkpoint_root(ns.checkpoint_root)
    ns.random_runs_root = common.resolve_random_runs_root(ns.random_runs_root)
    if ns.hash_cache is None:
        ns.hash_cache = Path(ns.out_dir) / "hash_cache.json"
    if ns.verify_resume:
        return verify_resume(ns)
    if ns.profile:
        return run_profile(ns)
    return run_matrix(ns)


if __name__ == "__main__":
    raise SystemExit(main())
