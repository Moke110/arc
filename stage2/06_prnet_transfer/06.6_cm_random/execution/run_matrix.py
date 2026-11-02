#!/usr/bin/env python
"""Ticket-#6 queue wrapper: run the 45-unit matrix and record peaks/wall time.

    python execution/run_matrix.py --scope "systemd-run --user --scope ..." \\
        --workers 2 --omp-threads 2

What it adds around the ticket-#4 ``run_queue.py`` (which stays the only queue
implementation):

1. **Pre-flight** — refuse to start when ``MemAvailable`` is under ``--free-gib-min``
   (the ticket's 12 GiB rule) or when two other heavy GPU jobs are already running;
   record both readings.
2. **Sampling** — 1 Hz host process-tree RSS, scope cgroup memory and GPU used
   memory, so the record carries measured peaks, not estimates.
3. **Run-level record** — ``<out>/record/matrix_run.json``: wall, peaks, scope,
   pre-flight readings, queue summary and per-unit status, next to the queue's own
   ``queue_status.csv``.
4. **Resume check** (``--verify-resume``) — re-run the queue over the finished
   matrix and assert all 45 units are skipped as complete; then copy one finished
   run, delete its manifest/checkpoint to simulate an interrupted unit, and prove a
   real re-run completes it (and a second pass then skips it).

Optional ``--verify-resume`` is a separate mode and does not touch the matrix.
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

import common  # noqa: E402

RUN_SCHEMA = "cm_random.matrix_run/v1"
RESUME_SCHEMA = "cm_random.resume_check/v1"


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
    """Compute/Graphics processes holding GPU memory (``nvidia-smi`` accounting)."""
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
        rows.append({"pid": int(parts[0]), "used_mib": float(parts[1])})
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
            # comm may contain spaces/parentheses; fields[i] is /proc stat field i+3.
            # field 24 (rss) is in pages, so fields[21] * page size gives KiB.
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
        super().__init__(name="cm-random-sampler", daemon=True)
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
        # sample immediately so a fast failure still has a start reading
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
    """Invoke the ticket-#4 ``run_queue.py`` CLI; return (exit, summary-json, wall)."""
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
    # run_queue.py prints its summary JSON last; recover it from stdout or the log.
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


def _write_matrix_run(out_dir: Path, payload: dict[str, Any]) -> Path:
    path = out_dir / "record" / "matrix_run.json"
    common.write_json_atomic(path, payload)
    return path


def _per_unit_rows(status_csv: Path) -> list[dict[str, Any]]:
    if not status_csv.exists():
        return []
    rows: list[dict[str, Any]] = []
    with open(status_csv, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(dict(row))
    return rows


# --------------------------------------------------------------------------- #
# modes
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

    available = mem_available_gib()
    gpu = gpu_reading()
    other_gpu = [row for row in gpu_processes() if row["used_mib"] >= float(ns.heavy_gpu_mib)]
    preflight = {
        "created_utc": now_utc(),
        "host": socket.gethostname(),
        "nproc": os.cpu_count(),
        "python": sys.executable,
        "queue": str(queue_path),
        "queue_sha256": common.sha256_file(queue_path),
        "units": n_units,
        "workers": int(ns.workers),
        "omp_threads": int(ns.omp_threads),
        "scope": ns.scope,
        "mem_available_gib": round(available, 2),
        "gpu": gpu,
        "other_gpu_processes": other_gpu,
        "heavy_gpu_threshold_mib": float(ns.heavy_gpu_mib),
        "free_gib_min": float(ns.free_gib_min),
        "cuda_check": None,
    }
    # The queue's own interpreter must be the pinned CUDA env; prove it.
    check = subprocess.run(
        [sys.executable, "-c", "import torch;print(int(torch.cuda.is_available()))"],
        capture_output=True, text=True, check=False,
    )
    preflight["cuda_check"] = {
        "returncode": int(check.returncode),
        "stdout": check.stdout.strip()[:200],
        "stderr": check.stderr.strip()[-500:],
    }
    if available < float(ns.free_gib_min):
        print(f"abort: MemAvailable {available:.1f} GiB < {ns.free_gib_min} GiB rule", file=sys.stderr)
        _write_matrix_run(out_dir, {
            "schema": RUN_SCHEMA, "created_utc": now_utc(), "status": "aborted_preflight",
            "reason": f"MemAvailable {available:.2f} GiB below the {ns.free_gib_min} GiB floor",
            "preflight": preflight,
        })
        return 3
    if len(other_gpu) >= int(ns.max_other_heavy_gpu):
        print(
            f"abort: {len(other_gpu)} other heavy GPU process(es) running "
            f"(>= {ns.heavy_gpu_mib} MiB); rule allows {ns.max_other_heavy_gpu}",
            file=sys.stderr,
        )
        _write_matrix_run(out_dir, {
            "schema": RUN_SCHEMA, "created_utc": now_utc(), "status": "aborted_preflight",
            "reason": "another heavy GPU job is already running",
            "preflight": preflight,
        })
        return 3
    if preflight["cuda_check"]["stdout"] != "1":
        print("abort: queue interpreter has no usable CUDA", file=sys.stderr)
        _write_matrix_run(out_dir, {
            "schema": RUN_SCHEMA, "created_utc": now_utc(), "status": "aborted_preflight",
            "reason": f"CUDA unavailable to {sys.executable}",
            "preflight": preflight,
        })
        return 3

    print(json.dumps({"preflight": preflight}, indent=2)[:4000])
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
    record = {
        "schema": RUN_SCHEMA,
        "created_utc": now_utc(),
        "status": "dry_run" if ns.dry_run else ("completed" if code == 0 else "failed"),
        "preflight": preflight,
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
    }
    path = _write_matrix_run(out_dir, record)
    print(json.dumps({
        "matrix_run": str(path),
        "status": record["status"],
        "wall_s": record["run"]["wall_s"],
        "peaks": peaks,
        "queue_summary": summary,
    }, indent=2))
    if ns.dry_run:
        return 0
    # The record writer and checker own candidate-level verdicts; this wrapper's
    # exit code only says "no unit crashed".
    return code


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
        "n_units": len(units),
        "checks": [],
    }

    def check(cid: str, ok: bool, detail: Any) -> None:
        report["checks"].append({"id": cid, "ok": bool(ok), "detail": detail})
        print(f"  [{'ok' if ok else 'FAIL'}] {cid}: {json.dumps(detail)[:600]}")

    # --- 1. the whole finished matrix is skipped as complete ------------------
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

    # --- 2. an interrupted unit re-runs, then is skipped ----------------------
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
        # Simulate an interrupted run: no manifest, no fit history, no checkpoint.
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
        # --- 3. the recovered unit is skipped on the next pass ---------------
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
        description="Run the ticket-#6 45-candidate matrix through the ticket-#4 queue, with peaks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=common.DEFAULT_DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=common.DEFAULT_OUTPUT)
    parser.add_argument("--queue", type=Path, default=None)
    parser.add_argument("--table-cache", type=Path, default=None,
                        help="derived CSV table cache for the pinned env (default: none)")
    parser.add_argument("--workers", type=int, default=2, help="concurrent training processes (1..4)")
    parser.add_argument("--omp-threads", type=int, default=2)
    parser.add_argument("--expect-units", type=int, default=45)
    parser.add_argument("--free-gib-min", type=float, default=12.0,
                        help="abort when MemAvailable is below this many GiB")
    parser.add_argument("--heavy-gpu-mib", type=float, default=1024.0,
                        help="a GPU process using at least this much memory counts as heavy")
    parser.add_argument("--max-other-heavy-gpu", type=int, default=1,
                        help="abort when more than this many other heavy GPU processes run")
    parser.add_argument("--interval-s", type=float, default=1.0, help="sampler period")
    parser.add_argument("--scope", default="", help="the systemd-run scope string, recorded verbatim")
    parser.add_argument("--dry-run", action="store_true", help="ask the queue for its plan only")
    parser.add_argument("--verify-resume", action="store_true",
                        help="run the resume proof instead of the matrix")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    if ns.verify_resume:
        return verify_resume(ns)
    return run_matrix(ns)


if __name__ == "__main__":
    raise SystemExit(main())
