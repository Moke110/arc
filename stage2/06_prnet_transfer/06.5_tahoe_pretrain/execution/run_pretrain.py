#!/usr/bin/env python
"""Run the ticket-#5 seeded Tahoe pretrains through the 06.4 pipeline.

One invocation = one *sweep*: the frozen seeds are trained **sequentially**
(never concurrently — a formal pretrain is the heaviest job in the repository),
each through ``06.4_pipeline/execution/train.py --stage pretrain`` with the
frozen locked config and **no** budget override.  After every unit the wrapper
reads the pipeline's own ``manifest.json`` / ``fit_history.json`` /
``model_card.md``, verifies the config digest against the plan, measures the
run, and finally writes the substep record (see ``record.py``).

    # formal: the three frozen seeds, full locked budget (100 epochs cap,
    #         20-epoch floor, early stop after 20 stale epochs, lr 1e-3)
    python run_pretrain.py --mode formal

    # what would run?  (no data read, no GPU)
    python run_pretrain.py --mode formal --dry-run

    # one seed of the sweep (the record is written as a *partial* formal record)
    python run_pretrain.py --mode formal --seeds 20261021

    # short real-data profile: the 2-epoch bounded smoke + one full-width epoch
    python run_pretrain.py --mode profile

    # resume a sweep whose first seeds already completed
    python run_pretrain.py --mode formal --resume

The wrapper refuses to start when the frozen input, the derived memmap cache or
the frozen-table cache disagree with the 06.3 contract, when the locked recipe
differs from the ticket's locked recipe, when available RAM is below the
execution-environment floor, or when another pretrain runs in the same VM.
``run.sh`` (and ``make profile`` / ``make formal``) put this process inside a
``systemd-run --user --scope`` memory cap — use that, not a bare invocation.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

EXECUTION = Path(__file__).resolve().parent
STEP = EXECUTION.parent
PIPELINE_EXECUTION = STEP.parent / "06.4_pipeline" / "execution"

# Bind the 06.4 pipeline library (ADR 0032: #5–#11 import ptlib, never copy it).
sys.path.insert(0, str(PIPELINE_EXECUTION))

import plan  # noqa: E402  (needs the sys.path insertion above)
import record  # noqa: E402

FROZEN_INPUT = STEP.parent / "06.3_input" / "output"
PIPELINE_OUTPUT = STEP.parent / "06.4_pipeline" / "output"
DEFAULT_CACHE = PIPELINE_OUTPUT / "cache" / "tahoe_log1p_cp10k_v1.npy"
DEFAULT_TABLE_CACHE = PIPELINE_OUTPUT / "input_cache"
DEFAULT_OUT_ROOT = STEP / "output"

#: execution-environment floors (docs/agents/execution-environment.md)
MIN_AVAILABLE_GIB = 12.0
FORMAL_AVAILABLE_GIB = 16.0

#: RAM caps run.sh applies; recorded for reference.
DEFAULT_CAPS_GIB = {"profile": 12, "formal": 18}


# --------------------------------------------------------------------------- #
# small probes
# --------------------------------------------------------------------------- #
def available_gib() -> float:
    with open("/proc/meminfo", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024.0 / 1024.0
    raise RuntimeError("MemAvailable missing from /proc/meminfo")


def cgroup_info() -> dict[str, Any]:
    """Record the cgroup (and its memory cap) the sweep runs under."""
    info: dict[str, Any] = {"cgroup": None, "memory_max": None, "invocation": None}
    try:
        lines = Path("/proc/self/cgroup").read_text(encoding="utf-8").splitlines()
        for line in lines:
            if line.startswith("0::"):
                info["cgroup"] = line.split("::", 1)[1]
                break
    except OSError:
        return info
    info["invocation"] = os.environ.get("INVOCATION_ID")
    path = info["cgroup"]
    if path:
        try:
            info["memory_max"] = Path("/sys/fs/cgroup" + path + "/memory.max").read_text().strip()
        except OSError:
            info["memory_max"] = None
    return info


def other_pretrain_running() -> list[str]:
    """Command lines of *other* running ``train.py --stage pretrain`` processes."""
    found: list[str] = []
    own = {os.getpid(), os.getppid()}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        if pid in own:
            continue
        try:
            cmdline = (entry / "cmdline").read_bytes().decode("utf-8", "replace")
        except OSError:
            continue
        if "train.py" not in cmdline or "--stage" not in cmdline:
            continue
        parts = cmdline.split("\0")
        if "pretrain" in parts:
            found.append(f"pid {pid}: {' '.join(part for part in parts if part)}")
    return found


def host_block(device: str) -> dict[str, Any]:
    import torch  # noqa: PLC0415 — ptlib is already imported; record the exact stack

    block: dict[str, Any] = {
        "python": sys.executable,
        "python_version": sys.version.split()[0],
        "torch": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "requested_device": device,
        "cwd": str(Path.cwd()),
        "repo_root": str(record.find_repo_root(EXECUTION)),
        "mem_available_gib_at_start": round(available_gib(), 2),
        "cgroup": cgroup_info(),
    }
    if torch.cuda.is_available():
        block["cuda_device_name"] = torch.cuda.get_device_name(0)
        block["cuda_device_total_mib"] = round(
            torch.cuda.get_device_properties(0).total_memory / 1024**2, 0
        )
    return block


# --------------------------------------------------------------------------- #
# external resource sampler (GPU; host RSS cross-check)
# --------------------------------------------------------------------------- #
class ResourceSampler(threading.Thread):
    """Sample the trainer's host high-water RSS and the GPU memory it holds.

    The trainer's own manifest is authoritative for host peak RSS
    (``ru_maxrss`` of the training process); this sampler is the cross-check and
    the only source of the GPU peak (device-level and attributed to the trainer
    pid, so another process on the GPU cannot be mistaken for the run).
    """

    def __init__(self, pid: int, interval_s: float = 1.0) -> None:
        super().__init__(name="06.5-resource-sampler", daemon=True)
        self.pid = int(pid)
        self.interval_s = max(0.2, float(interval_s))
        self.host_peak_mb = 0.0
        self.gpu_device_peak_mib = 0.0
        self.gpu_process_peak_mib = 0.0
        self.gpu_baseline_mib = self._device_used_mib() or 0.0
        self.samples = 0
        self.gpu_process_supported: bool | None = None
        self.errors: list[str] = []
        self._stop_event = threading.Event()

    # -- probes ------------------------------------------------------------
    @staticmethod
    def _nvidia(*args: str) -> list[str]:
        out = subprocess.run(
            ["nvidia-smi", *args], capture_output=True, text=True, timeout=20, check=False
        )
        if out.returncode != 0:
            raise RuntimeError(out.stderr.strip() or f"nvidia-smi exited {out.returncode}")
        return [line.strip() for line in out.stdout.splitlines() if line.strip()]

    def _device_used_mib(self) -> float | None:
        try:
            values = self._nvidia("--query-gpu=memory.used", "--format=csv,noheader,nounits")
            return float(values[0])
        except (OSError, RuntimeError, ValueError, IndexError) as exc:
            self.errors.append(f"gpu device sample failed: {exc}")
            return None

    def _process_used_mib(self) -> float | None:
        """GPU memory held by the trainer pid, or ``None`` when the driver cannot tell.

        Some WSL2 drivers report ``[N/A]`` for ``used_memory`` in
        ``--query-compute-apps``; the device-level peak (minus the baseline) is
        then the attribution, and the record says so.
        """
        try:
            rows = self._nvidia(
                "--query-compute-apps=pid,used_memory", "--format=csv,noheader,nounits"
            )
        except (OSError, RuntimeError) as exc:
            self.errors.append(f"gpu process sample failed: {exc}")
            return None
        parsable = 0
        total = 0.0
        for row in rows:
            parts = [piece.strip() for piece in row.split(",")]
            if len(parts) < 2:
                continue
            try:
                pid = int(parts[0])
                used = float(parts[1])
            except ValueError:
                continue
            parsable += 1
            if pid == self.pid:
                total += used
        if rows and parsable == 0:
            self.gpu_process_supported = False
            return None
        return total

    def _host_peak_mb(self) -> float | None:
        try:
            status = Path(f"/proc/{self.pid}/status").read_text(encoding="utf-8")
        except OSError:
            return None
        for line in status.splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) / 1024.0
        return None

    # -- thread ------------------------------------------------------------
    def run(self) -> None:  # noqa: D401 — thread body
        while not self._stop_event.is_set():
            self.sample_once()
            self._stop_event.wait(self.interval_s)

    def sample_once(self) -> None:
        host = self._host_peak_mb()
        if host is not None:
            self.host_peak_mb = max(self.host_peak_mb, host)
        device = self._device_used_mib()
        if device is not None:
            self.gpu_device_peak_mib = max(self.gpu_device_peak_mib, device)
        process = self._process_used_mib()
        if process is not None:
            self.gpu_process_peak_mib = max(self.gpu_process_peak_mib, process)
        self.samples += 1

    def stop(self) -> dict[str, Any]:
        self._stop_event.set()
        try:
            self.join(timeout=30)
        except RuntimeError:  # pragma: no cover — thread never started
            pass
        after = self._device_used_mib()
        return record.resource_report(
            host_peak_rss_mb_pipeline=None,  # filled from the trainer manifest
            host_peak_rss_mb_sampler=round(self.host_peak_mb, 1) or None,
            gpu_peak_mib_device=round(self.gpu_device_peak_mib, 1) or None,
            gpu_peak_mib_process=round(self.gpu_process_peak_mib, 1) or None,
            gpu_baseline_mib=round(self.gpu_baseline_mib, 1),
            gpu_after_mib=None if after is None else round(after, 1),
            samples=self.samples,
            interval_s=self.interval_s,
            gpu_process_attribution=(
                "nvidia-smi compute-apps pid,used_memory"
                if self.gpu_process_supported
                else "unavailable: the driver reports [N/A] for per-process GPU memory; "
                     "the device peak minus the baseline is the run's GPU peak when the "
                     "run is alone on the GPU"
            ),
        )


# --------------------------------------------------------------------------- #
# preflight: frozen input, caches, locked config, memory, concurrency
# --------------------------------------------------------------------------- #
def _read_sidecar(cache: Path) -> dict[str, Any]:
    sidecar = Path(str(cache) + ".json")
    if not sidecar.exists():
        raise SystemExit(f"derived Tahoe cache sidecar is missing: {sidecar}")
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"derived Tahoe cache sidecar is not valid JSON: {sidecar} ({exc})")
    if payload.get("schema") != "ptlib.tahoe_log1p_cache/v1":
        raise SystemExit(f"derived Tahoe cache sidecar schema {payload.get('schema')!r}, expected "
                         "'ptlib.tahoe_log1p_cache/v1'")
    return payload


def cache_identity(cache: Path, identity: Mapping[str, Any]) -> dict[str, Any]:
    """Verify the derived memmap is the frozen layer of *this* frozen input."""
    cache = Path(cache)
    if not cache.exists():
        raise SystemExit(
            f"derived Tahoe cache is missing: {cache}\n"
            "build it once (one heavy job at a time):\n"
            f"  {sys.executable} {PIPELINE_EXECUTION}/tools/build_cache.py "
            f"--data-dir {identity['data_dir']} --cache {cache}"
        )
    sidecar = _read_sidecar(cache)
    frozen_sha = (identity.get("artifacts") or {}).get("tahoe_transfer.h5ad")
    if sidecar.get("source_sha256") != frozen_sha:
        raise SystemExit(
            "derived Tahoe cache was not built from the recorded frozen tahoe_transfer.h5ad:\n"
            f"  cache source_sha256 = {sidecar.get('source_sha256')}\n"
            f"  frozen sha256       = {frozen_sha}\n"
            "Fix (keeps the bytes, proves equality, costs ~20 s + one 16 GB read):\n"
            f"  {sys.executable} {PIPELINE_EXECUTION}/tools/build_cache.py "
            f"--data-dir {identity['data_dir']} --cache {cache} --rebind-frozen\n"
            "Or rebuild from the frozen input:\n"
            f"  {sys.executable} {PIPELINE_EXECUTION}/tools/build_cache.py "
            f"--data-dir {identity['data_dir']} --cache {cache}"
        )
    n_obs = int((identity.get("counts_csv", {}).get("tahoe") or {}).get("rows", -1))
    n_genes = int((identity.get("gene_axis") or {}).get("n_genes", -1))
    shape = sidecar.get("shape")
    if [int(x) for x in (shape or [])] != [n_obs, n_genes]:
        raise SystemExit(
            f"derived Tahoe cache shape {shape} != frozen view {(n_obs, n_genes)}; rebuild it"
        )
    if sidecar.get("dtype") != "float32":
        raise SystemExit(f"derived Tahoe cache dtype {sidecar.get('dtype')!r}, expected float32")
    if int(sidecar.get("bytes", -1)) != cache.stat().st_size:
        raise SystemExit("derived Tahoe cache size does not match its sidecar; rebuild it")
    return {
        "path": str(cache),
        "source_sha256": sidecar.get("source_sha256"),
        "layer_sha256": sidecar.get("layer_sha256"),
        "file_sha256": sidecar.get("file_sha256"),
        "shape": [int(x) for x in shape],
        "dtype": sidecar.get("dtype"),
        "bytes": int(sidecar.get("bytes")),
        "sample_max_abs_diff": sidecar.get("sample_max_abs_diff"),
        "schema": sidecar.get("schema"),
        "built_utc": sidecar.get("created_utc"),
    }


def table_cache_identity(table_cache: Path, identity: Mapping[str, Any]) -> dict[str, Any]:
    """Verify the derived CSV cache of the frozen parquet tables (pinned env has no engine)."""
    table_cache = Path(table_cache)
    manifest_path = table_cache / "manifest.json"
    if not manifest_path.exists():
        raise SystemExit(
            f"derived frozen-table CSV cache is missing: {manifest_path}\n"
            "build it once with the builder environment:\n"
            "  /home/chang/miniconda3/envs/arc/bin/python "
            f"{PIPELINE_EXECUTION}/tools/convert_frozen_tables.py "
            f"--data-dir {identity['data_dir']} --cache-dir {table_cache}"
        )
    manifest = record.read_json(manifest_path)
    if manifest.get("schema") != "ptlib.table_cache/v1":
        raise SystemExit(f"table cache schema {manifest.get('schema')!r}, expected ptlib.table_cache/v1")
    tables = manifest.get("tables") or {}
    needed = "tahoe_control_groups.parquet"
    if needed not in tables:
        raise SystemExit(f"table cache has no entry for {needed}")
    frozen_sha = (identity.get("artifacts") or {}).get(needed)
    entry = tables[needed]
    if entry.get("source_sha256") != frozen_sha:
        raise SystemExit(
            f"table cache {needed} source_sha256 {entry.get('source_sha256')} "
            f"!= frozen {frozen_sha}; rebuild with convert_frozen_tables.py"
        )
    csv_path = table_cache / str(entry.get("csv"))
    if not csv_path.exists():
        raise SystemExit(f"table cache CSV is missing: {csv_path}")
    return {
        "path": str(table_cache),
        "manifest_sha256": record.sha256_file(manifest_path),
        "tables": {
            name: {"csv": value.get("csv"), "source_sha256": value.get("source_sha256")}
            for name, value in tables.items()
        },
    }


# --------------------------------------------------------------------------- #
# command construction
# --------------------------------------------------------------------------- #
def unit_command(
    *,
    python: str,
    data_dir: Path,
    run_dir: Path,
    seed: int,
    budget: str,
    device: str,
    cache: Path,
    table_cache: Path,
    verify_cache: str,
    hash_cache: Path | None,
) -> list[str]:
    command = [
        python,
        str(PIPELINE_EXECUTION / "train.py"),
        "--stage", "pretrain",
        "--data-dir", str(data_dir),
        "--out", str(run_dir),
        "--seed", str(int(seed)),
        "--device", device,
        "--cache-mode", "require",
        "--tahoe-cache", str(cache),
        "--table-cache", str(table_cache),
        "--verify-cache", verify_cache,
    ]
    if hash_cache is not None:
        command += ["--hash-cache", str(hash_cache)]
    if budget == "smoke":
        command.append("--smoke")
    elif budget == "one-epoch":
        command += ["--epochs", "1", "--min-epochs", "1"]
    elif budget != "formal":
        raise SystemExit(f"unknown budget {budget!r}")
    return command


def _read_log_tail(path: Path, lines: int = 80) -> str:
    if not path.exists():
        return ""
    try:
        return "\n".join(path.read_text(encoding="utf-8", errors="replace").splitlines()[-lines:])
    except OSError:
        return ""


def _resources_path(run_dir: Path) -> Path:
    return run_dir / "06.5_unit_resources.json"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the ticket-#5 seeded Tahoe pretrains through the 06.4 pipeline: one "
            "--stage pretrain invocation per seed, locked config, no recipe overrides."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mode", choices=("profile", "formal"), required=True,
                        help="profile = short real-data check; formal = the three frozen seeds")
    parser.add_argument("--seeds", default=None,
                        help="comma-separated seeds (default: the frozen seeds for formal, "
                             "the first frozen seed for profile)")
    parser.add_argument("--profile-budget", choices=("smoke", "one-epoch", "all"), default="all",
                        help="profile mode only: 2-epoch bounded smoke, one full-width epoch, or both")
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT,
                        help="substep output root (the record and run trees land under <out-root>/<mode>/)")
    parser.add_argument("--data-dir", type=Path, default=FROZEN_INPUT,
                        help="frozen 06.3 input directory")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE,
                        help="derived Tahoe log1p memmap (.npy); its sidecar must bind the frozen input")
    parser.add_argument("--table-cache", type=Path, default=DEFAULT_TABLE_CACHE,
                        help="derived CSV cache of the frozen parquet tables")
    parser.add_argument("--hash-cache", type=Path, default=None,
                        help="optional 06.4 artifact-hash memo (default: fresh verification per run)")
    parser.add_argument("--python", default=sys.executable,
                        help="pinned GPU python that runs 06.4 train.py")
    parser.add_argument("--device", default="cuda", choices=("cuda", "cuda:0"),
                        help="torch device (cuda is a hard failure when unavailable)")
    parser.add_argument("--threads", type=int, default=2,
                        help="OMP/MKL/OPENBLAS threads for the trainer (2 while another job may run)")
    parser.add_argument("--interval", type=float, default=1.0,
                        help="resource sampler interval, seconds")
    parser.add_argument("--resume", action="store_true",
                        help="reuse a run directory whose manifest matches the planned config")
    parser.add_argument("--allow-concurrent", action="store_true",
                        help="do not refuse when another --stage pretrain is running")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan (commands, budgets, config digests) and exit")
    return parser


def _parse_seeds(raw: str | None, default: Sequence[int]) -> list[int]:
    if raw is None:
        return [int(seed) for seed in default]
    seeds = [int(piece) for piece in raw.replace(" ", "").split(",") if piece]
    if not seeds:
        raise SystemExit("--seeds is empty")
    return seeds


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = record.find_repo_root(EXECUTION)
    pipeline = plan.pipeline_paths()
    if pipeline["root"].resolve() != root.resolve():
        raise SystemExit(
            f"06.4 ptlib resolved a different repository root ({pipeline['root']}) than this "
            f"wrapper ({root}); the substep is not inside its pipeline"
        )
    if not (PIPELINE_EXECUTION / "train.py").exists():
        raise SystemExit(f"06.4 pipeline train.py is missing: {PIPELINE_EXECUTION / 'train.py'}")
    python = str(args.python)
    if not Path(python).exists():
        raise SystemExit(f"python interpreter not found: {python}")

    data_dir = Path(args.data_dir).resolve()
    cache = Path(args.cache).resolve()
    table_cache = Path(args.table_cache).resolve()
    out_root = Path(args.out_root).resolve()
    base = out_root / args.mode
    base.mkdir(parents=True, exist_ok=True)

    # ---- frozen contract ------------------------------------------------------
    raw_locked, frozen_config, locked_notes = plan.load_locked(data_dir)
    reduced = record.normalize_locked_config(raw_locked)
    contract_errors = record.ticket_expectation_errors(reduced)
    if contract_errors:
        raise SystemExit(
            "REFUSING: the frozen locked config does not match the ticket's locked recipe:\n  "
            + "\n  ".join(contract_errors)
        )
    identity = record.frozen_input_identity(data_dir)
    cache_rec = cache_identity(cache, identity)
    table_rec = table_cache_identity(table_cache, identity)

    seeds_planned = [int(seed) for seed in frozen_config.get("seeds") or []]
    if not seeds_planned:
        raise SystemExit("the frozen locked config declares no seeds")
    if args.mode == "formal":
        seeds = _parse_seeds(args.seeds, seeds_planned)
        budgets = ["formal"]
    else:
        seeds = _parse_seeds(args.seeds, [seeds_planned[0]])
        budgets = ["smoke", "one-epoch"] if args.profile_budget == "all" else [args.profile_budget]

    planned = plan.plan_units(frozen_config, seeds, budgets, data_dir=data_dir, device=args.device)

    # ---- unit layout ---------------------------------------------------------
    units: list[dict[str, Any]] = []
    for unit in planned:
        if args.mode == "formal":
            unit_id = f"tahoe_seed{unit['seed']}"
        else:
            unit_id = f"tahoe_seed{unit['seed']}_{unit['budget']}"
        run_dir = base / "runs" / unit_id
        command = unit_command(
            python=python,
            data_dir=data_dir,
            run_dir=run_dir,
            seed=unit["seed"],
            budget=unit["budget"],
            device=args.device,
            cache=cache,
            table_cache=table_cache,
            verify_cache="deep" if unit["budget"] == "formal" else "metadata",
            hash_cache=Path(args.hash_cache).resolve() if args.hash_cache else None,
        )
        units.append({**unit, "unit_id": unit_id, "run_dir": run_dir, "command": command})
    expectations = plan.expectation_block(
        frozen_config, identity["locked_config_sha256"], units
    )

    host = host_block(args.device)
    memory_now = round(available_gib(), 2)
    if memory_now < MIN_AVAILABLE_GIB:
        raise SystemExit(
            f"REFUSING: only {memory_now:.1f} GiB RAM available (< {MIN_AVAILABLE_GIB:.0f} GiB floor); "
            "free memory first (docs/agents/execution-environment.md)"
        )
    if args.mode == "formal" and memory_now < FORMAL_AVAILABLE_GIB:
        print(
            f"WARNING: {memory_now:.1f} GiB available; a formal pretrain should start with "
            f">= {FORMAL_AVAILABLE_GIB:.0f} GiB free. Continuing because the cap is enforced.",
            file=sys.stderr,
        )
    busy = other_pretrain_running()
    if busy and not args.allow_concurrent:
        raise SystemExit(
            "REFUSING: another `--stage pretrain` process is running (one training job at a time):\n  "
            + "\n  ".join(busy)
            + "\nPass --allow-concurrent only if you have re-read the execution-environment rules."
        )

    config_overrides = {str(unit["unit_id"]): dict(unit["overrides"]) for unit in units}
    if args.mode == "formal":
        record_budget = {
            "kind": "locked",
            "smoke": False,
            "row_subset": False,
            "profile_budgets": None,
            "overrides": {},
            "overrides_by_unit": {},
            **{key: expectations["budget"][key] for key in expectations["budget"]},
        }
    else:
        record_budget = {
            "kind": "profile",
            "smoke": all(unit["budget"] == "smoke" for unit in units),
            "row_subset": any(unit["budget"] in ("smoke", "one-epoch") for unit in units),
            "profile_budgets": list(budgets),
            "overrides": {},
            "overrides_by_unit": config_overrides,
            "reference": dict(expectations["budget"]),
        }

    plan_payload = {
        "mode": args.mode,
        "repo_root": str(root),
        "data_dir": str(data_dir),
        "cache": cache_rec,
        "table_cache": table_rec,
        "frozen_seeds": seeds_planned,
        "seeds_requested": seeds,
        "memory_available_gib": memory_now,
        "budget": record_budget,
        "expected_config_sha256": expectations["config_sha256"],
        "units": [
            {
                "unit_id": unit["unit_id"],
                "seed": unit["seed"],
                "budget": unit["budget"],
                "run_dir": str(unit["run_dir"]),
                "config_sha256": unit["config_sha256"],
                "budget_config": unit["budget_config"],
                "overrides": unit["overrides"],
                "command": unit["command"],
            }
            for unit in units
        ],
    }
    if args.dry_run:
        print(json.dumps(plan_payload, indent=2, default=str))
        return 0

    # ---- run ---------------------------------------------------------------
    print(f"06.5 Tahoe pretrain — mode={args.mode} seeds={seeds} budgets={budgets}", flush=True)
    print(f"  frozen input: {data_dir} ({len(identity['artifacts'])} artifacts hashed)", flush=True)
    print(f"  cache: {cache} (source sha256 {cache_rec['source_sha256'][:16]}…)", flush=True)
    print(f"  out:   {base}", flush=True)
    print(f"  host:  {memory_now:.1f} GiB available, cgroup {host['cgroup'].get('memory_max')}",
          flush=True)
    for unit in units:
        print(f"  plan {unit['unit_id']}: {unit['budget']} cfg {unit['config_sha256'][:16]}…",
              flush=True)

    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = str(max(1, int(args.threads)))

    entries: list[dict[str, Any]] = []
    failed = False
    for index, unit in enumerate(units, start=1):
        run_dir: Path = unit["run_dir"]
        print(
            f"\n[{index}/{len(units)}] {unit['unit_id']} (seed {unit['seed']}, {unit['budget']})",
            flush=True,
        )
        print("  " + " ".join(unit["command"]), flush=True)
        resources: dict[str, Any] = {}
        outer_wall: float | None = None
        exit_code = 1
        manifest_file = run_dir / "manifest.json"
        if manifest_file.exists():
            try:
                existing = record.read_json(manifest_file)
            except record.RecordError:
                existing = {}
            if existing.get("config_sha256") == unit["config_sha256"]:
                if not args.resume:
                    raise SystemExit(
                        f"REFUSING to overwrite a completed run at {run_dir}; "
                        "pass --resume to reuse it, or move it away"
                    )
                resources = _read_resources(run_dir)
                exit_code = 0
                print("  resumed: an existing run matches this config digest", flush=True)
            else:
                raise SystemExit(
                    f"REFUSING: {run_dir} holds a run with a different config digest "
                    f"({existing.get('config_sha256')}); move it away before rerunning"
                )
        if exit_code != 0:
            run_dir.mkdir(parents=True, exist_ok=True)
            log_path = run_dir / "train_stdout.log"
            started = time.perf_counter()
            with open(log_path, "w", encoding="utf-8") as handle:
                handle.write("$ " + " ".join(unit["command"]) + "\n")
                handle.flush()
                process = subprocess.Popen(
                    unit["command"], stdout=handle, stderr=subprocess.STDOUT, env=env, cwd=root
                )
                sampler = ResourceSampler(process.pid, interval_s=args.interval)
                sampler.start()
                exit_code = process.wait()
                outer_wall = time.perf_counter() - started
                resources = sampler.stop()
                resources["host_peak_rss_mb_pipeline"] = None  # from the trainer manifest below
                res_path = _resources_path(run_dir)
                record.write_json(res_path, resources)
            print(f"  exit {exit_code} in {outer_wall:.1f}s", flush=True)
        entry = record.run_entry(
            run_dir=run_dir,
            unit_id=unit["unit_id"],
            seed=unit["seed"],
            budget=unit["budget"],
            command=unit["command"],
            exit_code=exit_code,
            expected_config_sha256=unit["config_sha256"],
            wall_s_outer=outer_wall,
            resources=resources,
            log_tail=_read_log_tail(run_dir / "train_stdout.log") if exit_code != 0 else None,
        )
        manifest_peak = (entry.get("resources") or {}).get("host_peak_rss_mb_pipeline")
        if entry["status"] == "completed":
            card = run_dir / "model_card.md"
            if card.exists():
                target = base / "model_cards" / unit["unit_id"] / "model_card.md"
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(card, target)
                entry["model_card"]["tracked_copy"] = str(target)
                entry["model_card"]["tracked_copy_rel"] = _rel(target, root)
            auprc = entry.get("best_valid_auprc")
            print(
                f"  OK {unit['unit_id']}: epochs={entry.get('epochs_run')} "
                f"stop={entry.get('stop_reason')} best_valid_AUPRC@200={auprc} "
                f"wall={entry['timings'].get('wall_s_pipeline')}s "
                f"host_peak={manifest_peak} MiB "
                f"gpu_peak={resources.get('gpu_peak_mib_attributed')} MiB "
                f"(device {resources.get('gpu_peak_mib_device')})",
                flush=True,
            )
        else:
            failed = True
            print(f"  FAILED {unit['unit_id']}: {entry.get('failures')}", flush=True)
            if entry.get("log_tail"):
                print("  --- last log lines ---", flush=True)
                print(entry["log_tail"], flush=True)
        entries.append(entry)
        if exit_code in (-9, 137):
            print(
                "  STOPPING the sweep: the trainer was killed (exit -9/137, likely the "
                "cgroup memory cap or the OOM killer). Fix the cause before rerunning.",
                flush=True,
            )
            break

    seeds_completed = sorted({int(entry["seed"]) for entry in entries
                              if entry.get("status") == "completed"})
    seeds_missing = [seed for seed in seeds_planned if seed not in seeds_completed]
    complete = not failed and bool(seeds_completed) and (args.mode != "formal" or not seeds_missing)
    rec = record.build_record(
        mode=args.mode,
        units=entries,
        budget=record_budget,
        expectations=expectations,
        input_identity={
            **identity,
            "cache": cache_rec,
            "table_cache": table_rec,
        },
        command=sys.argv,
        git=record.git_info(root),
        seeds_planned=seeds_planned if args.mode == "formal" else sorted({int(s) for s in seeds}),
        seeds_requested=seeds,
        notes=[
            "formal runs use the frozen locked recipe with no overrides; the only budget "
            "overrides in this record are the documented profile budgets"
            if args.mode == "profile"
            else "formal runs use the frozen locked recipe and no budget override",
            f"locked config notes: mapping={locked_notes.get('mapping')} "
            f"defaulted={locked_notes.get('defaulted')}",
        ],
        host=host,
    )
    rec["complete"] = complete
    stem = (
        "profile_record"
        if args.mode == "profile"
        else ("tahoe_pretrain_record" if complete else "tahoe_pretrain_record_partial")
    )
    paths = record.write_record(rec, base, stem=stem)
    print(f"\nrecord: {paths['json']}", flush=True)
    print(f"curve:  {paths['csv']}", flush=True)

    sweep_is_complete_by_design = args.mode == "formal" and set(seeds) == set(seeds_planned)
    errors, warnings = record.check_record(
        rec,
        locked=reduced,
        expect_mode=args.mode,
        expect_complete=sweep_is_complete_by_design,
    )
    for warning in warnings:
        print(f"check warning: {warning}", flush=True)
    if errors:
        print(f"record check FAILED ({len(errors)} error(s)):", flush=True)
        for error in errors:
            print(f"  - {error}", flush=True)
        return 1
    print(f"record check: PASS ({len(rec['units'])} unit(s), mode {args.mode})", flush=True)
    if args.mode == "formal" and not complete:
        print(
            "PARTIAL formal sweep: the record is written as "
            f"{paths['json'].name}; run the remaining seeds before treating it as evidence.",
            flush=True,
        )
        return 0
    return 0


def _read_resources(run_dir: Path) -> dict[str, Any]:
    path = _resources_path(run_dir)
    if not path.exists():
        return {}
    try:
        return record.read_json(path)
    except record.RecordError:
        return {}


def _rel(path: Path, root: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(root))
    except ValueError:
        return str(path)


if __name__ == "__main__":
    raise SystemExit(main())
