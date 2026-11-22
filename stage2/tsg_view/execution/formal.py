#!/usr/bin/env python3
"""Drive the sealed 45-fit `tsg_view` matrix: plan, run, verify.

One bounded GPU scope, one fit at a time, one fresh process per phase. The driver:

1. preflights (available RAM >= 12 GiB, exactly one usable GPU) and seals the plan --
   schedule, freeze/truth/code digests, device, caps -- before the first fit;
2. groups the schedule by ``(source, fold)`` and, for each group, completes every
   validation (train -> fresh-process replay) before sealing the fold's validation choice
   receipt, then runs every test -> fresh-process replay;
3. records 10-second RAM/GPU telemetry per phase, a per-run record (seed, realized
   epochs, best epoch, artifact hashes, device, runtime, peak RSS/GPU) and a resumable
   state machine keyed on the run's ``terminal.json``;
4. ``verify`` re-checks the sealed plan, every run manifest, the fold receipts, the
   separate-process tokens, and reports PASS/DRIFT for each ``full``-arm origin against
   the sealed no-dose cold-drug-dose anchors.

Use ``plan`` to write/print the plan only, ``run`` to execute (``--only-source`` and
``--only-fold`` stage it; a filter that would split a fold group is refused because the
choice receipt is per fold), and ``verify`` to check a completed or partial candidate.

Every phase runs inside the one capped scope: the driver re-execs itself under
``systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G`` (``--no-scope`` is a
test hook only).
"""
import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import learner  # noqa: E402
import matrix  # noqa: E402
from views import repo_root, sha256, write_json  # noqa: E402

RUN_REPORT_SCHEMA = "tsg-view-formal-run-v1"
VERIFY_SCHEMA = "tsg-view-formal-verify-v1"
CAP_SCHEMA = "tsg-view-formal-cap-v1"
DEFAULT_MIN_RAM_KIB = 12 * 1024 * 1024
MIN_GPU_FREE_MIB = 4096
TELEMETRY_SECONDS = 10.0
MEMINFO_KEYS = {"MemTotal": "mem_total_kib", "MemAvailable": "mem_available_kib",
                "SwapFree": "swap_free_kib"}
RUN_REPORT_FILE = "formal_run_report.json"
#: Orchestration-only code. A driver bug fix must be able to resume an otherwise
#: complete candidate without rewriting its sealed plan; the drift is recorded in the
#: run/verify report, never hidden, and any scientific/freeze/truth/runtime change
#: still fails closed.
DRIVER_CODE_FILES = ("formal.py",)


#: The plan's two code-digest blocks: the driver/verifier files (``code``) and the worker
#: plus reused shared modules the runs were fitted with (``learner_runtime``).
CODE_SECTIONS = ("code", "learner_runtime")


def code_drift(sealed, current):
    """``(entries, inputs_identical)`` for every formal code digest that moved.

    ``inputs_identical`` is ``True`` only when the two plans agree on everything except
    the two code blocks, i.e. the plan body, freezes, truth, native source, runtime
    environment and schedule are unchanged. That distinction keeps the fail-closed rule
    about *inputs* while letting a verified refactor be reported instead of refused.
    """
    sealed_copy = json.loads(json.dumps(sealed))
    current_copy = json.loads(json.dumps(current))
    entries = []
    for section in CODE_SECTIONS:
        before_block = sealed_copy.pop(section, None) or {}
        after_block = current_copy.pop(section, None) or {}
        for name in sorted(set(before_block) | set(after_block)):
            before, after = before_block.get(name), after_block.get(name)
            if before != after:
                entries.append({"section": section, "file": name,
                                "sealed_sha256": before, "current_sha256": after})
    return entries, sealed_copy == current_copy


def driver_code_drift(sealed, current):
    """Return the driver-only code drift, ``[]`` if identical, or ``None`` (fail closed)."""
    if sealed == current:
        return []
    entries, identical = code_drift(sealed, current)
    if not identical:
        return None
    driver = [entry for entry in entries
              if entry["section"] == "code" and entry["file"] in DRIVER_CODE_FILES]
    if len(driver) != len(entries):
        return None
    return driver


def code_identity(sealed, current):
    """The transparency block for a verify pass, or ``None`` when an input moved.

    ``verify`` must be able to run on newer code than the one that produced the runs --
    otherwise no refactor could ever be verified -- so it records every code file that
    moved (``code_drift``) instead of refusing. A change to anything the runs were
    fitted against still fails closed. ``run``/resume keeps the stricter driver-only rule.
    """
    entries, identical = code_drift(sealed, current)
    if not identical:
        return None
    return {"sealed": sealed.get("code"), "current": current.get("code"),
            "sealed_runtime": sealed.get("learner_runtime"),
            "current_runtime": current.get("learner_runtime"),
            "status": "MATCH" if not entries else "CODE_DRIFT",
            "code_drift": entries,
            "driver_code_drift": [entry for entry in entries
                                  if entry["section"] == "code"
                                  and entry["file"] in DRIVER_CODE_FILES],
            "note": ("verify records every code digest that moved since the plan was "
                     "sealed (both the driver/verifier files and the learner runtime); the "
                     "refactor-only drift is reported here and never hidden. Any non-code "
                     "input change still fails closed, a replay still refuses a moved "
                     "learner runtime, and run/resume keeps the driver-file-only allowance")}


def emit_json(value, *, pretty=False):
    """Best-effort stdout emission.

    A blocked or broken stdout pipe (``BlockingIOError``/``OSError``/``BrokenPipeError``)
    must never decide the driver's exit status: the durable report files are the
    contract. Returns whether the text was written.
    """
    text = json.dumps(value, indent=2 if pretty else None, sort_keys=True) + "\n"
    try:
        sys.stdout.write(text)
        sys.stdout.flush()
    except (BlockingIOError, OSError):
        return False
    return True


class _NullStdout:
    """A stdout stand-in used only when the real stream cannot be flushed."""

    def write(self, *args, **kwargs):
        return 0

    def flush(self):
        return None

    def isatty(self):
        return False


def quiet_stdout():
    """Make interpreter-shutdown's stdout flush harmless after best-effort emission.

    CPython exits 120 when flushing ``sys.stdout`` fails at shutdown, so a blocked pipe
    could still change the exit status after ``main`` returned 0. Once writes have
    failed we either flush cleanly or replace the stream with a no-op.
    """
    try:
        sys.stdout.flush()
    except (BlockingIOError, OSError):
        sys.stdout = _NullStdout()


def write_run_report(candidate, report):
    """Write the durable final run report; this file, not stdout, is the evidence."""
    path = Path(candidate) / RUN_REPORT_FILE
    write_json(path, report)
    return path
TIME_BINARY = Path("/usr/bin/time")
NVIDIA = shutil.which("nvidia-smi") or "/usr/lib/wsl/lib/nvidia-smi"


# --------------------------------------------------------------------------- #
# resource evidence
# --------------------------------------------------------------------------- #
def parse_meminfo(text):
    """Parse the three sampled /proc/meminfo fields (kB); refuse a partial sample."""
    values = {}
    for line in text.splitlines():
        name, _, rest = line.partition(":")
        if name in MEMINFO_KEYS and rest.strip():
            try:
                values[MEMINFO_KEYS[name]] = int(rest.split()[0])
            except ValueError as error:
                raise ValueError(f"malformed /proc/meminfo line: {line!r}") from error
    missing = set(MEMINFO_KEYS.values()) - set(values)
    if missing:
        raise ValueError(f"/proc/meminfo sample misses {sorted(missing)}")
    return values


def parse_gpu_line(text):
    """Parse one ``nvidia-smi --query-gpu=memory.used,memory.free,utilization.gpu`` row."""
    fields = [field.strip() for field in str(text).strip().splitlines()[0].split(",")]
    if len(fields) != 3:
        raise ValueError(f"malformed nvidia-smi sample: {text!r}")
    try:
        used, free, util = (int(fields[0]), int(fields[1]), int(fields[2]))
    except ValueError as error:
        raise ValueError(f"malformed nvidia-smi sample: {text!r}") from error
    return {"gpu_used_mib": used, "gpu_free_mib": free, "gpu_util_pct": util}


def sample_resources(*, meminfo_path=Path("/proc/meminfo"), nvidia=NVIDIA,
                     runner=subprocess.run):
    sample = {"utc_unix_seconds": round(time.time(), 3),
              **parse_meminfo(Path(meminfo_path).read_text())}
    try:
        completed = runner([nvidia, "--query-gpu=memory.used,memory.free,utilization.gpu",
                            "--format=csv,noheader,nounits"], capture_output=True, text=True,
                           timeout=15, check=True)
        sample.update(parse_gpu_line(completed.stdout.splitlines()[0]))
    except (OSError, subprocess.SubprocessError, ValueError, IndexError):
        sample.update({"gpu_used_mib": None, "gpu_free_mib": None, "gpu_util_pct": None})
    return sample


def summarize_telemetry(samples):
    """Peak/min reductions; missing GPU fields are None, never silently zero."""
    if not samples:
        return {"samples": 0, "min_mem_available_kib": None, "peak_gpu_used_mib": None,
                "max_gpu_util_pct": None, "mean_gpu_util_pct": None}
    available = [sample["mem_available_kib"] for sample in samples]
    gpu_used = [sample["gpu_used_mib"] for sample in samples
                if sample.get("gpu_used_mib") is not None]
    gpu_util = [sample["gpu_util_pct"] for sample in samples
                if sample.get("gpu_util_pct") is not None]
    return {
        "samples": len(samples),
        "min_mem_available_kib": min(available),
        "max_mem_available_kib": max(available),
        "peak_gpu_used_mib": max(gpu_used) if gpu_used else None,
        "max_gpu_util_pct": max(gpu_util) if gpu_util else None,
        "mean_gpu_util_pct": round(sum(gpu_util) / len(gpu_util), 2) if gpu_util else None,
    }


class Telemetry:
    """10-second /proc/meminfo + nvidia-smi sampler writing one JSONL block per phase.

    A sample is written at the start and at the stop of every phase in addition to the
    10-second cadence, so even a short bounded phase carries a non-empty resource trace.
    """

    def __init__(self, path, interval=TELEMETRY_SECONDS, sampler=sample_resources):
        self.path = Path(path)
        self.interval = float(interval)
        self._sampler = sampler
        self.samples = []
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._thread = None
        self._handle = None

    def _sample(self):
        try:
            sample = self._sampler()
        except (OSError, ValueError):
            return
        with self._lock:
            self.samples.append(sample)
            if self._handle is not None:
                self._handle.write(json.dumps(sample, sort_keys=True) + "\n")
                self._handle.flush()

    def _run(self):
        while not self._stop.wait(self.interval):
            self._sample()

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("w")
        self._sample()
        self._thread = threading.Thread(target=self._run, name="resource-telemetry",
                                        daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=20)
        self._sample()
        with self._lock:
            self._handle.close()
            self._handle = None
        return False

    @property
    def summary(self):
        return summarize_telemetry(self.samples)


# --------------------------------------------------------------------------- #
# preflight / cap
# --------------------------------------------------------------------------- #
def memory_preflight(min_ram_kib=DEFAULT_MIN_RAM_KIB, meminfo_path=Path("/proc/meminfo")):
    available = parse_meminfo(Path(meminfo_path).read_text())["mem_available_kib"]
    if available < int(min_ram_kib):
        raise RuntimeError(
            f"available memory {available} KiB is below {int(min_ram_kib)} KiB; stop "
            "before the next formal phase")
    return available


def gpu_preflight(device, nvidia=NVIDIA, runner=subprocess.run):
    if device != "cuda":
        return None
    try:
        completed = runner([nvidia, "--query-gpu=name,memory.total,memory.free",
                            "--format=csv,noheader,nounits"], capture_output=True, text=True,
                           timeout=20, check=True)
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError(f"no usable NVIDIA GPU: {error}") from error
    rows = [row.split(",") for row in completed.stdout.strip().splitlines() if row.strip()]
    if not rows:
        raise RuntimeError("nvidia-smi returned no GPU")
    name, total, free = (field.strip() for field in rows[0][:3])
    if int(free) < MIN_GPU_FREE_MIB:
        raise RuntimeError(f"GPU free memory {free} MiB is below {MIN_GPU_FREE_MIB} MiB")
    return {"count": len(rows), "name": name, "total_mib": int(total), "free_mib": int(free)}


def preflight(device="cuda", min_ram_kib=DEFAULT_MIN_RAM_KIB, **arguments):
    return {"available_ram_kib": memory_preflight(min_ram_kib, **{
        key: value for key, value in arguments.items() if key == "meminfo_path"}),
        "gpu": gpu_preflight(device, **{key: value for key, value in arguments.items()
                                        if key in ("nvidia", "runner")})}


def scope_command(argv, unit):
    return ["systemd-run", "--user", "--scope", "--unit", unit,
            "-p", "MemoryMax=10G", "-p", "MemorySwapMax=8G", "--", *map(str, argv)]


def record_cap(candidate, unit):
    """Record the systemd scope identity the run is executing inside."""
    evidence = {"schema": CAP_SCHEMA, "unit": unit,
                "memory_max": matrix.CAPS["memory_max"],
                "memory_swap_max": matrix.CAPS["memory_swap_max"],
                "threads": matrix.CAPS["threads"]}
    try:
        completed = subprocess.run(
            ["systemctl", "--user", "show", f"{unit}.scope", "-p", "Id", "-p",
             "InvocationID", "-p", "MemoryMax", "-p", "MemorySwapMax"],
            capture_output=True, text=True, timeout=20, check=False)
        if completed.returncode == 0:
            evidence["scope"] = dict(
                line.split("=", 1) for line in completed.stdout.splitlines() if "=" in line)
    except (OSError, subprocess.SubprocessError):
        pass
    write_json(Path(candidate) / "formal_cap.json", evidence)
    return evidence


# --------------------------------------------------------------------------- #
# phase invocation
# --------------------------------------------------------------------------- #
def _child_environment():
    environment = dict(os.environ)
    environment.update({"CUBLAS_WORKSPACE_CONFIG": ":4096:8",
                        "OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2",
                        "OPENBLAS_NUM_THREADS": "2", "NUMEXPR_NUM_THREADS": "2",
                        "PYTHONUNBUFFERED": "1"})
    return environment


def phase_command(context, origin, action, selection):
    learner_phase = {"train": "train", "replay": "replay", "evaluate": "test"}[action]
    run = matrix.run_path(context["candidate"], origin)
    freeze = matrix.freeze_path(context["freezes"], origin["view"], origin["source"])
    command = [sys.executable, str(HERE / "learner.py"), str(freeze),
               str(context["native_source"]), str(run), learner_phase,
               "--device", context["device"], "--fold", str(int(origin["fold"])),
               "--seed", str(int(origin["seed"])), "--mode", context["mode"]]
    if context.get("repo_root_path") is not None:
        command += ["--repo-root", str(context["repo_root_path"])]
    if context.get("sealed") is not None:
        command += ["--sealed", str(context["sealed"])]
    if context["mode"] != "smoke":
        command += ["--matrix-root", str(context["candidate"]),
                    "--truth", str(context["truth"])]
        if selection is not None:
            command += ["--selection", str(selection)]
        if context["mode"] == "probe" and action == "train":
            command += ["--rows", str(context["probe"]["rows"]),
                        "--epochs", str(context["probe"]["epochs"])]
    return command


def _time_parser(path):
    if not Path(path).is_file():
        return {}
    values = {}
    for line in Path(path).read_text().splitlines():
        name, _, rest = line.partition(":")
        if name.strip() == "Maximum resident set size (kbytes)":
            values["peak_rss_kib"] = int(rest.strip())
        elif name.strip() == "Exit status":
            values["time_exit_status"] = int(rest.strip())
    return values


def invoke(context, origin, action, phase, *, selection, min_ram_kib,
           telemetry_seconds, runner=subprocess.run):
    """Run exactly one fresh process for one phase, under telemetry and the RAM floor."""
    available = memory_preflight(min_ram_kib)
    candidate = Path(context["candidate"])
    attempt = matrix.resources_path(candidate, origin, phase, action, time.time_ns())
    attempt.mkdir(parents=True)
    plan = matrix.load_plan(candidate)
    command = phase_command(context, origin, action, selection)
    started = time.time()
    write_json(attempt / "phase.json", {
        "origin": origin, "phase": phase, "action": action,
        "available_ram_kib": available,
        "started_utc_unix_seconds": started,
        "plan_sha256": sha256(candidate / "formal_plan.json"),
        "memory_max": matrix.CAPS["memory_max"],
        "memory_swap_max": matrix.CAPS["memory_swap_max"],
        "threads": matrix.CAPS["threads"]})
    wrapped = list(command)
    if TIME_BINARY.is_file():
        wrapped = [str(TIME_BINARY), "-v", "-o", str(attempt / "time.txt")] + wrapped
    emit_json({"event": "start", "origin": origin, "action": action, "phase": phase})
    telemetry = Telemetry(attempt / "telemetry.jsonl", interval=telemetry_seconds)
    with telemetry, (attempt / "subprocess.log").open("w") as log:
        completed = runner(wrapped, stdout=log, stderr=subprocess.STDOUT,
                           env=_child_environment(), check=False)
    elapsed = time.time() - started
    time_values = _time_parser(attempt / "time.txt")
    result = {
        "returncode": completed.returncode,
        "elapsed_seconds": round(elapsed, 3),
        "peak_rss_kib": time_values.get("peak_rss_kib"),
        "time_exit_status": time_values.get("time_exit_status"),
        "telemetry": telemetry.summary,
        "plan_mode": plan["mode"]}
    write_json(attempt / "result.json", result)
    completed.check_returncode()
    emit_json({"event": "complete", "origin": origin, "action": action, "phase": phase,
               "elapsed_seconds": result["elapsed_seconds"]})
    return result


def attempt_records(candidate, origin):
    """All recorded phase attempts of one origin, for the per-run record."""
    relative = matrix.run_relative(origin).relative_to("runs")
    base = Path(candidate) / "resources" / relative
    records = []
    if base.is_dir():
        for attempt in sorted(base.iterdir()):
            result, phase = attempt / "result.json", attempt / "phase.json"
            if not result.is_file() or not phase.is_file():
                continue
            value = json.loads(result.read_text())
            value.update(json.loads(phase.read_text()))
            value["peak_gpu_used_mib"] = (value.get("telemetry") or {}).get(
                "peak_gpu_used_mib")
            value["attempt"] = attempt.name
            records.append(value)
    return records


# --------------------------------------------------------------------------- #
# schedule grouping
# --------------------------------------------------------------------------- #
def groups_of(origins):
    groups, order = {}, []
    for origin in origins:
        key = (origin["source"], int(origin["fold"]))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(origin)
    return [groups[key] for key in order]


def check_group_completeness(full_plan, selected):
    """A folded choice receipt only exists for a complete ``(source, fold)`` group."""
    plan_groups = {(o["source"], int(o["fold"])) for o in full_plan["origins"]}
    for source, fold in plan_groups:
        whole = {matrix.origin_key(o) for o in full_plan["origins"]
                 if (o["source"], int(o["fold"])) == (source, fold)}
        chosen = {matrix.origin_key(o) for o in selected
                  if (o["source"], int(o["fold"])) == (source, fold)}
        if chosen and chosen != whole:
            raise ValueError(
                f"filters split the {source}/fold-{fold} group; the validation choice "
                "receipt is sealed per fold, so stage with --only-source/--only-fold or "
                "select whole groups")


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
def status(run):
    path = Path(run) / "terminal.json"
    if path.is_file():
        return json.loads(path.read_text()).get("status")
    return "interrupted" if Path(run).exists() else None


def build_context(candidate, *, freezes, truth, native_source, device, mode,
                  probe=None, repo_root_path=None, sealed=None):
    return {"candidate": Path(candidate), "freezes": Path(freezes), "truth": Path(truth),
            "native_source": Path(native_source), "device": str(device), "mode": mode,
            "probe": probe, "repo_root_path": repo_root_path, "sealed": sealed}


def plan_for(*, freezes, truth, native_source, device, mode, probe_rows=None,
             probe_epochs=None, only_source=(), only_view=(), only_fold=(),
             only_seed=(), repo_root_path=None):
    return matrix.frozen_plan(
        freezes=freezes, truth=truth, source_dir=native_source, device=device, mode=mode,
        probe_rows=probe_rows, probe_epochs=probe_epochs, only_source=only_source,
        only_view=only_view, only_fold=only_fold, only_seed=only_seed,
        repo_root_path=repo_root_path)


def run_matrix(candidate, *, freezes, truth, native_source, device="cuda", mode="formal",
               probe_rows=None, probe_epochs=None, only_source=(), only_view=(),
               only_fold=(), only_seed=(), min_ram_kib=DEFAULT_MIN_RAM_KIB,
               telemetry_seconds=TELEMETRY_SECONDS, repo_root_path=None, sealed=None,
               runner=subprocess.run):
    """Execute (or resume) the selected schedule; returns the run report."""
    context = build_context(candidate, freezes=freezes, truth=truth,
                            native_source=native_source, device=device, mode=mode,
                            probe={"rows": probe_rows, "epochs": probe_epochs}
                            if mode == "probe" else None,
                            repo_root_path=repo_root_path, sealed=sealed)
    candidate = context["candidate"]
    plan = plan_for(freezes=freezes, truth=truth, native_source=native_source,
                    device=device, mode=mode, probe_rows=probe_rows,
                    probe_epochs=probe_epochs,
                    only_source=only_source if mode == "probe" else (),
                    only_view=only_view if mode == "probe" else (),
                    only_fold=only_fold if mode == "probe" else (),
                    only_seed=only_seed if mode == "probe" else (),
                    repo_root_path=repo_root_path)
    if (candidate / "formal_plan.json").is_file():
        sealed_plan = matrix.load_plan(candidate)
        code_drift = driver_code_drift(sealed_plan, plan)
        if code_drift is None:
            raise ValueError("sealed formal plan, runtime or input identity changed")
    else:
        matrix.seal_plan(candidate, plan)
        sealed_plan = matrix.load_plan(candidate)
        code_drift = []
    selected = matrix.filter_origins(plan["origins"], only_source=only_source,
                                     only_view=only_view, only_fold=only_fold,
                                     only_seed=only_seed)
    check_group_completeness(plan, selected)
    records = []
    with (candidate / ".driver.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for group in groups_of(selected):
            source, fold = group[0]["source"], int(group[0]["fold"])
            for origin in group:
                while True:
                    action = matrix.next_action(status(matrix.run_path(candidate, origin)),
                                                "validation")
                    if action is None:
                        break
                    invoke(context, origin, action, "validation", selection=None,
                           min_ram_kib=min_ram_kib,
                           telemetry_seconds=telemetry_seconds, runner=runner)
            selection = matrix.lock_selection(candidate, source, fold)
            for origin in group:
                while True:
                    action = matrix.next_action(status(matrix.run_path(candidate, origin)),
                                                "test")
                    if action is None:
                        break
                    invoke(context, origin, action, "test", selection=selection,
                           min_ram_kib=min_ram_kib,
                           telemetry_seconds=telemetry_seconds, runner=runner)
            for origin in group:
                record = matrix.run_record(candidate, origin,
                                           attempt_records(candidate, origin))
                matrix.write_record(candidate, origin, record)
                records.append(record)
            write_json(candidate / "run_records.json", {
                "schema": RUN_REPORT_SCHEMA,
                "plan_sha256": matrix.plan_sha256(candidate),
                "mode": sealed_plan["mode"],
                "selected": len(selected),
                "complete": len(records),
                "runs": sorted(records, key=lambda item: matrix.origin_key(item["origin"]))})
    return {"schema": RUN_REPORT_SCHEMA, "plan_sha256": matrix.plan_sha256(candidate),
            "mode": sealed_plan["mode"], "selected": len(selected),
            "complete": len(records), "runs": records,
            "driver_code_drift": code_drift}


# --------------------------------------------------------------------------- #
# verify
# --------------------------------------------------------------------------- #
#: The sealed Step 04 anchors are equal-fold means at this training seed only. The matrix
#: runs ``full`` at both seeds, so only the seed-42 arms may claim anchor reproduction.
ANCHOR_SEED = 42


def sensitivity_record(source, metrics, *, level, seed, anchor_seed=ANCHOR_SEED):
    """A non-anchor ``full`` arm: its distance to the anchor, with no PASS/DRIFT claim."""
    if int(seed) == int(anchor_seed):
        raise ValueError(f"the anchor seed {anchor_seed} cannot be labelled a sensitivity")
    record = matrix.compare_full_arm(source, metrics, level=level)
    record.update({
        "status": "sensitivity", "anchor_seed": int(anchor_seed), "seed": int(seed),
        "note": (f"the sealed anchor was measured at seed {int(anchor_seed)}; this "
                 f"seed-{int(seed)} arm is reported for transparency and claims no anchor "
                 "reproduction")})
    return record


def full_arm_evidence(measured):
    """Split the ``full`` arm into the sealed-anchor check and a labelled sensitivity.

    ``measured`` is a list of ``(origin, primary_metrics)``. Every outcome is recorded --
    drift is never hidden -- but only the anchor seed's origins and fold means are
    compared with the sealed anchors (PASS/DRIFT). A non-anchor seed is a sensitivity:
    it carries the anchor, the measured value and the drift, and its status is
    ``"sensitivity"`` so it can never be read as a reproduced sealed score.

    ``comparisons`` is the number of anchor-comparable records returned alongside the
    block: it is zero when the selection carries no seed-anchor ``full`` arm, which
    ``verify`` refuses instead of reporting a vacuous PASS.
    """
    origins, fold_mean, sensitivity, comparisons = [], [], [], []
    grouped = {}
    for origin, primary in measured:
        source, seed = origin["source"], int(origin["seed"])
        grouped.setdefault((source, seed), []).append(primary)
        if seed == ANCHOR_SEED:
            record = matrix.compare_full_arm(source, primary)
            record.update({"origin": origin, "anchor_seed": ANCHOR_SEED})
            origins.append(record)
            comparisons.append(record)
        else:
            record = sensitivity_record(source, primary, level="origin", seed=seed)
            record["origin"] = origin
            sensitivity.append(record)
    for (source, seed), values in sorted(grouped.items()):
        mean = {key: sum(item[key] for item in values) / len(values)
                for key in ("pcc", "rmse")}
        if seed == ANCHOR_SEED:
            record = matrix.compare_full_arm(source, mean, level="fold_mean")
            record.update({"source": source, "seed": seed, "folds": len(values),
                           "anchor_seed": ANCHOR_SEED})
            fold_mean.append(record)
            comparisons.append(record)
        else:
            record = sensitivity_record(source, mean, level="fold_mean", seed=seed)
            record.update({"source": source, "folds": len(values)})
            sensitivity.append(record)
    return ({"anchor_seed": ANCHOR_SEED, "comparisons": len(comparisons),
             "origins": origins, "fold_mean": fold_mean,
             "sensitivity": sensitivity}, comparisons)


def _verify_origin(candidate, origin, sealed_plan, selection_sha256):
    run = matrix.run_path(candidate, origin)
    manifest = learner.verify(run)
    config = json.loads((run / "config.json").read_text())
    terminal = json.loads((run / "terminal.json").read_text())
    binding = sealed_plan["freezes"].get(f"{origin['view']}/{origin['source']}")
    expected_mode = sealed_plan["mode"]
    if (terminal.get("status") != "complete" or terminal.get("mode") != expected_mode
            or terminal.get("smoke") or not terminal.get("weights_updated")):
        raise ValueError(f"run is not a complete {expected_mode} fit: {run}")
    if config.get("origin") != origin or config.get("mode") != expected_mode:
        raise ValueError(f"run config identity mismatch: {run}")
    if ((config.get("plan") or {}).get("sha256") != matrix.plan_sha256(candidate)
            or config.get("freeze", {}).get("manifest_sha256")
            != binding["manifest_sha256"]):
        raise ValueError(f"run is bound to another plan or freeze: {run}")
    if (config.get("truth") != (config.get("plan") or {}).get("truth")
            or terminal.get("truth") != config.get("truth")):
        raise ValueError(f"run truth identity mismatch: {run}")
    identity = config["identity"]
    if (identity.get("mode") != expected_mode
            or identity.get("device") != sealed_plan["device"]):
        raise ValueError(f"run device or mode mismatch: {run}")
    if expected_mode == "formal" and int(identity.get("max_epochs", -1)) != int(
            sealed_plan["max_epochs"]):
        raise ValueError(f"run epoch budget mismatch: {run}")
    if terminal.get("test_selection_sha256") != selection_sha256:
        raise ValueError(f"run used another validation choice receipt: {run}")
    validation = json.loads((run / "validation_manifest.json").read_text())
    if any(manifest["outputs"].get(name) != digest
           for name, digest in validation["outputs"].items()):
        raise ValueError(f"sealed validation artifacts changed before test: {run}")
    tokens = [terminal.get("training_process_token"),
              terminal.get("validation_replay_token"),
              terminal.get("test_evaluation_process_token"),
              terminal.get("test_replay_token")]
    if any(not token for token in tokens) or len(set(tokens)) != 4:
        raise ValueError(f"formal phases were not four separate processes: {run}")
    return {"origin": origin, "terminal": terminal, "config": config,
            "run_manifest_sha256": sha256(run / "manifest.json"),
            "selection_sha256": selection_sha256,
            "metrics": json.loads((run / "test_metrics.json").read_text())}


def verify(candidate, *, freezes, truth, native_source, device="cuda", mode="formal",
           probe_rows=None, probe_epochs=None, only_source=(), only_view=(),
           only_fold=(), only_seed=(), repo_root_path=None, record=True):
    """Verify the sealed plan, runs, receipts, tokens and the full-arm tolerance.

    A verify whose filtered selection leaves no seed-42 ``full`` origin or fold mean to
    compare with the sealed anchors fails closed: an empty comparison set would otherwise
    make ``status`` a vacuous PASS. ``full_arm.comparisons`` records the count.
    """
    candidate = Path(candidate)
    sealed_plan = matrix.load_plan(candidate)
    if mode is not None and sealed_plan["mode"] != mode:
        raise ValueError(f"the sealed plan is a {sealed_plan['mode']} plan, not {mode}")
    identity = None
    if sealed_plan["mode"] == "formal":
        current = plan_for(freezes=freezes, truth=truth, native_source=native_source,
                           device=device, mode="formal", repo_root_path=repo_root_path)
        identity = code_identity(sealed_plan, current)
        if identity is None:
            raise ValueError("sealed formal plan, runtime or input identity changed")
    selected = matrix.filter_origins(sealed_plan["origins"], only_source=only_source,
                                     only_view=only_view, only_fold=only_fold,
                                     only_seed=only_seed)
    evidence, measured = [], []
    for group in groups_of(selected):
        source, fold = group[0]["source"], int(group[0]["fold"])
        receipt = matrix.selection_path(candidate, source, fold)
        if not receipt.is_file():
            raise ValueError(f"the fold choice receipt is missing: {receipt}")
        expected = matrix.selection_evidence(candidate, source, fold)
        if json.loads(receipt.read_text()) != expected:
            raise ValueError(f"the fold choice receipt changed: {receipt}")
        choice_sha256 = sha256(receipt)
        for origin in group:
            verified = _verify_origin(candidate, origin, sealed_plan, choice_sha256)
            evidence.append(verified)
            if origin["view"] == "full":
                measured.append((origin, verified["metrics"]["primary"]))
    full_arm, comparisons = full_arm_evidence(measured)
    if not comparisons:
        # A selection that excludes the seed-42 ``full`` arm has nothing to compare, and
        # an empty comparison set would make ``status`` a vacuous PASS. The rule is the
        # formal evidence rule; a probe verify observes it too so the two modes cannot
        # disagree about what a PASS means.
        raise ValueError(
            f"the filtered verify has no anchor-comparable full origin or fold mean: the "
            f"sealed anchors are measured at seed {ANCHOR_SEED} only, so a selection that "
            f"excludes the seed-{ANCHOR_SEED} full arm cannot reproduce them")
    status = "PASS" if all(item["status"] == "PASS" for item in comparisons) else "DRIFT"
    report = {"schema": VERIFY_SCHEMA, "status": status,
              "mode": sealed_plan["mode"], "plan_sha256": matrix.plan_sha256(candidate),
              "origins": len(evidence), "full_arm": full_arm, "code_identity": identity,
              "runs": [{"origin": item["origin"],
                        "run_manifest_sha256": item["run_manifest_sha256"],
                        "selection_sha256": item["selection_sha256"]} for item in evidence]}
    if record:
        write_json(candidate / "formal_verify.json", report)
    return report


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _plan_summary(plan_value):
    counts = {}
    for origin in plan_value["origins"]:
        counts[f"{origin['source']}/{origin['view']}"] = counts.get(
            f"{origin['source']}/{origin['view']}", 0) + 1
    return {"schema": plan_value["schema"], "mode": plan_value["mode"],
            "origins": len(plan_value["origins"]), "counts": counts,
            "device": plan_value["device"], "max_epochs": plan_value["max_epochs"],
            "freezes": sorted(plan_value["freezes"]),
            "truth": sorted(plan_value["truth"]),
            "code": plan_value["code"], "caps": plan_value["caps"]}


def main(argv=None):
    try:
        return _main(argv)
    finally:
        quiet_stdout()


def _main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=("plan", "run", "verify"))
    parser.add_argument("--candidate", type=Path, required=True,
                        help="the formal candidate root (plan, runs, selection, resources)")
    parser.add_argument("--views", type=Path,
                        default=Path(repo_root()) / "stage2/tsg_view/output/views")
    parser.add_argument("--freezes", type=Path, default=HERE / "freezes")
    parser.add_argument("--truth", type=Path,
                        default=Path(repo_root()) / "stage2/tsg_view/output/truth")
    parser.add_argument("--native-source", type=Path,
                        default=Path(repo_root()) / "pert_model/TranSiGen/source")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--mode", choices=learner.MODES, default="formal")
    parser.add_argument("--probe-rows", type=int, default=64)
    parser.add_argument("--probe-epochs", type=int, default=2)
    for name in ("source", "view", "fold", "seed"):
        parser.add_argument(f"--only-{name}", action="append", default=[])
    parser.add_argument("--sealed", type=Path)
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--min-ram-kib", type=int, default=DEFAULT_MIN_RAM_KIB)
    parser.add_argument("--telemetry-seconds", type=float, default=TELEMETRY_SECONDS)
    parser.add_argument("--no-scope", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--fail-on-drift", action="store_true")
    arguments = parser.parse_args(argv)
    filters = {"only_source": tuple(arguments.only_source),
               "only_view": tuple(arguments.only_view),
               "only_fold": tuple(int(value) for value in arguments.only_fold),
               "only_seed": tuple(int(value) for value in arguments.only_seed)}
    common = dict(freezes=arguments.freezes, truth=arguments.truth,
                  native_source=arguments.native_source, device=arguments.device,
                  mode=arguments.mode, probe_rows=arguments.probe_rows,
                  probe_epochs=arguments.probe_epochs,
                  repo_root_path=arguments.repo_root, **filters)
    if arguments.action == "plan":
        plan_value = plan_for(**common)
        path = matrix.seal_plan(arguments.candidate, plan_value)
        emit_json({"plan": str(path), **_plan_summary(plan_value)}, pretty=True)
        return 0
    if arguments.action == "verify":
        report = verify(arguments.candidate, freezes=arguments.freezes,
                        truth=arguments.truth,
                        native_source=arguments.native_source, device=arguments.device,
                        mode=arguments.mode, probe_rows=arguments.probe_rows,
                        probe_epochs=arguments.probe_epochs,
                        repo_root_path=arguments.repo_root, **filters)
        emit_json(report, pretty=True)
        return 1 if arguments.fail_on_drift and report["status"] == "DRIFT" else 0
    # run
    if not arguments.no_scope and os.environ.get("TSG_VIEW_SCOPED") != "1":
        unit = f"arc-tsg-view-formal-{os.getpid()}"
        environment = dict(os.environ, TSG_VIEW_SCOPED="1", TSG_VIEW_SCOPE_UNIT=unit)
        completed = subprocess.run(
            scope_command([sys.executable, str(Path(__file__).resolve()), *argv], unit),
            env=environment)
        return completed.returncode
    if os.environ.get("TSG_VIEW_SCOPED") == "1":
        record_cap(arguments.candidate, os.environ.get(
            "TSG_VIEW_SCOPE_UNIT", "arc-tsg-view-formal"))
    preflight_value = preflight(arguments.device, arguments.min_ram_kib)
    emit_json({"preflight": preflight_value})
    report = run_matrix(arguments.candidate, freezes=arguments.freezes,
                        truth=arguments.truth, native_source=arguments.native_source,
                        device=arguments.device, mode=arguments.mode,
                        probe_rows=arguments.probe_rows,
                        probe_epochs=arguments.probe_epochs,
                        min_ram_kib=arguments.min_ram_kib,
                        telemetry_seconds=arguments.telemetry_seconds,
                        repo_root_path=arguments.repo_root, sealed=arguments.sealed,
                        **filters)
    write_run_report(arguments.candidate, report)
    emit_json(report, pretty=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
