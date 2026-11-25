"""Preregister and replay the 36 formal Stage 2 Step 05 fits."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from tempfile import TemporaryDirectory
from threading import Event, Thread
import time

from evaluate_loss import coordinates, read_run
from objective import Objective
from run_loss import AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING, PHASES, STEP, axis_run, bind

SCHEMA = "stage2-tenv-transigen-loss-matrix-v1"
GIB = 1024**3
SAMPLE_SECONDS = 10


def run_paths(root):
    return {(loss, alpha, fold): root / "runs" / loss / f"alpha-{alpha.replace('/', '-')}" / f"fold-{fold}"
            for loss, alpha, fold in coordinates()}


def digest_values(values):
    return sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def preflight(training, representations, source):
    expected = {}
    inputs = {"training": str(training), "representations": str(representations), "source": str(source)}
    with axis_run.open_inputs(training, representations) as (data, vectors):
        folds = {
            fold: bind({"axis": AXIS, "fold": fold, "loss": "mse", "alpha": "1/3",
                        "smoke": False, "device": "cuda", **inputs}, data, vectors)
            for fold in range(3)
        }
    for loss, alpha, fold in coordinates():
        expected[loss, alpha, fold] = {**folds[fold], "objective": Objective(loss, Fraction(alpha)).identity()}
    evidence = {"input_identity": digest_values(folds), "seed": axis_run.SEED, "scenario": "cold-drug",
                "axis": AXIS, "folds": [
                    {"fold": fold, "counts": {part: len(ids) for part, ids in config["conditions"].items()},
                     "ordered_conditions_sha256": {part: digest_values(ids)
                                                   for part, ids in config["conditions"].items()}}
                    for fold, config in folds.items()
                ]}
    return expected, evidence


def check_registrations(paths, expected):
    """Refuse to fit unless all ordered populations and registered inputs match preflight."""
    if paths.keys() != expected.keys():
        raise ValueError("formal registration plan is incomplete")
    nonces, verified = set(), {}
    for key, run in paths.items():
        path = run / "registration.json"
        config = axis_run.read_json(path)
        args = config["arguments"]
        nonce = config["nonce"]
        if ((args["loss"], args["alpha"], args["fold"]) != key or type(args["fold"]) is not int
                or args["smoke"] is not False or (args["axis"], args["device"]) != (AXIS, "cuda")
                or nonce in nonces
                or {field: value for field, value in config.items() if field not in {"arguments", "nonce"}}
                != expected[key]):
            raise ValueError(f"formal registration differs from preflight: {key}")
        nonces.add(nonce)
        verified[key] = axis_run.digest(path)
    return verified


def enforce_record(path, record):
    if path.exists():
        if axis_run.read_json(path) != record:
            raise ValueError(f"sealed matrix record differs: {path.name}")
    else:
        with TemporaryDirectory(prefix=f".{path.name}-", dir=path.parent) as scratch:
            staged = Path(scratch) / path.name
            axis_run.write_json(staged, record)
            staged.rename(path)


def register_coordinate(run, scratch_parent, register):
    """Make a completed preregistration visible only after its writer exits."""
    with TemporaryDirectory(prefix=".register-", dir=scratch_parent) as scratch:
        staged = Path(scratch) / "coordinate"
        register(staged)
        axis_run.read_json(staged / "registration.json")
        run.parent.mkdir(parents=True, exist_ok=True)
        staged.rename(run)


def scope_directory():
    group = Path("/proc/self/cgroup").read_text().strip().split(":", 2)[2]
    scope = Path("/sys/fs/cgroup") / Path(group).relative_to("/")
    for filename, limit in (("memory.max", 10 * GIB), ("memory.swap.max", 8 * GIB)):
        if (scope / filename).read_text().strip() != str(limit):
            raise ValueError("formal fit requires the 10 GiB RAM and 8 GiB swap systemd scope")
    return scope


def monitor(logs, scope, phase):
    meminfo = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    available = int(meminfo["MemAvailable"].split()[0])
    total, used, utilization = map(int, subprocess.check_output(
        ["nvidia-smi", "--id=0", "--query-gpu=memory.total,memory.used,utilization.gpu",
         "--format=csv,noheader,nounits"], text=True, timeout=5).strip().split(","))
    disk_free = shutil.disk_usage(logs).free
    sample = {"unix_seconds": time.time(), "phase": phase, "wsl_available_kib": available,
              "wsl_total_kib": int(meminfo["MemTotal"].split()[0]), "gpu_total_mib": total,
              "gpu_used_mib": used, "gpu_utilization_pct": utilization, "disk_free_bytes": disk_free,
              "scope_memory_bytes": int((scope / "memory.current").read_text())}
    with (logs / "resource-samples.jsonl").open("a") as stream:
        print(json.dumps(sample, sort_keys=True), file=stream)
    if available < 12 * GIB // 1024 or disk_free < GIB or total - used < 512:
        raise RuntimeError(f"resource headroom exhausted during {phase}: {sample}")


class ResourceMonitor:
    """Sample on a monotonic schedule across both child phases and file checks."""

    def __init__(self, logs, scope):
        self.logs = logs
        self.scope = scope
        self.phase = "preflight"
        self.error = None
        self.stop = Event()
        self.thread = Thread(target=self.sample_periodically)

    def sample_periodically(self):
        deadline = self.started + SAMPLE_SECONDS
        while not self.stop.wait(max(0, deadline - time.monotonic())):
            try:
                monitor(self.logs, self.scope, self.phase)
            except Exception as error:
                self.error = error
                self.stop.set()
                return
            deadline += SAMPLE_SECONDS

    def check(self):
        if self.error is not None:
            raise self.error

    def __enter__(self):
        self.started = time.monotonic()
        monitor(self.logs, self.scope, self.phase)
        self.thread.start()
        return self

    def __exit__(self, exc_type, *_):
        self.stop.set()
        self.thread.join()
        if exc_type is None:
            self.check()


def invoke(logs, resources, coordinate, phase, args):
    resources.phase = phase
    resources.check()
    log_dir = logs / coordinate
    log_dir.mkdir(parents=True, exist_ok=True)
    base = log_dir / f"{phase}-{time.time_ns()}"
    command = ["/usr/bin/time", "-v", "-o", str(base.with_suffix(".resources")),
               sys.executable, str(STEP / "execution/run_loss.py"), phase, *args]
    with base.with_suffix(".log").open("xb") as output:
        with subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, start_new_session=True) as process:
            try:
                while True:
                    try:
                        code = process.wait(timeout=SAMPLE_SECONDS)
                        resources.check()
                        break
                    except subprocess.TimeoutExpired:
                        resources.check()
            except BaseException:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=8)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                raise
    if code:
        raise RuntimeError(f"{coordinate} {phase} failed, see {base.with_suffix('.log')}")



def execute(root, logs, training, representations, source):
    if logs.is_relative_to(root) or root.is_relative_to(logs):
        raise ValueError("matrix and audit logs require separate directories")
    logs.mkdir(parents=True, exist_ok=True)
    scope = scope_directory()
    with (ResourceMonitor(logs, scope) as resources,
          axis_run.publication_lock(root) as output):
        output.mkdir(exist_ok=True)
        paths = run_paths(output)
        expected, evidence = preflight(training, representations, source)
        preflight_path = output / "preflight.json"
        enforce_record(preflight_path, {"schema": SCHEMA, **evidence})
        inputs = ["--training", str(training), "--representations", str(representations), "--source", str(source)]
        for (loss, alpha, fold), path in paths.items():
            if not (path / "registration.json").exists():
                if path.is_dir() and not any(path.iterdir()):
                    path.rmdir()
                def register(staged):
                    invoke(logs, resources, path.relative_to(output), "register",
                           ["--loss", loss, "--alpha", alpha, "--fold", str(fold),
                            "--output", str(staged), *inputs])

                register_coordinate(path, output, register)
        registrations = check_registrations(paths, expected)
        print(json.dumps({"phase": "registered", "coordinates": len(registrations)}), flush=True)
        fits = []
        for key, path in paths.items():
            relative = path.relative_to(output)
            for phase, receipt in PHASES.items():
                if not (path / receipt).is_file():
                    invoke(logs, resources, relative, phase, ["--output", str(path), *inputs])
            resources.phase = "verify"
            resources.check()
            verified_key, (_, _, manifest, training_record, _) = read_run(path, False)
            if verified_key != key or manifest["registration_sha256"] != registrations[key]:
                raise ValueError(f"sealed run differs from the formal registration: {relative}")
            fits.append({"path": relative.as_posix(), "loss": key[0], "alpha": key[1], "fold": key[2],
                         "registration_sha256": registrations[key],
                         "run_manifest_sha256": axis_run.digest(path / "manifest.json"),
                         "best_epoch": training_record["best_epoch"], "stopping_epoch": training_record["epochs"],
                         "train_wall_seconds": training_record["wall_seconds"],
                         "train_peak_rss_kib": training_record["peak_rss_kib"]})
            print(json.dumps({"phase": "complete", "completed": len(fits), **fits[-1]}), flush=True)
        resources.check()
        resources.phase = "matrix-complete"
        scope_peak = int((scope / "memory.peak").read_text())
        enforce_record(logs / f"scope-{time.time_ns()}.json", {"memory_peak_bytes": scope_peak})
        manifest_path = output / "manifest.json"
        enforce_record(manifest_path, {"schema": SCHEMA, "preflight_sha256": axis_run.digest(preflight_path),
                                       "fits": fits})
        return {"result": "PASS", "fits": len(fits), "manifest_sha256": axis_run.digest(manifest_path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    args = parser.parse_args()
    result = execute(args.output.resolve(), args.logs.resolve(), args.training.resolve(),
                     args.representations.resolve(), args.source.resolve())
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
