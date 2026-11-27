"""Run one capped CUDA coordinate at a time and seal the matrix."""

import argparse
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

from chem_runs import ORIGINAL_CODE, PHASES, add_inputs, arguments, axis_run, bind, coordinates, open_inputs, read_run
from matrix import SAMPLE_SECONDS, ResourceMonitor, digest_values, enforce_record, register_coordinate
from matrix_measured import MEMORY_CAP, SWAP_CAP, scope_directory
from matrix_sweep import signal_process_group
from util.training.provenance import historical_gate_configuration, legacy_configuration

SCHEMA = "cpst-alpha-chemical-matrix-v1"
RUNNER = Path(__file__).with_name("chem_runs.py")


def run_paths(root, roster):
    return {(arm, fold, seed): root / "runs" / arm / f"fold-{fold}" / f"seed-{seed}"
            for arm, fold, seed in roster}


def invoke(logs, resources, coordinate, phase, args):
    resources.phase = phase
    resources.check()
    directory = logs / coordinate
    directory.mkdir(parents=True, exist_ok=True)
    attempt = directory / f"{phase}-{time.time_ns()}"
    command = ["/usr/bin/time", "-v", "-o", str(attempt.with_suffix(".resources")),
               sys.executable, str(RUNNER), phase, *args]
    with attempt.with_suffix(".log").open("xb") as output:
        with subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, start_new_session=True) as process:
            try:
                while True:
                    try:
                        status = process.wait(timeout=SAMPLE_SECONDS)
                        resources.check()
                        break
                    except subprocess.TimeoutExpired:
                        resources.check()
            except BaseException:
                if process.poll() is None:
                    signal_process_group(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=8)
                    except subprocess.TimeoutExpired:
                        signal_process_group(process.pid, signal.SIGKILL)
                        process.wait()
                raise
    if status:
        raise RuntimeError(f"{coordinate}: {phase} failed, see {attempt.with_suffix('.log')}")


def preflight(training, representations, source, smoke, *, legacy=False):
    roster = coordinates(smoke)
    with open_inputs(training, representations) as (data, vectors):
        expected = {key: bind(key, data, vectors, source, smoke) for key in roster}
        if legacy:
            expected = {key: legacy_configuration(config, axis_run.REPO, ORIGINAL_CODE)
                        for key, config in expected.items()}
        identities = list(expected.values())
        if legacy:
            identities = [historical_gate_configuration(config) for config in identities]
        populations = []
        for fold in range(3):
            ids = {phase: list(data.condition_ids(data.rows(fold, phase))) for phase in axis_run.PARTITIONS}
            populations.append({"fold": fold, "counts": {phase: len(values) for phase, values in ids.items()},
                                "conditions_sha256": {phase: digest_values(values) for phase, values in ids.items()}})
        evidence = {"schema": SCHEMA, "smoke": smoke, "formal_roster": [list(key) for key in coordinates()],
                    "frozen_representations": {arm: {"width": vectors.width(arm),
                                                  "sha256": vectors.manifest["files"][f"{arm}.npy"]}
                                               for arm, _, _ in coordinates(True)},
                    "parent_count": len(vectors.parents), "condition_count": len(data.parents),
                    "populations": populations, "registered_roster": [list(key) for key in roster],
                    "expected_identity_sha256": digest_values(identities)}
    return expected, evidence


def verify_registrations(paths, expected, smoke):
    nonces, registrations = set(), {}
    for key, path in paths.items():
        record = path / "registration.json"
        config = axis_run.read_json(record)
        nonce, args = config.pop("nonce"), config.pop("arguments")
        if nonce in nonces or args != arguments(key, smoke) or config != expected[key]:
            raise ValueError(f"registration differs from frozen preflight: {key}")
        nonces.add(nonce)
        registrations[key] = axis_run.digest(record)
    return registrations


def fit_evidence(output, paths, registrations, smoke):
    fits = []
    for key, path in paths.items():
        verified_key, (_, _, replay, training, _) = read_run(path, smoke)
        if key != verified_key or replay["registration_sha256"] != registrations[key]:
            raise ValueError(f"sealed run differs from the registered coordinate: {key}")
        fits.append({"coordinate": list(key), "path": path.relative_to(output).as_posix(),
                     "registration_sha256": registrations[key], "manifest_sha256": axis_run.digest(path / "manifest.json"),
                     "initialization": training["initialization"], "best_epoch": training["best_epoch"],
                     "stopping_epoch": training["epochs"], "train_peak_rss_kib": training["peak_rss_kib"],
                     "train_cuda_peak_allocated_bytes": training["cuda_peak_allocated_bytes"]})
    return fits


def execute(root, logs, training, representations, source, smoke=False, fresh=None):
    """Register the entire planned roster before fitting the selected fresh coordinates."""
    root, logs = Path(root).resolve(), Path(logs).resolve()
    if logs.is_relative_to(root) or root.is_relative_to(logs):
        raise ValueError("matrix and audit log directories must be separate")
    roster = coordinates(smoke)
    fresh = roster if fresh is None else tuple(fresh)
    if len(fresh) != len(set(fresh)) or not set(fresh) <= set(roster) or not fresh:
        raise ValueError("fresh fits must be a unique nonempty subset of the registered roster")
    if smoke and fresh != roster:
        raise ValueError("the bounded gate must execute all seven arms")
    scope = scope_directory()
    inputs = ["--training", str(training), "--representations", str(representations), "--source", str(source)]
    logs.mkdir(parents=True, exist_ok=True)
    with axis_run.publication_lock(root) as output, ResourceMonitor(logs, scope) as resources:
        output.mkdir(exist_ok=True)
        paths = run_paths(output, fresh)
        expected, evidence = preflight(training, representations, source, smoke)
        preflight_path = output / "preflight.json"
        enforce_record(preflight_path, {**evidence, "fresh_fits": [list(key) for key in fresh]})
        for (arm, fold, seed), path in paths.items():
            if not (path / "registration.json").is_file():
                def register(staged):
                    invoke(logs, resources, path.relative_to(output), "register",
                           ["--arm", arm, "--fold", str(fold), "--seed", str(seed), "--output", str(staged), *inputs,
                            *(["--smoke"] if smoke else [])])
                register_coordinate(path, output, register)
        registrations = verify_registrations(paths, expected, smoke)
        print(json.dumps({"phase": "registered", "coordinates": len(registrations)}), flush=True)
        for index, (key, path) in enumerate(paths.items(), 1):
            for phase, receipt in PHASES.items():
                if not (path / receipt).is_file():
                    invoke(logs, resources, path.relative_to(output), phase, ["--output", str(path), *inputs])
            print(json.dumps({"phase": "replayed", "completed": index, "coordinate": list(key)}), flush=True)
        resources.phase = "matrix-verification"
        resources.check()
        fits = fit_evidence(output, paths, registrations, smoke)
        peak = int((scope / "memory.peak").read_text().strip())
        enforce_record(logs / f"scope-{time.time_ns()}.json",
                       {"memory_peak_bytes": peak, "memory_max_bytes": MEMORY_CAP, "memory_swap_max_bytes": SWAP_CAP})
        manifest = output / "manifest.json"
        enforce_record(manifest, {"schema": SCHEMA, "preflight_sha256": axis_run.digest(preflight_path), "fits": fits})
        return {"status": "PASS", "fits": len(fits), "manifest_sha256": axis_run.digest(manifest)}


def verify_matrix(root, smoke):
    """Read verified fresh fits without requiring their original checkout location."""
    root = Path(root)
    manifest = axis_run.read_json(root / "manifest.json")
    preflight = axis_run.read_json(root / "preflight.json")
    fresh = tuple(map(tuple, preflight["fresh_fits"]))
    if (manifest["schema"] != SCHEMA or preflight["schema"] != SCHEMA or preflight["smoke"] is not smoke
            or preflight["registered_roster"] != [list(key) for key in coordinates(smoke)]
            or manifest["preflight_sha256"] != axis_run.digest(root / "preflight.json")
            or len(fresh) != len(set(fresh)) or not fresh or not set(fresh) <= set(coordinates(smoke))
            or (smoke and fresh != coordinates(True))):
        raise ValueError("invalid recorded Step 08 matrix")
    paths = run_paths(root, fresh)
    registrations = {key: axis_run.digest(path / "registration.json") for key, path in paths.items()}
    if fit_evidence(root, paths, registrations, smoke) != manifest["fits"]:
        raise ValueError("sealed Step 08 matrix differs from its ledger")
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_inputs(parser)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    print(json.dumps(execute(args.output, args.logs, args.training.resolve(), args.representations.resolve(),
                             args.source.resolve(), args.smoke), sort_keys=True))


if __name__ == "__main__":
    main()
