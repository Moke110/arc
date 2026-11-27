"""Preregister and verify the complete Step 10 CUDA fit matrix."""

import argparse
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

from dose_inputs import open_inputs
from dose_protocol import REFERENCE, STEP, TRANSIGEN, coordinates
from dose_runs import PHASES, add_inputs, axis_run, bind, read_run, verify_release
from matrix import ResourceMonitor, SAMPLE_SECONDS, digest_values, enforce_record, register_coordinate
from matrix_measured import MEMORY_CAP, SWAP_CAP, scope_directory
from matrix_sweep import signal_process_group

SCHEMA = "carps-vae-dose-matrix-v1"
RUNNER = Path(__file__).with_name("dose_runs.py")


def run_paths(root, bounded):
    return {key: root / "runs" / key[0] / f"fold-{key[1]}" / f"seed-{key[2]}"
            for key in coordinates(bounded)}


def invoke(logs, resources, path, phase, options):
    resources.phase = phase
    resources.check()
    directory = logs / path
    directory.mkdir(parents=True, exist_ok=True)
    attempt = directory / f"{phase}-{time.time_ns()}"
    command = ["/usr/bin/time", "-v", "-o", str(attempt.with_suffix(".resources")),
               sys.executable, str(RUNNER), phase, *options]
    with attempt.with_suffix(".log").open("xb") as output:
        with subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, start_new_session=True) as process:
            try:
                while True:
                    try:
                        status = process.wait(timeout=SAMPLE_SECONDS)
                    except subprocess.TimeoutExpired:
                        resources.check()
                        continue
                    resources.check()
                    break
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
        raise RuntimeError(f"{path}: {phase} failed, see {attempt.with_suffix('.log')}")


def verify_pairing(runs):
    """Check same-shaped initialized tensors and shared RNG draws by fold and seed."""
    groups = {}
    for (arm, fold, seed), (_, _, _, training, _) in runs.items():
        if arm != TRANSIGEN:
            groups.setdefault((fold, seed), {})[arm] = training
    for key, group in groups.items():
        reference = max(group.values(), key=lambda run: len(run["epoch_streams"]))
        streams = {name: value for name, value in group[REFERENCE]["initialization"].items()
                   if name.endswith("_rng_sha256")}
        tensors = {}
        for arm, run in group.items():
            initial = run["initialization"]
            epochs = run["epoch_streams"]
            if any(initial[name] != value for name, value in streams.items()) or (
                    epochs != reference["epoch_streams"][:len(epochs)]):
                raise ValueError(f"{key}: {arm} does not share the registered training stream")
            for name, record in initial["parameters"].items():
                identity = name, tuple(record["shape"])
                if tensors.setdefault(identity, record["sha256"]) != record["sha256"]:
                    raise ValueError(f"{key}: same-shaped initialized parameter differs for {arm}: {name}")


def preflight(step03, representations, source, bounded):
    verify_release()
    with open_inputs(step03, representations) as (data, vectors):
        formal = {key: bind(key, data, vectors, source, False) for key in coordinates()}
        selected = ({key: bind(key, data, vectors, source, True) for key in coordinates(True)}
                    if bounded else formal)
    folds = []
    for fold in range(3):
        config = formal[(REFERENCE, fold, 17)]
        folds.append({**config["population"], "fold": fold,
                      "ordered_conditions_sha256": {phase: digest_values(ids)
                                                    for phase, ids in config["conditions"].items()}})
    return selected, {
        "schema": SCHEMA, "bounded": bounded, "formal_roster": [list(key) for key in formal],
        "registered_roster": [list(key) for key in selected],
        "formal_configs_sha256": digest_values([json.dumps(config, sort_keys=True) for config in formal.values()]),
        "registered_configs_sha256": digest_values([json.dumps(config, sort_keys=True) for config in selected.values()]),
        "folds": folds,
    }


def read_fits(root, bounded):
    runs, fits = {}, []
    paths = run_paths(root, bounded)
    observed_paths = {path.parent for path in (root / "runs").glob("*/*/*/registration.json")}
    if observed_paths != set(paths.values()):
        raise ValueError("Step 10 matrix registration paths differ from its full roster")
    for key, path in paths.items():
        observed, run = read_run(path, bounded)
        if observed != key:
            raise ValueError(f"sealed run does not match its coordinate: {key}")
        runs[key] = run
        _, _, replay, trained, _ = run
        fits.append({
            "coordinate": list(key), "path": path.relative_to(root).as_posix(),
            "registration_sha256": replay["registration_sha256"],
            "manifest_sha256": axis_run.digest(path / "manifest.json"),
            "best_epoch": trained["best_epoch"], "epochs": trained["epochs"],
            "train_peak_rss_kib": trained["peak_rss_kib"],
            "train_cuda_peak_allocated_bytes": trained["cuda_peak_allocated_bytes"],
        })
    if len({run[1]["nonce"] for run in runs.values()}) != len(runs):
        raise ValueError("different Step 10 coordinates share a registration nonce")
    verify_pairing(runs)
    return runs, fits


def verify_bounded_gate(evidence):
    """Require the completed ten-arm CUDA gate for this formal roster."""
    receipt = axis_run.read_json(STEP / "log/377-bounded.json")
    if (receipt["status"] != "PASS" or receipt["ticket"] != 377
            or receipt["bounded_coordinates"] != [list(key) for key in coordinates(True)]
            or [row["coordinate"] for row in receipt["bounded_fit_lineage"]] != receipt["bounded_coordinates"]
            or any(set(row["replay"].values()) != {"exact"} for row in receipt["bounded_fit_lineage"])
            or receipt["figures"]["status"] != "PASS"
            or receipt["formal_coordinates"] != len(coordinates())
            or receipt["formal_configs_sha256"] != evidence["formal_configs_sha256"]
            or receipt["full_folds"] != evidence["folds"]):
        raise ValueError("formal fit does not match the accepted Step 10 CUDA gate")


def execute(root, logs, step03, representations, source, bounded=False):
    root, logs = root.resolve(), logs.resolve()
    if root.is_relative_to(logs) or logs.is_relative_to(root):
        raise ValueError("matrix and audit logs require separate directories")
    scope = scope_directory()
    logs.mkdir(parents=True, exist_ok=True)
    locations = ["--step03", str(step03), "--representations", str(representations), "--source", str(source)]
    with axis_run.publication_lock(root) as output, ResourceMonitor(logs, scope) as resources:
        output.mkdir(exist_ok=True)
        expected, evidence = preflight(step03, representations, source, bounded)
        if not bounded:
            verify_bounded_gate(evidence)
        preflight_path = output / "preflight.json"
        enforce_record(preflight_path, evidence)
        paths = run_paths(output, bounded)
        nonces = set()
        for key, path in paths.items():
            arm, fold, seed = key
            if not (path / "registration.json").is_file():
                def register(staged):
                    invoke(logs, resources, path.relative_to(output), "register",
                           ["--arm", arm, "--fold", str(fold), "--seed", str(seed), "--output", str(staged),
                            *locations, *(["--bounded"] if bounded else [])])
                register_coordinate(path, output, register)
            registration = axis_run.read_json(path / "registration.json")
            nonce = registration.pop("nonce")
            if nonce in nonces or registration != expected[key]:
                raise ValueError(f"run registration differs from Step 10 preflight: {key}")
            nonces.add(nonce)
        for index, (key, path) in enumerate(paths.items(), 1):
            for phase, receipt in PHASES.items():
                if not (path / receipt).is_file():
                    invoke(logs, resources, path.relative_to(output), phase, ["--output", str(path), *locations])
            print(json.dumps({"phase": "replayed", "completed": index, "coordinate": list(key)}), flush=True)
        resources.phase = "matrix-verification"
        _, fits = read_fits(output, bounded)
        peak = int((scope / "memory.peak").read_text().strip())
        enforce_record(logs / f"scope-{time.time_ns()}.json",
                       {"memory_peak_bytes": peak, "memory_max_bytes": MEMORY_CAP, "memory_swap_max_bytes": SWAP_CAP})
        manifest = output / "manifest.json"
        enforce_record(manifest, {"schema": SCHEMA, "preflight_sha256": axis_run.digest(preflight_path), "fits": fits})
        return {"status": "PASS", "fits": len(fits), "manifest_sha256": axis_run.digest(manifest)}


def verify_matrix(root, bounded):
    """Verify fit roster, independent phase receipts, paired initialization and streams."""
    preflight_path = root / "preflight.json"
    preflight = axis_run.read_json(preflight_path)
    manifest = axis_run.read_json(root / "manifest.json")
    if (preflight["schema"] != SCHEMA or preflight["bounded"] is not bounded
            or preflight["formal_roster"] != [list(key) for key in coordinates()]
            or preflight["registered_roster"] != [list(key) for key in coordinates(bounded)]
            or manifest["schema"] != SCHEMA or manifest["preflight_sha256"] != axis_run.digest(preflight_path)):
        raise ValueError("invalid Step 10 matrix or preflight record")
    runs, fits = read_fits(root, bounded)
    selected = [{name: value for name, value in run[1].items() if name != "nonce"} for run in runs.values()]
    # JSON registration files are sorted by key, so normalize the field order
    # before comparing configuration identities across independently read files.
    selected = [json.dumps(config, sort_keys=True) for config in selected]
    if (manifest["fits"] != fits or preflight["registered_configs_sha256"] != digest_values(selected)):
        raise ValueError("sealed Step 10 fits differ from the preregistered matrix")
    return runs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_inputs(parser)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--bounded", action="store_true")
    args = parser.parse_args()
    print(json.dumps(execute(args.output, args.logs, args.step03.resolve(), args.representations.resolve(),
                             args.source.resolve(), args.bounded), sort_keys=True))


if __name__ == "__main__":
    main()
