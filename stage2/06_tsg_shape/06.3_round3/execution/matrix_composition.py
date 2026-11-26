"""Preregister, fit and independently replay the Step 06.3 CUDA coordinates."""

import argparse
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

from verify_composition import read_run
from run_composition import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                       PHASES, axis_run, bind)
from composition_model import STEP
from composition_roster import coordinates, fit_coordinates
from matrix import (SAMPLE_SECONDS, ResourceMonitor, digest_values, enforce_record,
                    register_coordinate, scope_directory)
from matrix_sweep import signal_process_group

SCHEMA = "stage2-tenv-transigen-architecture-composition-matrix-v1"



def run_paths(root, smoke=False):
    return {(shape, fold, seed): root / "runs" / shape / f"fold-{fold}" / f"seed-{seed}"
            for shape, fold, seed in fit_coordinates(smoke)}


def preflight(training, representations, source, smoke):
    inputs = {"training": str(training), "representations": str(representations), "source": str(source)}
    expected, folds = {}, {}
    with axis_run.open_inputs(training, representations) as (data, vectors):
        for shape, fold, seed in fit_coordinates(smoke):
            config = bind({"axis": AXIS, "fold": fold, "shape": shape, "seed": seed,
                           "smoke": smoke, "device": "cuda", **inputs}, data, vectors)
            expected[shape, fold, seed] = config
            folds.setdefault(fold, {
                "fold": fold,
                "ordered_conditions_sha256": {part: digest_values(ids)
                                               for part, ids in config["conditions"].items()},
                "counts": {part: len(ids) for part, ids in config["conditions"].items()},
            })
    return expected, {"schema": SCHEMA, "axis": AXIS, "smoke": smoke,
                      "coordinates": [list(key) for key in coordinates(smoke)],
                      "new_fits": [list(key) for key in fit_coordinates(smoke)],
                      "expected_config_sha256": digest_values(list(expected.values())),
                      "folds": list(folds.values())}


def verify_matrix(root, smoke):
    paths = run_paths(root, smoke)
    receipt = axis_run.read_json(root / "manifest.json")
    preflight = axis_run.read_json(root / "preflight.json")
    if (receipt["schema"] != SCHEMA
            or receipt["preflight_sha256"] != axis_run.digest(root / "preflight.json")
            or preflight["coordinates"] != [list(key) for key in coordinates(smoke)]
            or preflight["new_fits"] != [list(key) for key in paths]
            or [(fit["shape"], fit["fold"], fit["seed"]) for fit in receipt["fits"]] != list(paths)):
        raise ValueError("matrix receipt has an invalid preflight or coordinate roster")
    fits = dict(zip(paths, receipt["fits"], strict=True))
    for key, path in paths.items():
        fit = fits[key]
        registration = path / "registration.json"
        if (fit["path"] != path.relative_to(root).as_posix()
                or fit["run_manifest_sha256"] != axis_run.digest(path / "manifest.json")
                or fit["registration_sha256"] != axis_run.digest(registration)):
            raise ValueError(f"sealed matrix fit differs from its ledger: {key}")
    return paths, fits


def check_registrations(paths, expected, smoke):
    if tuple(paths) != fit_coordinates(smoke) or paths.keys() != expected.keys():
        raise ValueError("formal registration plan is incomplete")
    seen, hashes = set(), {}
    for key, run in paths.items():
        path = run / "registration.json"
        config = axis_run.read_json(path)
        args, nonce = config["arguments"], config["nonce"]
        identity = {field: value for field, value in config.items() if field not in {"arguments", "nonce"}}
        if ((args["shape"], args["fold"], args["seed"]) != key
                or type(args["fold"]) is not int or type(args["seed"]) is not int
                or (args["axis"], args["device"], args["smoke"]) != (AXIS, "cuda", smoke)
                or args["smoke"] is not smoke or nonce in seen or identity != expected[key]):
            raise ValueError(f"formal registration differs from preflight: {key}")
        seen.add(nonce)
        hashes[key] = axis_run.digest(path)
    return hashes


def invoke(logs, resources, coordinate, phase, args):
    resources.phase = phase
    resources.check()
    directory = logs / coordinate
    directory.mkdir(parents=True, exist_ok=True)
    attempt = directory / f"{phase}-{time.time_ns()}"
    command = ["/usr/bin/time", "-v", "-o", str(attempt.with_suffix(".resources")),
               sys.executable, str(STEP / "execution/run_composition.py"), phase, *args]
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
        raise RuntimeError(f"{coordinate} {phase} failed, see {attempt.with_suffix('.log')}")


def execute(root, logs, training, representations, source, smoke):
    if logs.is_relative_to(root) or root.is_relative_to(logs):
        raise ValueError("matrix and audit logs require separate directories")
    logs.mkdir(parents=True, exist_ok=True)
    scope = scope_directory()
    inputs = ["--training", str(training), "--representations", str(representations), "--source", str(source)]
    with (ResourceMonitor(logs, scope) as resources, axis_run.publication_lock(root) as output):
        output.mkdir(exist_ok=True)
        paths = run_paths(output, smoke)
        expected, evidence = preflight(training, representations, source, smoke)
        preflight_path = output / "preflight.json"
        enforce_record(preflight_path, evidence)
        for (shape, fold, seed), path in paths.items():
            if not (path / "registration.json").exists():
                def register(staged):
                    invoke(logs, resources, path.relative_to(output), "register",
                           ["--shape", shape, "--fold", str(fold), "--seed", str(seed),
                            "--output", str(staged), *inputs, *(["--smoke"] if smoke else [])])
                register_coordinate(path, output, register)
        registrations = check_registrations(paths, expected, smoke)
        print(json.dumps({"phase": "registered", "coordinates": len(registrations)}), flush=True)

        fits = []
        for key, path in paths.items():
            relative = path.relative_to(output)
            for phase, receipt in PHASES.items():
                if not (path / receipt).is_file():
                    invoke(logs, resources, relative, phase, ["--output", str(path), *inputs])
            resources.phase = "verify"
            resources.check()
            verified_key, (_, _, manifest, training_record, _) = read_run(path, smoke)
            if verified_key != key or manifest["registration_sha256"] != registrations[key]:
                raise ValueError(f"sealed run differs from the formal registration: {relative}")
            fits.append({"path": relative.as_posix(), "shape": key[0], "fold": key[1], "seed": key[2],
                         "registration_sha256": registrations[key],
                         "run_manifest_sha256": axis_run.digest(path / "manifest.json"),
                         "best_epoch": training_record["best_epoch"], "stopping_epoch": training_record["epochs"],
                         "train_wall_seconds": training_record["wall_seconds"],
                         "train_peak_rss_kib": training_record["peak_rss_kib"],
                         "train_cuda_peak_allocated_bytes": training_record["cuda_peak_allocated_bytes"]})
            print(json.dumps({"phase": "complete", "completed": len(fits), **fits[-1]}), flush=True)
        resources.check()
        resources.phase = "matrix-complete"
        enforce_record(logs / f"scope-{time.time_ns()}.json",
                       {"memory_peak_bytes": int((scope / "memory.peak").read_text())})
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
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    result = execute(args.output.resolve(), args.logs.resolve(), args.training.resolve(),
                     args.representations.resolve(), args.source.resolve(), args.smoke)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
