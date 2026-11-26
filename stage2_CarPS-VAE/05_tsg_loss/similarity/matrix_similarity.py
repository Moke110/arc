"""Preregister, fit and independently replay both formal Delta similarity sweeps."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

from evaluate_similarity import coordinates, read_run
from matrix import (SAMPLE_SECONDS, ResourceMonitor, digest_values, enforce_record,
                    register_coordinate, scope_directory)
from matrix_reconstruction import PHASE_OUTPUTS, signal_process_group
from run_similarity import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                            EXECUTION, PHASES, axis_run, bind)
from similarity_objective import METRICS, SimilarityObjective

SCHEMA = "stage2-tenv-transigen-bounded-delta-similarity-matrix-v1"


def run_paths(root):
    return {(metric, alpha, fold): root / "runs" / metric / f"alpha-{alpha}" / f"fold-{fold}"
            for metric in METRICS for alpha, fold in coordinates()}


def preflight(training, representations, source):
    inputs = {"training": str(training), "representations": str(representations), "source": str(source)}
    with axis_run.open_inputs(training, representations) as (data, vectors):
        folds = {(metric, fold): bind({"axis": AXIS, "metric": metric, "fold": fold, "alpha": "0",
                                       "smoke": False, "device": "cuda", **inputs}, data, vectors)
                 for metric in METRICS for fold in range(3)}
    expected = {(metric, alpha, fold): {
        **folds[metric, fold], "objective": SimilarityObjective(metric, Fraction(alpha)).identity(),
    } for metric, alpha, fold in run_paths(Path("."))}
    evidence = {
        "schema": SCHEMA, "driver_sha256": axis_run.digest(Path(__file__)),
        "input_identity": digest_values(list(folds.values())), "seed": axis_run.SEED,
        "scenario": "cold-drug", "axis": AXIS, "metrics": list(METRICS),
        "folds": [
            {"fold": fold, "counts": {part: len(ids) for part, ids in config["conditions"].items()},
             "ordered_conditions_sha256": {part: digest_values(ids) for part, ids in config["conditions"].items()}}
            for (metric, fold), config in folds.items() if metric == METRICS[0]
        ],
    }
    return expected, evidence


def check_registrations(paths, expected):
    nonces, verified = set(), {}
    for key, run in paths.items():
        record = run / "registration.json"
        config = axis_run.read_json(record)
        args, nonce = config["arguments"], config["nonce"]
        identity = {field: value for field, value in config.items() if field not in {"arguments", "nonce"}}
        if ((args["metric"], args["alpha"], args["fold"]) != key or type(args["fold"]) is not int
                or args["smoke"] is not False or (args["axis"], args["device"]) != (AXIS, "cuda")
                or nonce in nonces or identity != expected[key]):
            raise ValueError(f"formal registration differs from preflight: {key}")
        nonces.add(nonce)
        verified[key] = axis_run.digest(record)
    return verified


def invoke(logs, resources, coordinate, phase, args):
    resources.phase = phase
    resources.check()
    directory = logs / coordinate
    directory.mkdir(parents=True, exist_ok=True)
    base = directory / f"{phase}-{time.time_ns()}"
    command = ["/usr/bin/time", "-v", "-o", str(base.with_suffix(".resources")),
               sys.executable, str(EXECUTION / "run_similarity.py"), phase, *args]
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
                    signal_process_group(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=8)
                    except subprocess.TimeoutExpired:
                        signal_process_group(process.pid, signal.SIGKILL)
                        process.wait()
                raise
    if code:
        raise RuntimeError(f"{coordinate} {phase} failed, see {base.with_suffix('.log')}")


def execute(root, logs, training, representations, source):
    if logs.is_relative_to(root) or root.is_relative_to(logs):
        raise ValueError("matrix and audit logs require separate directories")
    logs.mkdir(parents=True, exist_ok=True)
    scope = scope_directory()
    with ResourceMonitor(logs, scope) as resources, axis_run.publication_lock(root) as output:
        output.mkdir(exist_ok=True)
        paths = run_paths(output)
        expected, evidence = preflight(training, representations, source)
        preflight_path = output / "preflight.json"
        enforce_record(preflight_path, evidence)
        inputs = ["--training", str(training), "--representations", str(representations), "--source", str(source)]

        for (metric, alpha, fold), path in paths.items():
            if not (path / "registration.json").exists():
                if path.is_dir() and not any(path.iterdir()):
                    path.rmdir()

                def register(staged):
                    invoke(logs, resources, path.relative_to(output), "register",
                           ["--metric", metric, "--alpha", alpha, "--fold", str(fold),
                            "--output", str(staged), *inputs])

                register_coordinate(path, output, register)
        registrations = check_registrations(paths, expected)
        print(json.dumps({"phase": "registered", "coordinates": len(registrations)}), flush=True)

        fits = []
        for (metric, alpha, fold), path in paths.items():
            key = metric, alpha, fold
            relative = path.relative_to(output)
            for phase, receipt in PHASES.items():
                if (path / receipt).is_file():
                    axis_run.verified_receipt(path, receipt, registrations[key], PHASE_OUTPUTS[phase])
                else:
                    invoke(logs, resources, relative, phase, ["--metric", metric, "--output", str(path), *inputs])
            resources.phase = "verify"
            resources.check()
            verified_key, (_, _, manifest, training_record, _) = read_run(path, False, metric)
            if verified_key != (alpha, fold) or manifest["registration_sha256"] != registrations[key]:
                raise ValueError(f"sealed run differs from the formal registration: {relative}")
            fits.append({"path": relative.as_posix(), "metric": metric, "alpha": alpha, "fold": fold,
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
