"""Preregister and run verified CarPS-VAE prediction-axis CUDA coordinates."""

import argparse
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

from axis_protocol import AXES, REFERENCE_AXIS, coordinates
from axis_runs import PHASES, add_inputs, arguments, axis_run, bind, read_run, step08_config
from chem_runs import open_inputs
from matrix import SAMPLE_SECONDS, ResourceMonitor, digest_values, enforce_record, register_coordinate
from matrix_measured import MEMORY_CAP, SWAP_CAP, scope_directory
from matrix_sweep import signal_process_group
from tenv_axes import tenv_axis_output

SCHEMA = "carps-vae-unimol-axis-matrix-v1"
RUNNER = Path(__file__).with_name("axis_runs.py")


def run_paths(output, keys):
    return {(axis, fold, seed): output / "runs" / axis / f"fold-{fold}" / f"seed-{seed}"
            for axis, fold, seed in keys}


def preflight(step03, representations, source, bounded):
    """Bind all 72 formal coordinates, even for the eight-coordinate bounded gate."""
    formal = coordinates()
    bounded_keys = set(coordinates(True))
    references = [key for key in formal if key[0] == REFERENCE_AXIS]
    reference_configs = {
        (fold, seed, False): step08_config(fold, seed, step03, representations, source, False)
        for _, fold, seed in references
    }
    if bounded:
        reference_configs[0, 17, True] = step08_config(0, 17, step03, representations, source, True)
    enrolled, registered = {}, {}
    membership = None
    for axis in AXES:
        training = tenv_axis_output(axis, step03)
        with open_inputs(training, representations) as (data, vectors):
            outputs = data.release.manifest["outputs"]
            current = {field: outputs[field] for field in ("members.csv", "exclusions.csv")}
            if membership is None:
                membership = current
            if current != membership:
                raise ValueError(f"{axis} differs from registered membership and exclusions")
            for key in (key for key in formal if key[0] == axis):
                _, fold, seed = key
                enrolled[key] = bind(key, data, vectors, source, False, reference_configs[fold, seed, False])
                if bounded and key in bounded_keys:
                    registered[key] = bind(key, data, vectors, source, True, reference_configs[fold, seed, True])
    configs = registered if bounded else enrolled
    populations = [{part: {"count": len(ids), "conditions_sha256": digest_values(ids)}
                    for part, ids in enrolled[REFERENCE_AXIS, fold, 17]["conditions"].items()}
                   for fold in range(3)]
    return configs, {
        "schema": SCHEMA, "bounded": bounded,
        "formal_roster": [list(key) for key in formal],
        "step08_reference_roster": [list(key) for key in references],
        "registered_roster": [list(key) for key in configs],
        "formal_configs_sha256": digest_values(list(enrolled.values())),
        "registered_configs_sha256": digest_values(list(configs.values())),
        "folds": populations, "membership_sha256": membership,
        "training_manifests_sha256": {axis: enrolled[axis, 0, 17]["training_sha256"] for axis in AXES},
        "unimol": enrolled[REFERENCE_AXIS, 0, 17]["representation"],
    }


def invoke(logs, resources, coordinate, phase, options):
    resources.phase = phase
    resources.check()
    directory = logs / coordinate
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


def verify_pairing(fits):
    groups = {}
    for fit in fits:
        _, fold, seed = fit["coordinate"]
        groups.setdefault((fold, seed), []).append(fit)
    for key, group in groups.items():
        reference = max(group, key=lambda fit: len(fit["epoch_streams"]))
        expected = {field: value for field, value in reference["initialization"].items()
                    if field != "axis_projection_shapes"}
        for fit in group:
            observed = {field: value for field, value in fit["initialization"].items()
                        if field != "axis_projection_shapes"}
            streams = fit["epoch_streams"]
            if observed != expected or streams != reference["epoch_streams"][:len(streams)]:
                raise ValueError(f"paired initialized parameters or training streams differ at {key}")


def fit_evidence(output, paths, bounded):
    fits = []
    for key, path in paths.items():
        observed, (_, config, replay, trained, _) = read_run(path, bounded)
        if observed != key:
            raise ValueError(f"sealed run does not match its coordinate: {key}")
        width = config["model"]["n_genes"]
        shapes = {"encoder_x1.0.weight": [1200, width], "encoder_delta.0.weight": [1200, width],
                  "decoder_x1.4.weight": [width, 800], "decoder_x1.4.bias": [width],
                  "decoder_delta.4.weight": [width, 800], "decoder_delta.4.bias": [width]}
        if trained["initialization"]["axis_projection_shapes"] != shapes:
            raise ValueError(f"{key} does not resize the four registered gene projections")
        fits.append({
            "coordinate": list(key), "path": path.relative_to(output).as_posix(),
            "registration_sha256": replay["registration_sha256"],
            "manifest_sha256": axis_run.digest(path / "manifest.json"),
            "initialization": trained["initialization"], "epoch_streams": trained["epoch_streams"],
            "best_epoch": trained["best_epoch"], "epochs": trained["epochs"],
            "train_peak_rss_kib": trained["peak_rss_kib"],
            "train_cuda_peak_allocated_bytes": trained["cuda_peak_allocated_bytes"],
        })
    verify_pairing(fits)
    return fits


def execute(root, logs, step03, representations, source, bounded=False, fresh=None):
    """Register the study before fitting and replaying each selected fresh coordinate."""
    root, logs = Path(root).resolve(), Path(logs).resolve()
    if logs.is_relative_to(root) or root.is_relative_to(logs):
        raise ValueError("matrix and audit log directories must be separate")
    roster = coordinates(bounded)
    fresh = roster if fresh is None else tuple(fresh)
    if (len(fresh) != len(set(fresh)) or not set(fresh) <= set(roster) or not fresh
            or bounded and fresh != roster):
        raise ValueError("fresh fits must be a unique nonempty subset, or all eight bounded coordinates")
    scope = scope_directory()
    inputs = ["--step03", str(step03), "--representations", str(representations), "--source", str(source)]
    logs.mkdir(parents=True, exist_ok=True)
    with axis_run.publication_lock(root) as output, ResourceMonitor(logs, scope) as resources:
        output.mkdir(exist_ok=True)
        paths = run_paths(output, fresh)
        config, evidence = preflight(step03, representations, source, bounded)
        preflight_path = output / "preflight.json"
        enforce_record(preflight_path, {**evidence, "fresh_fits": [list(key) for key in fresh]})
        nonces = set()
        for (axis, fold, seed), path in paths.items():
            if not (path / "registration.json").is_file():
                def register(staged):
                    invoke(logs, resources, path.relative_to(output), "register",
                           ["--axis", axis, "--fold", str(fold), "--seed", str(seed), "--output", str(staged),
                            *inputs, *(["--bounded"] if bounded else [])])
                register_coordinate(path, output, register)
            observed = axis_run.read_json(path / "registration.json")
            nonce, args = observed.pop("nonce"), observed.pop("arguments")
            key = axis, fold, seed
            if nonce in nonces or args != arguments(key, bounded) or observed != config[key]:
                raise ValueError(f"registration differs from frozen preflight: {key}")
            nonces.add(nonce)
        for index, (key, path) in enumerate(paths.items(), 1):
            for phase, receipt in PHASES.items():
                if not (path / receipt).is_file():
                    invoke(logs, resources, path.relative_to(output), phase, ["--output", str(path), *inputs])
            print(json.dumps({"phase": "replayed", "completed": index, "coordinate": list(key)}), flush=True)
        resources.phase = "matrix-verification"
        resources.check()
        fits = fit_evidence(output, paths, bounded)
        peak = int((scope / "memory.peak").read_text().strip())
        enforce_record(logs / f"scope-{time.time_ns()}.json",
                       {"memory_peak_bytes": peak, "memory_max_bytes": MEMORY_CAP, "memory_swap_max_bytes": SWAP_CAP})
        manifest = output / "manifest.json"
        enforce_record(manifest, {"schema": SCHEMA, "preflight_sha256": axis_run.digest(preflight_path), "fits": fits})
        return {"status": "PASS", "fits": len(fits), "manifest_sha256": axis_run.digest(manifest)}


def verify_matrix(output, bounded):
    """Verify the registered fresh fit roster, lineage and independent replay."""
    output = Path(output)
    manifest = axis_run.read_json(output / "manifest.json")
    preflight = axis_run.read_json(output / "preflight.json")
    fresh = tuple(map(tuple, preflight["fresh_fits"]))
    if (manifest["schema"] != SCHEMA or preflight["schema"] != SCHEMA or preflight["bounded"] is not bounded
            or preflight["formal_roster"] != [list(key) for key in coordinates()]
            or preflight["step08_reference_roster"] != [list(key) for key in coordinates() if key[0] == REFERENCE_AXIS]
            or preflight["registered_roster"] != [list(key) for key in coordinates(bounded)]
            or manifest["preflight_sha256"] != axis_run.digest(output / "preflight.json")
            or len(fresh) != len(set(fresh)) or not fresh or not set(fresh) <= set(coordinates(bounded))
            or bounded and fresh != coordinates(True)):
        raise ValueError("invalid Step 09 matrix record")
    paths = run_paths(output, fresh)
    if fit_evidence(output, paths, bounded) != manifest["fits"]:
        raise ValueError("sealed Step 09 fits differ from their matrix record")
    return dict(read_run(path, bounded) for path in paths.values())


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
