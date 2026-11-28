"""Audit the Step 05 objective traces and the capped CUDA resource record."""

import argparse
from itertools import groupby
import json
import math
from pathlib import Path

import pandas as pd

from matrix_measured import GIB, MEMORY_CAP, SWAP_CAP, verify_matrix
from measured_model import axis_run
from measured_roster import CONTROL_ONLY
from module_roster import fit_coordinates, protocol_roster
from util.training.resource_samples import gaps_summary, read_samples

RELATIVE_TOLERANCE = 1e-6
ABSOLUTE_TOLERANCE = 1e-4
EVALUATION_SEEDS = (22617, 100042)


def audit_objectives(matrix, smoke):
    paths, _ = verify_matrix(matrix, smoke)
    count, max_error = 0, 0.
    condition_orders, identities = {}, {}
    for (shape, fold), group in groupby(paths, key=lambda key: key[:2]):
        prediction_hashes = set()
        for coordinate in group:
            path = paths[coordinate]
            config = axis_run.read_json(path / "registration.json")
            if (config["seed"], config["validation_seed"], config["test_seed"]) != (coordinate[2], *EVALUATION_SEEDS):
                raise ValueError(f"training or evaluation seed differs in {coordinate}")
            order = config["conditions"]
            identity = (
                config["partitions"], config["model"]["transform"], config["training_sha256"],
                config["kpgt_manifest_sha256"], config["kpgt_vectors_sha256"], config["axis_sha256"],
            )
            if (condition_orders.setdefault(fold, order) != order or identities.setdefault(fold, identity) != identity):
                raise ValueError(f"registered populations, input identity or dose transform differ in {coordinate}")
            prediction_hashes.add(axis_run.digest(path / "prediction.npy"))

            trace = pd.read_csv(path / "learning_curve.csv", float_precision="round_trip")
            if smoke and trace.epoch.tolist() != [1, 2]:
                raise ValueError(f"expected two bounded epochs in {coordinate}")
            control_only = shape == CONTROL_ONLY
            for partition in ("train", "validation"):
                for row in trace.to_dict(orient="records"):
                    def component(name):
                        return row[f"{partition}_{name}"]

                    k1, k2, kp = [component(name) for name in ("control_kl", "treatment_kl", "prediction_kl")]
                    measured_mse = component("treatment_mse")
                    if (k1 <= 0 or (control_only and any((k2, kp, measured_mse)))
                            or (not control_only and min(k2, kp, measured_mse) <= 0)):
                        raise ValueError(f"objective roles differ in {coordinate}, {partition}")
                    objective = (0.1 * math.fsum((k1, k2, kp)) + component("control_mse") + measured_mse
                                 + component("delta_mse"))
                    if component("delta_l1") != 0:
                        raise ValueError(f"Unexpected L1 contribution in {coordinate}, {partition}")
                    observed = component("loss")
                    if not math.isclose(observed, objective, rel_tol=RELATIVE_TOLERANCE,
                                        abs_tol=ABSOLUTE_TOLERANCE):
                        raise ValueError(f"weighted objective differs in {coordinate}, {partition}")
                    max_error = max(max_error, abs(observed - objective) / objective)
                    count += 1
        if len(prediction_hashes) != (1 if smoke else len(protocol_roster()[1])):
            raise ValueError(f"training seeds yield duplicate predictions for {shape}, fold {fold}")
    return {
        "status": "PASS", "matrix_manifest_sha256": axis_run.digest(matrix / "manifest.json"),
        "fits": len(paths), "objective_rows": count, "objective_max_relative_error": max_error,
        "objective_relative_tolerance": RELATIVE_TOLERANCE, "objective_absolute_tolerance": ABSOLUTE_TOLERANCE,
        "evaluation_seeds": {"validation": EVALUATION_SEEDS[0], "test": EVALUATION_SEEDS[1]},
        "population_invariant_across_arms_and_seeds": True,
        "prediction_hashes_distinct_across_training_seeds": True,
    }


def audit_resources(audit, coordinator, destination, smoke):
    matrix_samples = read_samples(audit / "resource-samples.jsonl")
    scope_files = sorted(audit.glob("scope-*.json"))
    if not scope_files or not any(row["phase"] == "train" for row in matrix_samples):
        raise ValueError("resource audit requires completed scope and samples during training")
    scope_records = [axis_run.read_json(path) for path in scope_files]
    if any((scope["memory_max_bytes"], scope["memory_swap_max_bytes"]) != (MEMORY_CAP, SWAP_CAP)
           for scope in scope_records):
        raise ValueError("matrix resource scope differs from the registered limits")
    matrix_gaps = gaps_summary(matrix_samples)
    start = matrix_samples[0]["monotonic_seconds"]
    final_sample = matrix_samples[-1]
    finished = int(scope_files[-1].stem.removeprefix("scope-")) / 1e9
    end = final_sample["monotonic_seconds"] + max(0, finished - final_sample["unix_seconds"])
    coordinator_samples = read_samples(coordinator)
    before = [row for row in coordinator_samples if row["monotonic_seconds"] <= start]
    between = [row for row in coordinator_samples if start < row["monotonic_seconds"] < end]
    after = [row for row in coordinator_samples if row["monotonic_seconds"] >= end]
    if not before or not after:
        raise ValueError("coordinator samples must bracket the completed matrix")
    coverage = [before[-1], *between, after[0]]
    coordinator_gaps = gaps_summary(coverage)
    # The matrix timestamps follow its GPU query. The coordinator timestamps
    # precede that query and establish coverage independently of sensor latency.
    if coordinator_gaps["gaps_over_12_seconds"]:
        raise ValueError("coordinator records require ten-second coverage")
    host = [row["windows_host_total_and_available_kib"] for row in coverage]
    if any(not isinstance(value, str) for value in host):
        raise ValueError("coordinator samples require Windows host memory readings")
    host_free = min(int(value.split(",")[1]) for value in host)
    wsl_free = min(row["wsl"]["MemAvailable_kib"] for row in coverage)
    peak = max(record["memory_peak_bytes"] for record in scope_records)
    ram = min(row["wsl_available_kib"] for row in matrix_samples)
    gpu = min(row["gpu_total_mib"] - row["gpu_used_mib"] for row in matrix_samples)
    disk = min(row["disk_free_bytes"] for row in matrix_samples)
    if min(ram, wsl_free) < 12 * GIB // 1024 or gpu < 512 or disk < GIB or peak > MEMORY_CAP:
        raise ValueError("matrix exceeded resource scope or available headroom")

    evidence = list(scope_files)
    phase_peaks = {}
    for shape, fold, seed in fit_coordinates(smoke):
        directory = audit / "runs" / shape / f"fold-{fold}" / f"seed-{seed}"
        for phase in ("register", "train", "predict", "replay"):
            measurements = sorted(directory.glob(f"{phase}-*.resources"))
            if not measurements:
                raise ValueError(f"missing {phase} resource measurements for {(shape, fold, seed)}")
            for measurement in measurements:
                evidence.extend((measurement, measurement.with_suffix(".log")))
                for line in measurement.read_text().splitlines():
                    if "Maximum resident set size (kbytes):" in line:
                        peak_rss = int(line.rsplit(": ", 1)[1])
                        phase_peaks[phase] = max(phase_peaks.get(phase, 0), peak_rss)
    destination.write_text("".join(json.dumps(row) + "\n" for row in coverage), encoding="utf-8")
    return {
        "status": "PASS", **matrix_gaps,
        "scope_memory_max_bytes": MEMORY_CAP, "scope_swap_max_bytes": SWAP_CAP,
        "peak_scope_memory_bytes": peak, "phase_peak_rss_kib": phase_peaks,
        "minimum_wsl_available_kib": ram, "minimum_gpu_available_mib": gpu,
        "minimum_disk_free_bytes": disk,
        "matrix_samples_sha256": axis_run.digest(audit / "resource-samples.jsonl"),
        "coordinator_samples_sha256": axis_run.digest(destination),
        "phase_audit_sha256": axis_run.file_hashes(audit, [path.relative_to(audit).as_posix() for path in evidence]),
        "coordinator": {
            **coordinator_gaps, "minimum_windows_available_kib": host_free, "minimum_wsl_available_kib": wsl_free,
            "maximum_gpu_used_mib": max(int(row["gpu_total_used_mib_and_utilization_percent"].split(",")[1])
                                        for row in coverage),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    objective = commands.add_parser("objectives")
    objective.add_argument("--matrix", type=Path, required=True)
    resources = commands.add_parser("resources")
    resources.add_argument("--matrix-audit", type=Path, required=True)
    resources.add_argument("--coordinator", type=Path, required=True)
    resources.add_argument("--coverage-output", type=Path, required=True)
    for command in commands.choices.values():
        command.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    result = (audit_objectives(args.matrix, args.smoke) if args.command == "objectives"
              else audit_resources(args.matrix_audit, args.coordinator, args.coverage_output, args.smoke))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
