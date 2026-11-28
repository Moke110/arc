"""Seal the bounded CUDA gate and its input, replay, scoring and resource evidence."""

import argparse
from datetime import datetime
import json
from pathlib import Path

from dose_matrix import verify_matrix
from dose_protocol import coordinates
from dose_report import PLOT_TABLES, PRIMARY_METRICS, REPORT_METRICS
from dose_runs import axis_run
from cpst_util.comparison import FIGURES


def read_samples(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def verify_resource_samples(coordinator_path, matrix_path, final_measurement):
    coordinator, matrix = read_samples(coordinator_path), read_samples(matrix_path)
    if not coordinator or not matrix:
        raise ValueError("missing Step 10 coordinator or matrix resource samples")
    timestamps = [item["monotonic_seconds"] for item in coordinator]
    if (len(timestamps) < 2 or any(abs(stop - start - 10) > 1.0
                                  for start, stop in zip(timestamps, timestamps[1:]))
            or timestamps[0] > matrix[0]["monotonic_seconds"]
            or timestamps[-1] < matrix[-1]["monotonic_seconds"]
            or datetime.fromisoformat(coordinator[-1]["wall_time"]).timestamp() < final_measurement):
        raise ValueError("coordinator must sample every ten seconds from before fit through its last measurement")
    wsl = [item["wsl"]["MemAvailable_kib"] for item in coordinator]
    host = [int(item["windows_host_total_and_available_kib"].split(",")[1]) for item in coordinator]
    gpus = [list(map(int, item["gpu_total_used_mib_and_utilization_percent"].split(","))) for item in coordinator]
    if min(wsl) < 12 * 1024**2 or min(host) <= 0 or any(total - used < 512 for total, used, _ in gpus):
        raise ValueError("bounded fit used the last available memory headroom")
    return {"coordinator_samples": len(coordinator), "matrix_samples": len(matrix),
            "minimum_wsl_available_kib": min(wsl), "minimum_windows_available_kib": min(host),
            "maximum_gpu_used_mib": max(used for _, used, _ in gpus),
            "maximum_gpu_utilization_percent": max(utilization for _, _, utilization in gpus)}


def acceptance(bundle):
    matrix, report, figure = (bundle / name for name in ("matrix", "report", "figure"))
    runs = verify_matrix(matrix, True)
    preflight = axis_run.read_json(matrix / "preflight.json")
    report_receipt = axis_run.read_json(report / "manifest.json")
    figure_receipt = axis_run.read_json(bundle / "figure-replay.json")
    if (not report_receipt["bounded"] or report_receipt["scores"] != len(runs) * len(REPORT_METRICS)
            or report_receipt["primary_metrics"] != list(PRIMARY_METRICS)
            or report_receipt["outputs"] != axis_run.file_hashes(report, report_receipt["outputs"])
            or (figure_receipt["status"], figure_receipt["independent_csv_only_renders"]) != ("PASS", 2)
            or figure_receipt["png_sha256"] != axis_run.file_hashes(figure, FIGURES)
            or figure_receipt["csv_sha256"] != axis_run.file_hashes(figure / "data", PLOT_TABLES)):
        raise ValueError("bounded reports or figures lack verified evidence")
    scope_files = sorted((bundle / "matrix-audit").glob("scope-*.json"))
    scope = axis_run.read_json(scope_files[-1])
    if (scope["memory_max_bytes"], scope["memory_swap_max_bytes"]) != (12 * 1024**3, 8 * 1024**3):
        raise ValueError("bounded fit ran under a different resource envelope")
    resource_files = ("coordinator.jsonl", "matrix-audit/resource-samples.jsonl",
                      scope_files[-1].relative_to(bundle).as_posix(), "command.resources",
                      "report.resources", "figure.resources")
    resources = verify_resource_samples(
        bundle / resource_files[0], bundle / resource_files[1],
        (bundle / "command.resources").stat().st_mtime)
    lineage = []
    for key, (path, config, replay, training, _) in runs.items():
        lineage.append({
            "coordinate": list(key), "fit_sha256": axis_run.digest(path / "manifest.json"),
            "registration_sha256": replay["registration_sha256"], "best_epoch": training["best_epoch"],
            "validation_loss": training["validation_loss"], "test_conditions": replay["test_conditions"],
            "replay": {name: replay[name] for name in ("test_replay", "validation_replay",
                                                       "validation_objective_replay")},
            "model": config["model"],
        })
    return {
        "status": "PASS", "ticket": 377, "bounded_coordinates": [list(key) for key in coordinates(True)],
        "formal_coordinates": len(coordinates()), "formal_configs_sha256": preflight["formal_configs_sha256"],
        "bounded_fit_lineage": lineage, "full_folds": preflight["folds"],
        "matrix_manifest_sha256": axis_run.digest(matrix / "manifest.json"),
        "matrix_preflight_sha256": axis_run.digest(matrix / "preflight.json"),
        "report_manifest_sha256": axis_run.digest(report / "manifest.json"),
        "figures": figure_receipt,
        "resources": {**resources, **scope, "fit_peak_rss_kib": max(run[3]["peak_rss_kib"] for run in runs.values()),
                      "fit_peak_cuda_allocated_bytes": max(run[3]["cuda_peak_allocated_bytes"]
                                                           for run in runs.values()),
                      "files_sha256": axis_run.file_hashes(bundle, resource_files)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = acceptance(args.bundle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    axis_run.write_json(args.output, result)


if __name__ == "__main__":
    main()
