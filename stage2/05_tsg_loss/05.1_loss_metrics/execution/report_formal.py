"""Verify the loss matrix and prepare its accepted figures and evidence."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import time

import pandas as pd

from audit_learning_curves import audit as audit_curves
from evaluate_loss import REPORT_SCHEMA
from matrix import GIB, SAMPLE_SECONDS, SCHEMA, run_paths
from objective import NATIVE_ALPHA
from run_loss import AXIS_STEP, DEFAULT_TRAINING, STEP, axis_run
from training_data import TENV_DEFAULT

ANCHOR_ARRAYS = ("validation_prediction.npy", "prediction.npy", "predicted_delta.npy", "measured_delta.npy")
TRACE_FIELDS = ("epoch", "train_loss", "validation_loss", "validation_mse")
RECEIPT = Path("log/330-formal.json")
CURVE_RECEIPT = Path("log/330-curves.json")
METRIC_PNGS = ("delta-metrics.png", "delta-metrics-improvement.png")
CURVE_PNG = "learning-curves.png"
PLOT_FOLD_SCORES = "figure/data/fold-scores.csv"
PROJECTIONS = {"delta-metrics.csv": "figure/data/delta-metrics.csv",
               "epoch-objectives.csv": "figure/data/epoch-objectives.csv",
               "fold-scores.csv": "log/fold-scores.csv", "runs.csv": "log/runs.csv",
               "bandwidths.csv": "log/bandwidths.csv"}


def matrix_inputs(matrix):
    """Check the matrix's sealed input roster before scoring any prediction."""
    manifest = axis_run.read_json(matrix / "manifest.json")
    paths = run_paths(matrix)
    fits = {(fit["loss"], fit["alpha"], fit["fold"]): fit for fit in manifest["fits"]}
    if (manifest["schema"] != SCHEMA or len(manifest["fits"]) != len(paths) or fits.keys() != paths.keys()
            or manifest["preflight_sha256"] != axis_run.digest(matrix / "preflight.json")
            or set(paths.values()) != {path for path in (matrix / "runs").glob("*/*/*") if path.is_dir()}):
        raise ValueError("formal matrix must bind the exact 36-fit inventory")
    for key, path in paths.items():
        fit = fits[key]
        if (fit["path"] != path.relative_to(matrix).as_posix()
                or fit["registration_sha256"] != axis_run.digest(path / "registration.json")
                or fit["run_manifest_sha256"] != axis_run.digest(path / "manifest.json")):
            raise ValueError(f"matrix receipt differs for {key}")
    return paths, fits


def verify_native_control(candidate, previous):
    """Match the Step 04 native objective trace and its checkpoint predictions."""
    traces = [pd.read_csv(run / "learning_curve.csv", usecols=TRACE_FIELDS, float_precision="round_trip")
              for run in (candidate, previous)]
    training, original = (axis_run.read_json(run / "training.json") for run in (candidate, previous))
    epochs = {key: training[key] for key in ("best_epoch", "epochs")}
    if not traces[0].equals(traces[1]) or any(original[key] != value for key, value in epochs.items()):
        raise ValueError(f"Step 04 native control trace or selected epoch differs: {candidate.name}")
    hashes = axis_run.file_hashes(candidate, ANCHOR_ARRAYS)
    for name in ANCHOR_ARRAYS:
        if hashes[name] != axis_run.digest(previous / name):
            raise ValueError(f"Step 04 native control differs for {candidate.name}/{name}")
    return {**epochs, "arrays_sha256": hashes}


def resources(audit):
    sample_file = audit / "resource-samples.jsonl"
    samples = [json.loads(line) for line in sample_file.read_text().splitlines()]
    gaps = [later["monotonic_seconds"] - earlier["monotonic_seconds"]
            for earlier, later in zip(samples, samples[1:])
            if later["monitor_session"] == earlier["monitor_session"]]
    if any(gap <= 0 for gap in gaps):
        raise ValueError("formal resource samples are out of order")
    min_ram = min(sample["wsl_available_kib"] for sample in samples)
    min_disk = min(sample["disk_free_bytes"] for sample in samples)
    min_gpu = min(sample["gpu_total_mib"] - sample["gpu_used_mib"] for sample in samples)
    scope_peak = max(axis_run.read_json(path)["memory_peak_bytes"] for path in audit.glob("scope-*.json"))
    if min_ram < 12 * GIB // 1024 or min_disk < GIB or min_gpu < 512 or scope_peak > 10 * GIB:
        raise ValueError("formal resource audit exceeded its runtime limits")
    phase_logs = {}
    for path in run_paths(audit).values():
        for phase in ("register", "train", "predict", "replay"):
            measurements = sorted(path.glob(f"{phase}-*.resources"))
            if not measurements:
                raise ValueError(f"formal resource audit lacks {phase} measurements: {path}")
            for measurement in measurements:
                for record in (measurement, measurement.with_suffix(".log")):
                    phase_logs[record.relative_to(audit).as_posix()] = axis_run.digest(record)
    return {"sample_interval_seconds": SAMPLE_SECONDS, "samples": len(samples),
            "monitor_sessions": len({sample["monitor_session"] for sample in samples}),
            "maximum_observed_sample_gap_seconds": round(max(gaps, default=0), 3),
            "sample_gaps_over_12_seconds": sum(gap > 12 for gap in gaps),
            "minimum_wsl_available_kib": min_ram, "minimum_disk_free_bytes": min_disk,
            "maximum_gpu_used_mib": max(sample["gpu_used_mib"] for sample in samples),
            "minimum_gpu_available_mib": min_gpu, "scope_memory_peak_bytes": scope_peak,
            "phase_audit_sha256": phase_logs}


def command(script, args, activity, name):
    base = activity / f"{name}-{time.time_ns()}"
    with base.with_suffix(".log").open("xb") as output:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", str(base.with_suffix(".resources")),
                                 sys.executable, str(script), *args], stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{name} failed, see {base.with_suffix('.log')}")


def verify_evaluation(report, fits):
    evidence = axis_run.read_json(report / "manifest.json")
    if (evidence["schema"], evidence["status"], evidence["runs"], evidence["smoke"]) != (
            REPORT_SCHEMA, "PASS", len(fits), False):
        raise ValueError("formal scoring manifest differs from the formal matrix")
    runs = pd.read_csv(report / "runs.csv", dtype={"alpha": str})
    if len(runs) != len(fits):
        raise ValueError("formal scoring report has the wrong fit count")
    for row in runs.itertuples(index=False):
        fit = fits[row.loss, row.alpha, row.fold]
        if (row.run_sha256, row.manifest_sha256, row.best_epoch, row.stopping_epoch) != (
                fit["registration_sha256"], fit["run_manifest_sha256"], fit["best_epoch"], fit["stopping_epoch"]):
            raise ValueError(f"scoring/run identity mismatch: {fit['path']}")
    return evidence


def publish(args):
    paths, fits = matrix_inputs(args.matrix)
    args.activity.mkdir(parents=True, exist_ok=True)
    evaluate = STEP / "execution/evaluate_loss.py"
    render_metrics = STEP / "figure/script/render_metrics.py"
    render_curves = STEP / "figure/script/render_learning_curves.py"
    options = ["--training", str(args.training), "--bandwidth-training", str(args.bandwidth_training),
               "--output", str(args.report), "--scratch", str(args.scratch)]
    for path in paths.values():
        options.extend(("--run", str(path)))

    if args.command == "build":
        if args.publication.exists():
            raise ValueError("formal publication requires an unused output directory")
        if not args.report.exists():
            command(evaluate, ["score", *options], args.activity, "score")
        revision = args.code_revision
    else:
        revision = axis_run.read_json(args.publication / RECEIPT)["matrix_code_revision"]

    command(evaluate, ["verify", *options], args.activity, "score-verify")
    report = verify_evaluation(args.report, fits)
    bandwidth = axis_run.digest(args.report / "bandwidths.csv")
    if bandwidth != axis_run.digest(args.previous_bandwidths):
        raise ValueError("formal scoring bandwidths differ from the frozen Step 04 bandwidths")
    native = {str(fold): verify_native_control(paths["mse", NATIVE_ALPHA, fold], args.previous_runs / f"fold-{fold}")
              for fold in range(3)}

    with TemporaryDirectory(prefix="step05.1-formal-", dir=args.scratch) as scratch:
        scratch = Path(scratch)
        candidate = scratch / "publication"
        for name, destination in PROJECTIONS.items():
            projected = candidate / destination
            projected.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(args.report / name, projected)
        shutil.copyfile(args.report / "fold-scores.csv", candidate / PLOT_FOLD_SCORES)
        samples = "log/330-resources.jsonl"
        shutil.copyfile(args.matrix_audit / "resource-samples.jsonl", candidate / samples)
        data = candidate / "figure/data"
        figure = candidate / "figure"
        render_args = ["--data", str(data), "--output", str(figure)]
        command(render_metrics, render_args, args.activity, "metric-csv-render")
        command(render_curves, render_args, args.activity, "curve-csv-render")
        replay = scratch / "replay"
        command(render_metrics, ["--data", str(data), "--output", str(replay)],
                args.activity, "metric-csv-replay")
        for png in METRIC_PNGS:
            if axis_run.digest(figure / png) != axis_run.digest(replay / png):
                raise ValueError(f"independent CSV-only renders differ: {png}")
        published = axis_run.file_hashes(
            candidate, (*PROJECTIONS.values(), PLOT_FOLD_SCORES, *(f"figure/{name}" for name in (*METRIC_PNGS, CURVE_PNG)),
                        samples))
        code = (evaluate, render_metrics, render_curves, Path(__file__),
                STEP / "execution/audit_learning_curves.py")
        axis_run.write_json(candidate / RECEIPT, {
            "status": "PASS", "ticket": 330, "matrix_code_revision": revision, "fits": len(fits),
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "report_manifest_sha256": axis_run.digest(args.report / "manifest.json"),
            "report_outputs_sha256": report["outputs"], "frozen_step04_bandwidth_sha256": bandwidth,
            "step04_mse_alpha_1": native, "resources": resources(args.matrix_audit),
            "scientific_replay_tolerances": {"rtol": 0.0, "atol": 0.0},
            "published_files_sha256": published,
            "publication_code_sha256": {path.relative_to(STEP).as_posix(): axis_run.digest(path)
                                        for path in code},
        })
        axis_run.write_json(candidate / CURVE_RECEIPT, audit_curves(candidate, scratch))
        if args.command == "build":
            candidate.rename(args.publication)
        else:
            for path in (*published, RECEIPT, CURVE_RECEIPT):
                if axis_run.digest(candidate / path) != axis_run.digest(args.publication / path):
                    raise ValueError(f"formal publication replay differs: {path}")
    return {"status": "PASS", "fits": len(fits), "receipt_sha256": axis_run.digest(args.publication / RECEIPT)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    for name in ("matrix", "matrix-audit", "report", "publication", "activity", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--previous-runs", type=Path, default=AXIS_STEP / "output/formal-matched-control/runs/denoise-deg-5k")
    parser.add_argument("--previous-bandwidths", type=Path, default=AXIS_STEP / "output/report-matched-control/bandwidths.csv")
    parser.add_argument("--code-revision", help="pinned Git revision used by the matrix")
    args = parser.parse_args()
    if args.command == "build" and not args.code_revision:
        parser.error("build requires --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
