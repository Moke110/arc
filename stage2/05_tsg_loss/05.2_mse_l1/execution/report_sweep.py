"""Verify and publish the formal additive Delta L1 sweep."""

import argparse
from itertools import product
import json
import math
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import pandas as pd

from additive_objective import ALPHAS, LOSS_COMPARISON
from evaluate_sweep import METRICS, REPORT_SCHEMA
from matrix_sweep import SCHEMA, check_registrations, preflight, run_paths
from run_sweep import (DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                       STEP, axis_run)
from matrix import GIB, SAMPLE_SECONDS
from report_formal import ANCHOR_ARRAYS, command
from training_data import TENV_DEFAULT

RECEIPT = Path("log/333-formal.json")
SAMPLES = Path("log/333-resources.jsonl")
COORDINATOR_AUDIT = "coordinator-resource-samples.jsonl"
COORDINATOR_SAMPLES = Path("log/333-coordinator-resources.jsonl")
FIGURES = ("delta-metrics.png", "delta-metrics-improvement.png", "learning-curves.png")
PROJECTIONS = {
    "delta-metrics.csv": "figure/data/delta-metrics.csv",
    "epoch-objectives.csv": "figure/data/epoch-objectives.csv",
    "fold-scores.csv": "log/fold-scores.csv",
    "runs.csv": "log/runs.csv",
    "bandwidths.csv": "log/bandwidths.csv",
}
SCORE_FIELDS = ["fold", "seed", "metric", "value", "n_conditions", "n_eligible", "n_undefined"]
TRACE_FIELDS = ("epoch", "train_loss", "validation_loss", "validation_mse")
NATIVE_COMPONENTS = ("control_mse", "treatment_mse", "delta", "control_kl",
                     "treatment_kl", "prediction_kl")


def matrix_inputs(matrix, training, representations, source):
    manifest = axis_run.read_json(matrix / "manifest.json")
    paths = run_paths(matrix)
    fits = {(fit["alpha"], fit["fold"]): fit for fit in manifest["fits"]}
    preflight_record = matrix / "preflight.json"
    expected, evidence = preflight(training, representations, source)
    if (manifest["schema"] != SCHEMA or len(manifest["fits"]) != len(paths)
            or fits.keys() != paths.keys()
            or manifest["preflight_sha256"] != axis_run.digest(preflight_record)
            or axis_run.read_json(preflight_record) != evidence
            or set(paths.values()) != {path for path in (matrix / "runs").glob("*/*") if path.is_dir()}):
        raise ValueError("formal matrix must bind the 27 registered fits and frozen inputs")
    registration_digests = check_registrations(paths, expected)
    for key, path in paths.items():
        fit = fits[key]
        if (fit["path"] != path.relative_to(matrix).as_posix()
                or fit["registration_sha256"] != registration_digests[key]
                or fit["run_manifest_sha256"] != axis_run.digest(path / "manifest.json")):
            raise ValueError(f"matrix receipt differs for {key}")
    return paths, fits


def verify_native_control(candidate, previous):
    columns = [*TRACE_FIELDS, *(f"{partition}_{component}" for partition in ("train", "validation")
                                 for component in NATIVE_COMPONENTS)]
    rename = {f"{partition}_delta_mse": f"{partition}_delta"
              for partition in ("train", "validation")}
    fresh = pd.read_csv(candidate / "learning_curve.csv", float_precision="round_trip").rename(columns=rename)
    accepted = pd.read_csv(previous / "learning_curve.csv", float_precision="round_trip")
    training, original = (axis_run.read_json(path / "training.json") for path in (candidate, previous))
    epochs = {key: training[key] for key in ("best_epoch", "epochs")}
    if not fresh[columns].equals(accepted[columns]) or any(original[key] != value for key, value in epochs.items()):
        raise ValueError(f"Step 05.1 native objective trace or selected epoch differs: {candidate.name}")
    arrays = axis_run.file_hashes(candidate, ANCHOR_ARRAYS)
    if arrays != axis_run.file_hashes(previous, ANCHOR_ARRAYS):
        raise ValueError(f"Step 05.1 native prediction or Delta arrays differ: {candidate.name}")
    return {**epochs, "native_trace_columns": columns, "arrays_sha256": arrays}


def read_fold_scores(path):
    return pd.read_csv(path, dtype={"alpha": str}, float_precision="round_trip")


def verify_native_scores(current, prior):
    fresh = current.loc[current.alpha == "0", SCORE_FIELDS].reset_index(drop=True)
    accepted = prior.loc[prior.loss.eq("mse") & prior.alpha.eq("1"), SCORE_FIELDS].reset_index(drop=True)
    if len(fresh) != len(METRICS) * 3 or not fresh.equals(accepted):
        raise ValueError("Step 05.1 native fold scores or eligibility counts differ")
    return len(fresh)


def verify_predecessor(previous_step, report, paths):
    receipt_path = previous_step / "log/330-formal.json"
    receipt = axis_run.read_json(receipt_path)
    previous_matrix = previous_step / "output/formal/matrix"
    matrix_record = previous_matrix / "manifest.json"
    if (receipt["status"] != "PASS" or receipt["ticket"] != 330 or receipt["fits"] != 36
            or axis_run.digest(matrix_record) != receipt["matrix_manifest_sha256"]):
        raise ValueError("Step 05.1 native reference must be the accepted formal matrix")
    for filename in ("log/bandwidths.csv", "log/fold-scores.csv"):
        if axis_run.digest(previous_step / filename) != receipt["published_files_sha256"][filename]:
            raise ValueError(f"Step 05.1 published reference differs: {filename}")
    bandwidth = axis_run.digest(report / "bandwidths.csv")
    if bandwidth != receipt["frozen_step04_bandwidth_sha256"]:
        raise ValueError("formal scoring bandwidths differ from the frozen HVG-DEG 5k bandwidths")
    anchors = {fit["fold"]: fit for fit in axis_run.read_json(matrix_record)["fits"]
               if fit["loss"] == "mse" and fit["alpha"] == "1"}
    if set(anchors) != {0, 1, 2}:
        raise ValueError("Step 05.1 native reference must include all three folds")
    matched = {}
    for fold, anchor in sorted(anchors.items()):
        prior = previous_matrix / "runs" / "mse" / "alpha-1" / f"fold-{fold}"
        manifest_path = prior / "manifest.json"
        manifest = axis_run.read_json(manifest_path)
        if (axis_run.digest(manifest_path) != anchor["run_manifest_sha256"]
                or axis_run.digest(prior / "registration.json") != anchor["registration_sha256"]
                or axis_run.digest(prior / "training.json") != manifest["outputs"]["training.json"]):
            raise ValueError(f"Step 05.1 native reference manifest differs: fold {fold}")
        for filename in (*ANCHOR_ARRAYS, "learning_curve.csv"):
            if axis_run.digest(prior / filename) != manifest["outputs"][filename]:
                raise ValueError(f"Step 05.1 native reference file differs: fold {fold}/{filename}")
        matched[str(fold)] = verify_native_control(paths["0", fold], prior)
        matched[str(fold)]["step05_1_run_manifest_sha256"] = anchor["run_manifest_sha256"]
    scores = verify_native_scores(read_fold_scores(report / "fold-scores.csv"),
                                  read_fold_scores(previous_step / "log/fold-scores.csv"))
    return {"receipt_sha256": axis_run.digest(receipt_path), "folds": matched,
            "exact_fold_scores": scores, "frozen_bandwidths_sha256": bandwidth}


def verify_evaluation(report, fits):
    evidence = axis_run.read_json(report / "manifest.json")
    if (evidence["schema"], evidence["status"], evidence["runs"], evidence["smoke"], evidence["metrics"]) != (
            REPORT_SCHEMA, "PASS", len(fits), False, list(METRICS)):
        raise ValueError("formal scoring manifest differs from the formal matrix")
    runs = pd.read_csv(report / "runs.csv", dtype={"alpha": str})
    if list(zip(runs.alpha, runs.fold)) != list(fits):
        raise ValueError("formal scoring report has the wrong fit roster")
    for row in runs.itertuples(index=False):
        fit = fits[row.alpha, row.fold]
        if (row.run_sha256, row.manifest_sha256, row.best_epoch, row.stopping_epoch) != (
                fit["registration_sha256"], fit["run_manifest_sha256"], fit["best_epoch"], fit["stopping_epoch"]):
            raise ValueError(f"scoring/run identity mismatch: {fit['path']}")
    scores = read_fold_scores(report / "fold-scores.csv")
    roster = [(alpha, fold, metric) for (alpha, fold), metric in product(fits, METRICS)]
    if (list(zip(scores.alpha, scores.fold, scores.metric)) != roster
            or not scores.n_conditions.eq(scores.n_eligible + scores.n_undefined).all()
            or (scores[["n_conditions", "n_eligible", "n_undefined"]] < 0).any().any()):
        raise ValueError("formal scoring report must contain 243 ordered fold scores with eligibility counts")
    summary = read_fold_scores(report / "delta-metrics.csv")
    if list(zip(summary.alpha, summary.metric)) != list(product(ALPHAS, METRICS)):
        raise ValueError("formal scoring report must contain all 81 equal-fold aggregates")
    return evidence


def audit_resources(audit):
    samples = [json.loads(line) for line in (audit / "resource-samples.jsonl").read_text().splitlines()]
    gaps = [later["monotonic_seconds"] - earlier["monotonic_seconds"]
            for earlier, later in zip(samples, samples[1:])
            if later["monitor_session"] == earlier["monitor_session"]]
    if not gaps or min(gaps) <= 0 or not any(row["phase"] == "train" for row in samples):
        raise ValueError("formal matrix requires ordered resource samples during fitting")
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
            "maximum_observed_sample_gap_seconds": round(max(gaps), 3),
            "sample_gaps_over_12_seconds": sum(gap > 12 for gap in gaps),
            "minimum_wsl_available_kib": min_ram, "minimum_disk_free_bytes": min_disk,
            "maximum_gpu_used_mib": max(sample["gpu_used_mib"] for sample in samples),
            "minimum_gpu_available_mib": min_gpu, "scope_memory_peak_bytes": scope_peak,
            "phase_audit_sha256": phase_logs}


def coordinator_resource_summary(data):
    samples = [json.loads(line) for line in data.splitlines()]
    times = [sample["seconds_since_start"] for sample in samples]
    gaps = [later - earlier for earlier, later in zip(times, times[1:])]
    if not gaps or any(not math.isfinite(gap) or gap <= 0 for gap in gaps):
        raise ValueError("coordinator resource samples must be in monotonic order")
    return {"sample_interval_seconds": SAMPLE_SECONDS, "samples": len(samples),
            "maximum_observed_sample_gap_seconds": round(max(gaps), 3),
            "sample_gaps_over_12_seconds": sum(gap > 12 for gap in gaps),
            "minimum_wsl_available_kib": min(row["wsl_available_kib"] for row in samples),
            "maximum_gpu_used_mib": max(row["gpu_used_mib"] for row in samples)}


def render(data, figures, activity, attempt):
    args = ["--data", str(data), "--output", str(figures)]
    for filename in ("render_metrics", "render_learning_curves"):
        command(STEP / f"figure/script/{filename}.py", args, activity, f"{attempt}-{filename}")


def publish(args):
    paths, fits = matrix_inputs(args.matrix, args.training, args.representations, args.source)
    args.activity.mkdir(parents=True, exist_ok=True)
    evaluate = STEP / "execution/evaluate_sweep.py"
    options = ["--training", str(args.training), "--bandwidth-training", str(args.bandwidth_training),
               "--output", str(args.report), "--scratch", str(args.scratch)]
    for path in paths.values():
        options.extend(("--run", str(path)))
    if args.command == "build":
        if args.publication.exists():
            raise ValueError("formal publication requires an unused output directory")
        if not args.report.exists():
            command(evaluate, ["score", *options], args.activity, "score")
        matrix_revision = args.code_revision
    else:
        matrix_revision = axis_run.read_json(args.publication / RECEIPT)["matrix_code_revision"]
    command(evaluate, ["verify", *options], args.activity, "score-verify")
    report = verify_evaluation(args.report, fits)
    native = verify_predecessor(args.previous_step, args.report, paths)
    resources = audit_resources(args.matrix_audit)

    with TemporaryDirectory(prefix="step05.2-formal-", dir=args.scratch) as scratch:
        scratch = Path(scratch)
        candidate = scratch / "publication"
        for name, destination in PROJECTIONS.items():
            projected = candidate / destination
            projected.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(args.report / name, projected)
        shutil.copyfile(args.matrix_audit / "resource-samples.jsonl", candidate / SAMPLES)
        sample_paths = [str(SAMPLES)]
        coordinator_audit = args.matrix_audit / COORDINATOR_AUDIT
        if coordinator_audit.exists():
            coordinator_data = coordinator_audit.read_bytes()
            resources["coordinator"] = coordinator_resource_summary(coordinator_data)
            (candidate / COORDINATOR_SAMPLES).write_bytes(coordinator_data)
            sample_paths.append(str(COORDINATOR_SAMPLES))
        data, figures, replay = candidate / "figure/data", candidate / "figure", scratch / "replay"
        render(data, figures, args.activity, "render")
        render(data, replay, args.activity, "independent-render")
        for name in FIGURES:
            if axis_run.digest(figures / name) != axis_run.digest(replay / name):
                raise ValueError(f"independent CSV-only renders differ: {name}")
        published = axis_run.file_hashes(
            candidate, (*PROJECTIONS.values(), *(f"figure/{name}" for name in FIGURES), *sample_paths)
        )
        code = ("execution/additive_objective.py", "execution/evaluate_sweep.py", "execution/run_sweep.py",
                "execution/matrix_sweep.py", "execution/report_sweep.py",
                "figure/script/render_metrics.py", "figure/script/render_learning_curves.py")
        axis_run.write_json(candidate / RECEIPT, {
            "status": "PASS", "ticket": 333, "matrix_code_revision": matrix_revision,
            "fits": len(fits), "fold_scores": len(fits) * len(METRICS),
            "equal_fold_aggregates": len(ALPHAS) * len(METRICS),
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "preflight_sha256": axis_run.digest(args.matrix / "preflight.json"),
            "report_manifest_sha256": axis_run.digest(args.report / "manifest.json"),
            "report_outputs_sha256": report["outputs"], "step05_1_native_mse": native,
            "resources": resources, "scientific_replay_tolerances": {"rtol": 0.0, "atol": 0.0},
            "published_files_sha256": published,
            "publication_code_sha256": axis_run.file_hashes(STEP, code),
        })
        if args.command == "build":
            with axis_run.publication_lock(args.publication) as destination:
                if destination.exists():
                    raise ValueError("formal publication requires an unused output directory")
                candidate.rename(destination)
        else:
            for path in (*published, RECEIPT):
                if axis_run.digest(candidate / path) != axis_run.digest(args.publication / path):
                    raise ValueError(f"formal publication replay differs: {path}")
    return {"status": "PASS", "fits": len(fits), "fold_scores": len(fits) * len(METRICS),
            "receipt_sha256": axis_run.digest(args.publication / RECEIPT)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    for name in ("matrix", "matrix-audit", "report", "publication", "activity", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--previous-step", type=Path, default=LOSS_COMPARISON.parent)
    parser.add_argument("--code-revision", help="pinned Git revision used by the matrix")
    args = parser.parse_args()
    if args.command == "build" and not args.code_revision:
        parser.error("build requires --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
