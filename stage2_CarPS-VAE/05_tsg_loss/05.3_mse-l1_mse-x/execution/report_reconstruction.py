"""Verify and publish the formal reconstruction MSE sweep."""

import argparse
from itertools import product
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import pandas as pd

from evaluate_reconstruction import METRICS, REPORT_SCHEMA
from matrix_reconstruction import SCHEMA, check_registrations, preflight, run_paths
from run_reconstruction import (DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                                PREDECESSOR_EXECUTION, STEP, axis_run)
from sweep_protocol import ALPHAS, DELTA_L1_WEIGHT, REFERENCE_ALPHA
from additive_objective import ALPHA_LABELS as ADDITIVE_WEIGHT_LABELS
from evaluate_sweep import read_run as read_predecessor_run
from matrix import GIB, SAMPLE_SECONDS
from report_formal import ANCHOR_ARRAYS, command
from report_sweep import coordinator_resource_summary, read_fold_scores
from training_data import TENV_DEFAULT

RECEIPT = Path("log/336-formal.json")
SAMPLES = Path("log/336-resources.jsonl")
COORDINATOR_AUDIT = "coordinator-resource-samples.jsonl"
COORDINATOR_SAMPLES = Path("log/336-coordinator-resources.jsonl")
FIGURES = ("delta-metrics.png", "delta-metrics-improvement.png", "learning-curves.png")
PLOT_FOLD_SCORES = "figure/data/fold-scores.csv"
PROJECTIONS = {
    "delta-metrics.csv": "figure/data/delta-metrics.csv",
    "epoch-objectives.csv": "figure/data/epoch-objectives.csv",
    "fold-scores.csv": "log/fold-scores.csv",
    "runs.csv": "log/runs.csv",
    "bandwidths.csv": "log/bandwidths.csv",
}
SCORE_FIELDS = ["fold", "seed", "metric", "value", "n_conditions", "n_eligible", "n_undefined"]
PREDECESSOR_ALPHA = ADDITIVE_WEIGHT_LABELS[DELTA_L1_WEIGHT]


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


def verify_reference_scores(current, prior):
    fresh = current.loc[current.alpha.eq(REFERENCE_ALPHA), SCORE_FIELDS].reset_index(drop=True)
    accepted = prior.loc[prior.alpha.eq(PREDECESSOR_ALPHA), SCORE_FIELDS].reset_index(drop=True)
    expected = list(product(range(3), METRICS))
    if list(zip(fresh.fold, fresh.metric)) != expected or not fresh.equals(accepted):
        raise ValueError("Step 05.2 reference fold scores or eligibility counts differ")
    return len(fresh)


def verify_predecessor(previous_step, report, paths):
    receipt_path = previous_step / "log/333-formal.json"
    receipt = axis_run.read_json(receipt_path)
    previous_matrix = previous_step / "output/formal/matrix"
    matrix_record = previous_matrix / "manifest.json"
    matrix_digest = axis_run.digest(matrix_record)
    if (receipt["status"] != "PASS" or receipt["ticket"] != 333 or receipt["fits"] != 27
            or matrix_digest != receipt["matrix_manifest_sha256"]):
        raise ValueError("Step 05.2 reference must come from the accepted formal matrix")
    bandwidth_path, fold_scores_path = "log/bandwidths.csv", "log/fold-scores.csv"
    for name in (bandwidth_path, fold_scores_path):
        if axis_run.digest(previous_step / name) != receipt["published_files_sha256"][name]:
            raise ValueError(f"Step 05.2 published reference differs: {name}")
    bandwidth = axis_run.digest(report / "bandwidths.csv")
    if (bandwidth != receipt["published_files_sha256"][bandwidth_path]
            or bandwidth != receipt["step05_1_native_mse"]["frozen_bandwidths_sha256"]):
        raise ValueError("formal scoring bandwidths differ from the frozen HVG-DEG 5k bandwidths")
    anchors = [fit for fit in axis_run.read_json(matrix_record)["fits"]
               if fit["alpha"] == PREDECESSOR_ALPHA]
    if sorted(fit["fold"] for fit in anchors) != [0, 1, 2]:
        raise ValueError("Step 05.2 reference must include all three folds")

    matched = {}
    for anchor in anchors:
        fold = anchor["fold"]
        relative = f"runs/alpha-{PREDECESSOR_ALPHA}/fold-{fold}"
        prior = previous_matrix / relative
        if (anchor["path"] != relative
                or axis_run.digest(prior / "manifest.json") != anchor["run_manifest_sha256"]
                or axis_run.digest(prior / "registration.json") != anchor["registration_sha256"]):
            raise ValueError(f"Step 05.2 reference manifest differs: fold {fold}")
        key, (_, prior_config, prior_manifest, prior_training, _) = read_predecessor_run(prior, False)
        if key != (PREDECESSOR_ALPHA, fold):
            raise ValueError(f"Step 05.2 reference registration differs: fold {fold}")
        fresh_path = paths[REFERENCE_ALPHA, fold]
        fresh_config = axis_run.read_json(fresh_path / "registration.json")
        for field in ("conditions", "folds", "partitions", "seed", "max_epochs", "patience",
                      "training_sha256", "kpgt_vectors_sha256", "native_sha256", "environment",
                      "test_seed", "validation_seed"):
            if fresh_config[field] != prior_config[field]:
                raise ValueError(f"Step 05.2 reference {field} differs: fold {fold}")
        validated_models = ({field: value for field, value in config["model"].items()
                             if field != "validation_objective"}
                            for config in (fresh_config, prior_config))
        fresh_model, prior_model = validated_models
        if fresh_model != prior_model:
            raise ValueError(f"Step 05.2 reference model or training settings differ: fold {fold}")
        arrays_and_trace = (*ANCHOR_ARRAYS, "learning_curve.csv")
        artifacts = {name: prior_manifest["outputs"][name] for name in arrays_and_trace}
        if axis_run.file_hashes(fresh_path, arrays_and_trace) != artifacts:
            raise ValueError(f"Step 05.2 reference loss trace, predictions or Delta arrays differ: fold {fold}")
        fresh_training = axis_run.read_json(fresh_path / "training.json")
        epochs = {name: fresh_training[name] for name in ("best_epoch", "epochs")}
        if any(prior_training[name] != value for name, value in epochs.items()):
            raise ValueError(f"Step 05.2 reference selected or stopping epoch differs: fold {fold}")
        matched[str(fold)] = {
            **epochs, "exact_files_sha256": artifacts,
            "step05_2_run_manifest_sha256": anchor["run_manifest_sha256"],
            "step05_2_registration_sha256": anchor["registration_sha256"],
        }
    scores = verify_reference_scores(read_fold_scores(report / "fold-scores.csv"),
                                     read_fold_scores(previous_step / fold_scores_path))
    return {"receipt_sha256": axis_run.digest(receipt_path), "matrix_manifest_sha256": matrix_digest,
            "delta_l1_weight": PREDECESSOR_ALPHA, "folds": matched, "exact_fold_scores": scores,
            "frozen_bandwidths_sha256": bandwidth,
            "fold_scores_sha256": receipt["published_files_sha256"][fold_scores_path]}


def verify_evaluation(report, fits):
    evidence = axis_run.read_json(report / "manifest.json")
    if (evidence["schema"], evidence["status"], evidence["runs"], evidence["smoke"], evidence["metrics"],
            evidence["reference_alpha"], evidence["ordered_populations_verified"],
            evidence["selected_checkpoint_and_replay_verified"]) != (
                REPORT_SCHEMA, "PASS", len(fits), False, list(METRICS), REFERENCE_ALPHA, True, True):
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
    if not summary.loc[summary.alpha.eq(REFERENCE_ALPHA), "score_difference"].eq(0).all():
        raise ValueError("the nine unit-weight reference score differences must be exactly zero")
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
    phase_logs, phase_peak_rss = {}, {}
    for path in run_paths(audit).values():
        for phase in ("register", "train", "predict", "replay"):
            measurements = sorted(path.glob(f"{phase}-*.resources"))
            peaks = [int(line.rsplit(":", 1)[1]) for measurement in measurements
                     for line in measurement.read_text().splitlines()
                     if line.lstrip().startswith("Maximum resident set size (kbytes):")]
            if not peaks or max(peaks) <= 0:
                raise ValueError(f"formal resource audit lacks {phase} measurements: {path}")
            phase_peak_rss[phase] = max(phase_peak_rss.get(phase, 0), *peaks)
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
            "phase_peak_rss_kib": phase_peak_rss, "phase_audit_sha256": phase_logs}


def render(data, figures, activity, attempt):
    args = ["--data", str(data), "--output", str(figures)]
    for filename in ("render_metrics", "render_learning_curves"):
        command(STEP / f"figure/script/{filename}.py", args, activity, f"{attempt}-{filename}")


def publish(args):
    paths, fits = matrix_inputs(args.matrix, args.training, args.representations, args.source)
    args.activity.mkdir(parents=True, exist_ok=True)
    evaluate = STEP / "execution/evaluate_reconstruction.py"
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
    reference = verify_predecessor(args.previous_step, args.report, paths)
    resources = audit_resources(args.matrix_audit)

    with TemporaryDirectory(prefix="step05.3-formal-", dir=args.scratch) as scratch:
        scratch = Path(scratch)
        candidate = scratch / "publication"
        for name, destination in PROJECTIONS.items():
            projected = candidate / destination
            projected.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(args.report / name, projected)
        shutil.copyfile(args.report / "fold-scores.csv", candidate / PLOT_FOLD_SCORES)
        shutil.copyfile(args.matrix_audit / "resource-samples.jsonl", candidate / SAMPLES)
        sample_paths = [str(SAMPLES)]
        coordinator_audit = args.matrix_audit / COORDINATOR_AUDIT
        if coordinator_audit.exists():
            coordinator_data = coordinator_audit.read_bytes()
            resources["coordinator"] = coordinator_resource_summary(coordinator_data)
            if resources["coordinator"]["minimum_wsl_available_kib"] < 12 * GIB // 1024:
                raise ValueError("the coordinator observed less than 12 GiB of WSL headroom")
            (candidate / COORDINATOR_SAMPLES).write_bytes(coordinator_data)
            sample_paths.append(str(COORDINATOR_SAMPLES))
        data, figures, replay = candidate / "figure/data", candidate / "figure", scratch / "replay"
        render(data, figures, args.activity, "render")
        render(data, replay, args.activity, "independent-render")
        for name in FIGURES:
            if axis_run.digest(figures / name) != axis_run.digest(replay / name):
                raise ValueError(f"independent CSV-only renders differ: {name}")
        published = axis_run.file_hashes(
            candidate, (*PROJECTIONS.values(), PLOT_FOLD_SCORES, *(f"figure/{name}" for name in FIGURES), *sample_paths)
        )
        code = ("execution/reconstruction_objective.py", "execution/evaluate_reconstruction.py",
                "execution/run_reconstruction.py", "execution/sweep_protocol.py",
                "execution/matrix_reconstruction.py", "execution/report_reconstruction.py",
                "figure/script/render_metrics.py", "figure/script/render_learning_curves.py")
        axis_run.write_json(candidate / RECEIPT, {
            "status": "PASS", "ticket": 336, "matrix_code_revision": matrix_revision,
            "fits": len(fits), "fold_scores": len(fits) * len(METRICS),
            "equal_fold_aggregates": len(ALPHAS) * len(METRICS),
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "preflight_sha256": axis_run.digest(args.matrix / "preflight.json"),
            "report_manifest_sha256": axis_run.digest(args.report / "manifest.json"),
            "report_outputs_sha256": report["outputs"], "step05_2_reference": reference,
            "resources": resources, "reference_alpha": REFERENCE_ALPHA,
            "csv_only_render_replay": "byte-identical",
            "scientific_replay_tolerances": {"rtol": 0.0, "atol": 0.0},
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
    parser.add_argument("--previous-step", type=Path, default=PREDECESSOR_EXECUTION.parent)
    parser.add_argument("--code-revision", help="pinned Git revision used by the matrix")
    args = parser.parse_args()
    if args.command == "build" and not args.code_revision:
        parser.error("build requires --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
