"""Verify the sealed Step 06.1 matrix and build a CSV-reproducible publication."""

import argparse
import csv
from itertools import product
import json
import math
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import matplotlib
import numpy as np
import pandas as pd
import torch

from evaluate_shape import METRICS, REPORT_SCHEMA
from matrix_shape import SCHEMA as MATRIX_SCHEMA, run_paths
from report_formal import ANCHOR_ARRAYS, command
from report_sweep import SCORE_FIELDS
from run_shape import AXIS, DEFAULT_TRAINING
from shape_model import CONTROL, LOSS, SHAPES, STEP, axis_run
from training_data import TENV_DEFAULT

RECEIPT = "log/347-formal.json"
FIGURE_RECEIPT = "log/figure-render.json"
REFERENCE_RECEIPT = "log/333-formal.json"
REFERENCE_SHA256 = "bdc65b62eafb834cbe1cc498a8d63c7ab00bbe6fc1d9bd545f5f0ad2a606bcf8"
RENDERER = STEP / "figure/script/render_shape.py"
FIGURES = ("delta-metrics.png", "delta-metrics-paired-difference.png", "learning-curves.png")
EPOCH_OBJECTIVES = "figure/data/epoch-objectives.csv"
FIGURE_SOURCES = ("figure/data/delta-metrics.csv", "figure/data/fold-scores.csv", EPOCH_OBJECTIVES)
SCIENTIFIC_CODE = ("execution/matrix_shape.py", "execution/run_shape.py",
                   "execution/shape_model.py", "execution/evaluate_shape.py")
PROJECTIONS = {
    "runs.csv": "log/runs.csv",
    "bandwidths.csv": "log/bandwidths.csv",
    "kl-diagnostics.csv": "log/kl-diagnostics.csv",
    "fold-scores.csv": "figure/data/fold-scores.csv",
    "delta-metrics.csv": "figure/data/delta-metrics.csv",
}
RESOURCE_LOG = "log/347-resources.jsonl"
COORDINATOR_LOG = "log/347-coordinator-resources.jsonl"
GIB = 1024**3


def verify_checkpoint(candidate, reference, candidate_registration, reference_registration):
    fresh = torch.load(candidate, map_location="cpu", weights_only=True)
    accepted = torch.load(reference, map_location="cpu", weights_only=True)
    if (fresh["registration_sha256"], accepted["registration_sha256"]) != (
            candidate_registration, reference_registration):
        raise ValueError("selected checkpoint registration differs from its run")
    for field in ("epoch", "validation_loss", "validation_mse"):
        if fresh[field] != accepted[field]:
            raise ValueError(f"selected checkpoint {field} differs")
    states = fresh["state_dict"], accepted["state_dict"]
    if states[0].keys() != states[1].keys():
        raise ValueError("selected checkpoint parameter names differ")
    for name, value in states[0].items():
        prior = states[1][name]
        if value.dtype != prior.dtype or not torch.equal(value, prior):
            raise ValueError(f"selected checkpoint tensor differs: {name}")


def verify_matrix(matrix):
    paths = run_paths(matrix)
    receipt = axis_run.read_json(matrix / "manifest.json")
    if (receipt["schema"] != MATRIX_SCHEMA
            or receipt["preflight_sha256"] != axis_run.digest(matrix / "preflight.json")
            or [(fit["shape"], fit["fold"]) for fit in receipt["fits"]] != list(paths)):
        raise ValueError("formal matrix receipt has an invalid preflight or coordinate roster")
    fits = dict(zip(paths, receipt["fits"], strict=True))
    for key, path in paths.items():
        fit = fits[key]
        if (fit["path"] != path.relative_to(matrix).as_posix()
                or fit["run_manifest_sha256"] != axis_run.digest(path / "manifest.json")
                or fit["registration_sha256"] != axis_run.digest(path / "registration.json")):
            raise ValueError(f"sealed matrix run differs: {key}")
    return paths, fits


def verify_evaluation(report, fits):
    receipt = axis_run.read_json(report / "manifest.json")
    if (receipt["schema"], receipt["status"], receipt["axis"], receipt["smoke"], receipt["runs"],
            receipt["fold_reference"], receipt["metrics"], receipt["fold_scores"]) != (
            REPORT_SCHEMA, "PASS", AXIS, False, 24, CONTROL, list(METRICS), 216):
        raise ValueError("formal scoring receipt differs from the held-out shape comparison")
    axis_run.check_outputs(report, receipt)
    runs = pd.read_csv(report / "runs.csv")
    if list(zip(runs["shape"], runs.fold)) != list(fits):
        raise ValueError("formal scoring run roster differs from the matrix")
    for row in runs.itertuples(index=False):
        fit = fits[row.shape, row.fold]
        if (row.run_sha256, row.manifest_sha256, row.best_epoch, row.stopping_epoch) != (
                fit["registration_sha256"], fit["run_manifest_sha256"], fit["best_epoch"], fit["stopping_epoch"]):
            raise ValueError(f"formal scoring receipt differs for {fit['path']}")
    scores = pd.read_csv(report / "fold-scores.csv", float_precision="round_trip")
    eligibility = scores[["n_conditions", "n_eligible", "n_undefined"]]
    if (list(zip(scores["shape"], scores.fold, scores.metric)) !=
            [(shape, fold, metric) for shape, fold in fits for metric in METRICS]
            or (eligibility < 0).any().any()
            or not scores.n_conditions.eq(scores.n_eligible + scores.n_undefined).all()
            or not scores.value.notna().eq(scores.n_eligible.gt(0)).all()):
        raise ValueError("formal scores require all 216 fold/metric coordinates and valid eligibility")
    summary = pd.read_csv(report / "delta-metrics.csv", float_precision="round_trip")
    if list(zip(summary["shape"], summary.metric)) != list(product(SHAPES, METRICS)):
        raise ValueError("formal scoring requires all 72 equal-fold aggregates")
    return receipt, scores


def verify_reference_config(candidate, reference):
    fields = ("model", "environment", "objective", "scenario", "seed", "validation_seed", "test_seed",
              "axis_sha256", "training_sha256", "raw_sha256", "folds", "native_sha256",
              "kpgt_manifest_sha256", "kpgt_vectors_sha256", "conditions", "partitions", "patience", "max_epochs")
    for field in fields:
        if candidate[field] != reference[field]:
            raise ValueError(f"Step 05.2 reference protocol or runtime differs: {field}")


def verify_control(candidate, prior, fit):
    active = axis_run.read_json(candidate / "registration.json")
    reference = axis_run.read_json(prior / "registration.json")
    verify_reference_config(active, reference)
    if axis_run.digest(prior / "registration.json") != fit["registration_sha256"]:
        raise ValueError("Step 05.2 reference registration differs from the sealed matrix")
    run_receipt = prior / "manifest.json"
    if axis_run.digest(run_receipt) != fit["run_manifest_sha256"]:
        raise ValueError("Step 05.2 run receipt differs from the sealed matrix")
    recorded = axis_run.read_json(run_receipt)
    if (recorded["registration_sha256"] != fit["registration_sha256"]
            or [recorded[name] for name in ("validation_replay", "validation_objective_replay", "test_replay")]
            != ["exact"] * 3):
        raise ValueError("Step 05.2 reference run requires its exact replay")
    axis_run.check_outputs(prior, recorded)
    training, original = (axis_run.read_json(path / "training.json") for path in (candidate, prior))
    fields = ("best_epoch", "epochs", "validation_loss", "validation_mse", "parameters")
    if any(training[key] != original[key] for key in fields):
        raise ValueError("Step 05.2 trace selected a different checkpoint or objective")
    traces = [pd.read_csv(path / "learning_curve.csv", float_precision="round_trip") for path in (candidate, prior)]
    if not traces[0].equals(traces[1]):
        raise ValueError("Step 05.2 training or validation trace differs")
    arrays = axis_run.file_hashes(candidate, ANCHOR_ARRAYS)
    if any(value != recorded["outputs"][name] for name, value in arrays.items()):
        raise ValueError("Step 05.2 validation, test or Delta prediction arrays differ")
    verify_checkpoint(candidate / "checkpoint.pt", prior / "checkpoint.pt",
                      axis_run.digest(candidate / "registration.json"), fit["registration_sha256"])
    return {"best_epoch": training["best_epoch"], "stopping_epoch": training["epochs"],
            "trace_exact": True, "selected_checkpoint_state_exact": True, "arrays_sha256": arrays,
            "reference_run_manifest_sha256": fit["run_manifest_sha256"]}


def verify_predecessor(previous, report, scores, paths):
    receipt_path = previous / REFERENCE_RECEIPT
    if axis_run.digest(receipt_path) != REFERENCE_SHA256:
        raise ValueError("Step 05.2 publication evidence differs from the sealed reference")
    receipt = axis_run.read_json(receipt_path)
    matrix = previous / "output/formal/matrix"
    if axis_run.digest(matrix / "manifest.json") != receipt["matrix_manifest_sha256"]:
        raise ValueError("Step 05.2 matrix differs from the sealed reference")
    frozen = "log/bandwidths.csv"
    accepted_scores = "log/fold-scores.csv"
    for name in (frozen, accepted_scores):
        if axis_run.digest(previous / name) != receipt["published_files_sha256"][name]:
            raise ValueError(f"Step 05.2 reference table differs: {name}")
    bandwidth = receipt["published_files_sha256"][frozen]
    if axis_run.digest(report / "bandwidths.csv") != bandwidth:
        raise ValueError("Step 06.1 scoring bandwidths differ from the frozen Step 05.2 values")
    prior_scores = pd.read_csv(previous / accepted_scores, dtype={"alpha": str}, float_precision="round_trip")
    current = scores.loc[scores["shape"] == CONTROL, SCORE_FIELDS].reset_index(drop=True)
    accepted = prior_scores.loc[prior_scores.alpha.eq("1.75"), SCORE_FIELDS].reset_index(drop=True)
    if len(current) != 27 or not current.equals(accepted):
        raise ValueError("Step 05.2 weight-1.75 scores or eligibility differ")
    fits = [fit for fit in axis_run.read_json(matrix / "manifest.json")["fits"] if fit["alpha"] == "1.75"]
    if [fit["fold"] for fit in fits] != [0, 1, 2]:
        raise ValueError("Step 05.2 weight-1.75 must cover the three frozen folds")
    parity = {str(fit["fold"]): verify_control(paths[CONTROL, fit["fold"]], matrix / fit["path"], fit)
              for fit in fits}
    return {"receipt_sha256": REFERENCE_SHA256, "matrix_code_revision": receipt["matrix_code_revision"],
            "exact_fold_scores": len(current), "bandwidths_sha256": bandwidth, "folds": parity,
            "tolerance": {"rtol": 0, "atol": 0}}


def read_samples(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def gaps_summary(samples):
    gaps = [later["monotonic_seconds"] - earlier["monotonic_seconds"]
            for earlier, later in zip(samples, samples[1:])
            if earlier.get("monitor_session") == later.get("monitor_session")]
    if not gaps or min(gaps) <= 0 or not np.isfinite(gaps).all():
        raise ValueError("resource audit must record monotonic samples")
    return {"samples": len(samples), "sample_interval_seconds": 10,
            "maximum_gap_seconds": round(max(gaps), 3), "gaps_over_12_seconds": sum(gap > 12 for gap in gaps)}


def audit_resources(audit, coordinator):
    samples = read_samples(audit / "resource-samples.jsonl")
    if not any(row["phase"] == "train" for row in samples):
        raise ValueError("resource audit is missing samples during training")
    scope_files = sorted(audit.glob("scope-*.json"))
    peak = max(axis_run.read_json(path)["memory_peak_bytes"] for path in scope_files)
    ram = min(row["wsl_available_kib"] for row in samples)
    gpu = min(row["gpu_total_mib"] - row["gpu_used_mib"] for row in samples)
    disk = min(row["disk_free_bytes"] for row in samples)
    if ram < 12 * GIB // 1024 or gpu < 512 or disk < GIB or peak > 10 * GIB:
        raise ValueError("formal matrix exceeded its resource caps or headroom")
    files = list(scope_files)
    for path in run_paths(audit).values():
        for phase in ("register", "train", "predict", "replay"):
            measurements = sorted(path.glob(f"{phase}-*.resources"))
            if not measurements:
                raise ValueError(f"formal matrix requires resource evidence: {path}/{phase}")
            for measurement in measurements:
                files.extend((measurement, measurement.with_suffix(".log")))
    coordination = read_samples(coordinator)
    if (coordination[0]["monotonic_seconds"] > samples[0]["monotonic_seconds"]
            or coordination[-1]["monotonic_seconds"] < samples[-1]["monotonic_seconds"] - 12):
        raise ValueError("coordinator resource samples must cover the formal matrix")
    return {**gaps_summary(samples), "peak_scope_memory_bytes": peak,
            "minimum_wsl_available_kib": ram, "minimum_gpu_available_mib": gpu,
            "minimum_disk_free_bytes": disk,
            "phase_audit_sha256": axis_run.file_hashes(audit, [file.relative_to(audit).as_posix() for file in files]),
            "coordinator": {**gaps_summary(coordination),
                            "minimum_wsl_available_kib": min(row["wsl"]["MemAvailable_kib"] for row in coordination),
                            "maximum_sampled_gpu_used_mib": max(int(row["gpu_total_used_mib_and_utilization_percent"].split(",")[1])
                                                               for row in coordination
                                                               if isinstance(row["gpu_total_used_mib_and_utilization_percent"], str))}}


def write_epoch_objectives(destination, paths, fits):
    """Aggregate observed objectives over the full cohort at every plotted epoch."""
    hashes = {}
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("shape", "epoch", "mean_train_loss", "mean_validation_loss"))
        for shape in SHAPES:
            traces = []
            for fold in range(3):
                fit = fits[shape, fold]
                run = paths[shape, fold]
                trace_path = run / "learning_curve.csv"
                digest = axis_run.digest(trace_path)
                if digest != axis_run.read_json(run / "manifest.json")["outputs"][trace_path.name]:
                    raise ValueError(f"sealed learning curve differs: {fit['path']}")
                hashes[f"{fit['path']}/{trace_path.name}"] = digest
                with trace_path.open(newline="", encoding="utf-8") as trace:
                    rows = [(int(row["epoch"]), float(row["train_loss"]), float(row["validation_loss"]))
                            for row in csv.DictReader(trace)]
                if ([epoch for epoch, _, _ in rows] != list(range(1, fit["stopping_epoch"] + 1))
                        or any(not (math.isfinite(train) and math.isfinite(validation)
                                    and train > 0 and validation > 0) for _, train, validation in rows)):
                    raise ValueError(f"invalid observed objective trace: {fit['path']}")
                traces.append(rows)
            for aligned in zip(*traces):
                epoch = aligned[0][0]
                mean_train = math.fsum(row[1] for row in aligned) / len(aligned)
                mean_validation = math.fsum(row[2] for row in aligned) / len(aligned)
                writer.writerow((shape, epoch, mean_train, mean_validation))
    return hashes


def figure_evidence(candidate, trace_hashes):
    stage = STEP.parents[1]
    dependencies = (stage / "05_tsg_loss/figure/script/loss_plot_style.py",
                    stage / "util/output/api/stage2_util/palette.py",
                    stage / "util/output/api/stage2_util/plotting.py")
    return {
        "schema": "stage2-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": axis_run.platform.python_version(),
                        "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "independent_renders": 2,
        "scientific_receipt_sha256": axis_run.digest(candidate / RECEIPT),
        "learning_curve_source_sha256": trace_hashes,
        "csv_files_sha256": axis_run.file_hashes(candidate, FIGURE_SOURCES),
        "png_files_sha256": axis_run.file_hashes(candidate, [f"figure/{name}" for name in FIGURES]),
        "renderer_files_sha256": axis_run.file_hashes(stage, [
            path.relative_to(stage).as_posix() for path in (RENDERER, *dependencies)]),
        "publication_builder_sha256": axis_run.digest(Path(__file__)),
    }


def verify_scientific_receipt(candidate, evidence):
    """Compare the scientific fields of the original or current publication receipt."""
    recorded = axis_run.read_json(candidate / RECEIPT)
    if "rendering_environment" in recorded:
        # The original combined receipt seals now-historical PNGs and replay code.
        recorded.pop("rendering_environment")
        recorded.pop("independent_csv_render")
        recorded["outputs"] = {name: digest for name, digest in recorded["outputs"].items()
                               if not (name.startswith("figure/") and name.endswith(".png"))}
        recorded["publication_code_sha256"] = {
            name: digest for name, digest in recorded["publication_code_sha256"].items()
            if name not in ("execution/report_shape.py", "figure/script/render_shape.py")}
    if recorded != evidence:
        raise ValueError("formal scientific receipt differs from the verified replay")
    axis_run.check_outputs(candidate, recorded)


def publish(args):
    args.activity.mkdir(parents=True, exist_ok=True)
    paths, fits = verify_matrix(args.matrix)
    evaluate = STEP / "execution/evaluate_shape.py"
    options = ["--training", str(args.training), "--bandwidth-training", str(args.bandwidth_training),
               "--output", str(args.report), "--scratch", str(args.scratch)]
    for path in paths.values():
        options.extend(("--run", str(path)))
    with (axis_run.publication_lock(args.publication) as publication,
          TemporaryDirectory(prefix="step06.1-formal-", dir=args.scratch) as scratch):
        candidate = Path(scratch) / "publication"
        if args.command == "build":
            if publication.exists():
                raise ValueError("formal publication requires an unused output directory")
            if not args.report.exists():
                command(evaluate, ["score", *options], args.activity, "score")
            coordinator = args.coordinator_samples
            revision = args.code_revision
        else:
            coordinator = publication / COORDINATOR_LOG
            revision = axis_run.read_json(publication / RECEIPT)["matrix_code_revision"]
        command(evaluate, ["verify", *options], args.activity, "score-verify")
        report_receipt, scores = verify_evaluation(args.report, fits)
        parity = verify_predecessor(args.previous_step, args.report, scores, paths)

        projections = {destination: args.report / source for source, destination in PROJECTIONS.items()}
        projections.update({RESOURCE_LOG: args.matrix_audit / "resource-samples.jsonl",
                            COORDINATOR_LOG: coordinator})
        for destination, source in projections.items():
            target = candidate / destination
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        resources = audit_resources(args.matrix_audit, candidate / COORDINATOR_LOG)
        trace_hashes = write_epoch_objectives(candidate / EPOCH_OBJECTIVES, paths, fits)
        for directory in (candidate / "figure", Path(scratch) / "repeat-render"):
            command(RENDERER, ["--data", str(candidate / "figure/data"), "--output", str(directory)],
                    args.activity, directory.name)
        for name in FIGURES:
            if axis_run.digest(candidate / "figure" / name) != axis_run.digest(Path(scratch) / "repeat-render" / name):
                raise ValueError(f"CSV-only figure replay differs: {name}")
        evidence = {
            "status": "PASS", "ticket": 347, "matrix_code_revision": revision,
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "preflight_sha256": axis_run.digest(args.matrix / "preflight.json"),
            "report_manifest_sha256": axis_run.digest(args.report / "manifest.json"),
            "report_outputs_sha256": report_receipt["outputs"],
            "fits": len(fits), "independent_exact_replays": len(fits), "fold_scores": len(scores),
            "equal_fold_aggregates": len(SHAPES) * len(METRICS),
            "reference_parity": parity, "resources": resources,
            "publication_code_sha256": axis_run.file_hashes(STEP, SCIENTIFIC_CODE),
            "outputs": axis_run.file_hashes(candidate, projections),
        }
        if args.command == "build":
            axis_run.write_json(candidate / RECEIPT, evidence)
            axis_run.write_json(candidate / FIGURE_RECEIPT, figure_evidence(candidate, trace_hashes))
            if publication.exists():
                raise ValueError("publication destination appeared during acceptance")
            candidate.rename(publication)
        else:
            verify_scientific_receipt(publication, evidence)
            shutil.copyfile(publication / RECEIPT, candidate / RECEIPT)
            renders = figure_evidence(candidate, trace_hashes)
            if renders != axis_run.read_json(publication / FIGURE_RECEIPT):
                raise ValueError("figure receipt differs from the verified replay")
            for field in ("csv_files_sha256", "png_files_sha256"):
                axis_run.check_outputs(publication, {"outputs": renders[field]})
    return {"status": "PASS", "fits": len(fits), "fold_scores": len(scores),
            "receipt_sha256": axis_run.digest(args.publication / RECEIPT),
            "figure_receipt_sha256": axis_run.digest(args.publication / FIGURE_RECEIPT)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--previous-step", type=Path, default=LOSS / "05.2_mse_l1")
    for name in ("matrix", "matrix-audit", "report", "publication", "activity", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--coordinator-samples", type=Path)
    parser.add_argument("--code-revision")
    args = parser.parse_args()
    if args.command == "build" and (args.coordinator_samples is None or not args.code_revision):
        parser.error("build requires --coordinator-samples and --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
