"""Verify the formal Step 06.2 matrix and publish a reproducible comparison."""

import argparse
from itertools import product
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import matplotlib
import numpy as np
import pandas as pd

from architecture_model import STEP, axis_run
from architecture_roster import CONTROL, SHAPES, coordinates
from evaluate_architecture import METRICS, REPORT_SCHEMA, score_or_verify
from matrix_architecture import SCHEMA as MATRIX_SCHEMA, check_registrations, fit_coordinates, preflight
from prior_architecture import (PRIOR_MATRIX, PRIOR_STEP, PUBLISHED_BANDWIDTHS, PUBLISHED_SCORES,
                                load_originals, verify_original_scores)
from prior_shape import verify_matrix
from report_formal import command
from report_shape import gaps_summary, read_samples
from run_architecture import AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING, PHASES
from training_data import TENV_DEFAULT

RECEIPT = "log/353-formal.json"
FIGURE_RECEIPT = "log/figure-render.json"
MATRIX_SAMPLES = "log/353-resources.jsonl"
COORDINATOR_SAMPLES = "log/353-coordinator-resources.jsonl"
STAGE = STEP.parents[1]
RENDERER = STEP / "figure/script/render_architecture.py"
FIGURES = ("delta-metrics.png", "delta-metrics-paired-difference.png", "learning-curves.png")
FIGURE_TABLES = ("delta-metrics.csv", "fold-scores.csv", "epoch-objectives.csv")
PROJECTIONS = {**{name: f"figure/data/{name}" for name in FIGURE_TABLES},
               "bandwidths.csv": "log/bandwidths.csv", "runs.csv": "log/runs.csv"}
GIB = 1024 ** 3


def read_matrix(matrix):
    return verify_matrix(matrix, fit_coordinates(False), MATRIX_SCHEMA)


def verify_formal_configurations(args, paths):
    expected, evidence = preflight(args.training, args.representations, args.source, False)
    if axis_run.read_json(args.matrix / "preflight.json") != evidence:
        raise ValueError("Step 06.2 preflight differs from the current frozen inputs and run code")
    check_registrations(paths, expected, False)


def verify_report(report, fits):
    receipt = axis_run.read_json(report / "manifest.json")
    if (receipt["schema"], receipt["status"], receipt["axis"], receipt["smoke"], receipt["runs"],
            receipt["new_runs"], receipt["reused_runs"], receipt["fold_reference"], receipt["metrics"],
            receipt["fold_scores"]) != (
            REPORT_SCHEMA, "PASS", AXIS, False, 24, 21, 3, CONTROL, list(METRICS), 216):
        raise ValueError("formal scoring receipt differs from the Step 06.2 comparison")
    axis_run.check_outputs(report, receipt)

    runs = pd.read_csv(report / "runs.csv")
    if list(zip(runs["shape"], runs["fold"])) != list(coordinates()):
        raise ValueError("formal scoring roster differs from the Step 06.2 matrix")
    for row in runs.itertuples(index=False):
        fit = fits[row.shape, row.fold]
        if (row.run_sha256, row.manifest_sha256, row.best_epoch, row.stopping_epoch) != (
                fit["registration_sha256"], fit["run_manifest_sha256"], fit["best_epoch"], fit["stopping_epoch"]):
            raise ValueError(f"formal scoring fit identity differs: {row.shape}, fold {row.fold}")

    scores = pd.read_csv(report / "fold-scores.csv", float_precision="round_trip")
    eligibility = scores[["n_conditions", "n_eligible", "n_undefined"]]
    if (list(zip(scores["shape"], scores["fold"], scores["metric"])) !=
            [(shape, fold, metric) for shape, fold in coordinates() for metric in METRICS]
            or (eligibility < 0).any().any() or eligibility.isna().any().any()
            or not scores.n_conditions.eq(scores.n_eligible + scores.n_undefined).all()
            or not scores.value.notna().eq(scores.n_eligible.gt(0)).all()
            or not np.isfinite(scores.value.dropna()).all()):
        raise ValueError("formal scores require the full fold/metric roster and valid eligibility")
    summary = pd.read_csv(report / "delta-metrics.csv", float_precision="round_trip")
    if list(zip(summary["shape"], summary["metric"])) != list(product(SHAPES, METRICS)):
        raise ValueError("formal scoring requires all 72 equal-fold aggregates")

    verify_original_scores(scores, pd.read_csv(PRIOR_STEP / PUBLISHED_SCORES, float_precision="round_trip"))
    if axis_run.digest(report / "bandwidths.csv") != axis_run.digest(PRIOR_STEP / PUBLISHED_BANDWIDTHS):
        raise ValueError("formal metric bandwidths differ from the accepted Step 06.1 bandwidths")
    return receipt, {"fold_scores": len(scores), "undefined_fold_scores": int(scores.value.isna().sum()),
                     "equal_fold_aggregates": len(summary), "exact_original_fold_scores": 27}


def audit_resources(audit, coordinator, destination):
    matrix_samples = read_samples(audit / "resource-samples.jsonl")
    coordinator_samples = read_samples(coordinator)
    scope_files = sorted(audit.glob("scope-*.json"))
    start = matrix_samples[0]["monotonic_seconds"]
    last = matrix_samples[-1]
    # The scope receipt is written after all fit phases. Convert its clock
    # reading to the sampler's monotonic timeline using the last joint sample.
    scope_finished = int(scope_files[-1].stem.removeprefix("scope-")) / 1e9
    end = last["monotonic_seconds"] + max(0, scope_finished - last["unix_seconds"])
    before = [row for row in coordinator_samples if row["monotonic_seconds"] <= start]
    between = [row for row in coordinator_samples if start < row["monotonic_seconds"] < end]
    after = [row for row in coordinator_samples if row["monotonic_seconds"] >= end]
    if not before or not after:
        raise ValueError("coordinator resource samples must bracket the completed matrix")
    coverage = [before[-1], *between, after[0]]
    gap = gaps_summary(coverage)
    if gap["gaps_over_12_seconds"] or not any(row["phase"] == "train" for row in matrix_samples):
        raise ValueError("resource samples must cover the training phases at ten-second intervals")
    hosts = [row["windows_host_total_and_available_kib"] for row in coverage]
    if any(not isinstance(host, str) for host in hosts):
        raise ValueError("coordinator samples are missing Windows host RAM")
    host_free = [int(host.split(",")[1]) for host in hosts]

    peak = max(axis_run.read_json(path)["memory_peak_bytes"] for path in scope_files)
    ram = min(row["wsl_available_kib"] for row in matrix_samples)
    gpu = min(row["gpu_total_mib"] - row["gpu_used_mib"] for row in matrix_samples)
    disk = min(row["disk_free_bytes"] for row in matrix_samples)
    if ram < 12 * GIB // 1024 or gpu < 512 or disk < GIB or peak > 10 * GIB:
        raise ValueError("formal matrix exceeded its resource caps or headroom")

    files = list(scope_files)
    for shape, fold in fit_coordinates(False):
        path = audit / "runs" / shape / f"fold-{fold}"
        for phase in ("register", *PHASES):
            measurements = sorted(path.glob(f"{phase}-*.resources"))
            if not measurements:
                raise ValueError(f"missing resource evidence for {shape}, fold {fold}, {phase}")
            for measurement in measurements:
                files.extend((measurement, measurement.with_suffix(".log")))
    with destination.open("w", encoding="utf-8") as stream:
        for sample in coverage:
            print(json.dumps(sample), file=stream)
    return {**gaps_summary(matrix_samples), "peak_scope_memory_bytes": peak,
            "minimum_wsl_available_kib": ram, "minimum_gpu_available_mib": gpu,
            "minimum_disk_free_bytes": disk,
            "phase_audit_sha256": axis_run.file_hashes(audit, [path.relative_to(audit).as_posix() for path in files]),
            "coordinator": {**gap, "minimum_windows_available_kib": min(host_free),
                            "minimum_wsl_available_kib": min(row["wsl"]["MemAvailable_kib"] for row in coverage),
                            "maximum_sampled_gpu_used_mib": max(
                                int(row["gpu_total_used_mib_and_utilization_percent"].split(",")[1])
                                for row in coverage)}}


def render_figures(candidate, scratch, activity):
    csvs = [f"figure/data/{name}" for name in FIGURE_TABLES]
    pngs = [f"figure/{name}" for name in FIGURES]
    renders = []
    for attempt in ("first", "second"):
        draw = scratch / attempt
        data = draw / "data"
        data.mkdir(parents=True)
        for source in csvs:
            shutil.copyfile(candidate / source, data / Path(source).name)
        command(RENDERER, ["--data", str(data), "--output", str(draw / "figure")], activity, f"render-{attempt}")
        renders.append(axis_run.file_hashes(draw, pngs))
    if renders[0] != renders[1]:
        raise ValueError("independent CSV-only Step 06.2 renders differ")
    for image in pngs:
        shutil.copyfile(scratch / "first" / image, candidate / image)
    return {
        "schema": "stage2-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": axis_run.platform.python_version(),
                        "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "independent_renders": 2,
        "csv_files_sha256": axis_run.file_hashes(candidate, csvs),
        "png_files_sha256": renders[0],
        "renderer_files_sha256": axis_run.file_hashes(STAGE, [
            RENDERER.relative_to(STAGE).as_posix(), (STEP / "execution/architecture_roster.py").relative_to(STAGE).as_posix(),
            "util/output/api/stage2_util/palette.py", "util/output/api/stage2_util/plotting.py",
        ]),
    }


def score(args, verify):
    paths, _ = read_matrix(args.matrix)
    originals, _ = load_originals(args.prior_matrix)
    score_or_verify(args.report, args.scratch, paths.values(), args.training, args.bandwidth_training,
                    False, verify, originals)
    print(json.dumps({"status": "PASS", "report_sha256": axis_run.digest(args.report / "manifest.json")}, sort_keys=True))


def publish(args):
    paths, fits = read_matrix(args.matrix)
    verify_formal_configurations(args, paths)
    _, original = load_originals(args.prior_matrix)
    fits.update({(fit["shape"], fit["fold"]): fit for fit in original["coordinates"]})

    args.activity.mkdir(parents=True, exist_ok=True)
    options = ["--matrix", str(args.matrix), "--prior-matrix", str(args.prior_matrix),
               "--training", str(args.training), "--bandwidth-training", str(args.bandwidth_training),
               "--report", str(args.report), "--scratch", str(args.scratch)]
    build = args.command == "build"
    with (axis_run.publication_lock(args.publication) as release,
          TemporaryDirectory(prefix="step06.2-formal-", dir=args.scratch) as temporary):
        if build and release.exists():
            raise ValueError("formal publication requires an unused output directory")
        if build and not args.report.exists():
            command(Path(__file__), ["score", *options], args.activity, "score")
        command(Path(__file__), ["verify-score", *options], args.activity, "verify-score")
        score_receipt, outcomes = verify_report(args.report, fits)

        candidate = Path(temporary) / "publication"
        sources = {destination: args.report / source for source, destination in PROJECTIONS.items()}
        sources[MATRIX_SAMPLES] = args.matrix_audit / "resource-samples.jsonl"
        for name, source in sources.items():
            target = candidate / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        coordinator = args.coordinator_samples if build else release / COORDINATOR_SAMPLES
        resources = audit_resources(args.matrix_audit, coordinator, candidate / COORDINATOR_SAMPLES)
        figures = render_figures(candidate, Path(temporary), args.activity)
        revision = args.code_revision if build else axis_run.read_json(release / RECEIPT)["matrix_code_revision"]
        receipt = {
            "status": "PASS", "ticket": 353, "matrix_code_revision": revision,
            "new_fits": 21, "reused_original_fits": 3, "independent_exact_new_replays": 21,
            "independent_exact_score_replay": True, **outcomes, "prior_release": original,
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "preflight_sha256": axis_run.digest(args.matrix / "preflight.json"),
            "report_manifest_sha256": axis_run.digest(args.report / "manifest.json"),
            "report_outputs_sha256": score_receipt["outputs"],
            "scientific_code_sha256": axis_run.file_hashes(STEP / "execution", [
                "architecture_model.py", "architecture_roster.py", "run_architecture.py", "verify_architecture.py",
                "matrix_architecture.py", "evaluate_architecture.py", "prior_architecture.py", Path(__file__).name,
            ]),
            "resources": resources, "outputs": axis_run.file_hashes(candidate, [*sources, COORDINATOR_SAMPLES]),
        }
        axis_run.write_json(candidate / RECEIPT, receipt)
        figures["scientific_receipt_sha256"] = axis_run.digest(candidate / RECEIPT)
        axis_run.write_json(candidate / FIGURE_RECEIPT, figures)
        if build:
            candidate.rename(release)
        else:
            for name in (RECEIPT, FIGURE_RECEIPT):
                if axis_run.digest(candidate / name) != axis_run.digest(release / name):
                    raise ValueError(f"Step 06.2 receipt differs from independently verified artifacts: {name}")
            axis_run.check_outputs(release, receipt)
            for key in ("csv_files_sha256", "png_files_sha256"):
                axis_run.check_outputs(release, {"outputs": figures[key]})

    return {"status": "PASS", "new_fits": 21, "reused_original_fits": 3, **outcomes,
            "receipt_sha256": axis_run.digest(args.publication / RECEIPT),
            "figure_receipt_sha256": axis_run.digest(args.publication / FIGURE_RECEIPT)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify", "score", "verify-score"))
    for name in ("matrix", "report", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--prior-matrix", type=Path, default=PRIOR_MATRIX)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    for name in ("matrix-audit", "publication", "activity", "coordinator-samples"):
        parser.add_argument(f"--{name}", type=Path)
    parser.add_argument("--code-revision")
    args = parser.parse_args()
    if args.command in ("score", "verify-score"):
        score(args, args.command == "verify-score")
        return
    if not all((args.matrix_audit, args.publication, args.activity)):
        parser.error("build and verify require --matrix-audit, --publication and --activity")
    if args.command == "build" and not all((args.coordinator_samples, args.code_revision)):
        parser.error("build requires --coordinator-samples and --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
