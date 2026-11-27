"""Independently rescore and render a bounded or formal Stage 2 Step 06.3 study."""

import argparse
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import matplotlib
import numpy as np
import pandas as pd

from composition_model import STEP, axis_run
from composition_roster import coordinates, fit_coordinates
from evaluate_composition import METRICS, REPORT_SCHEMA
from matrix_composition import check_registrations, preflight, verify_matrix
from prior_composition import ROUND_ONE_MATRIX, ROUND_TWO_MATRIX
from prior_unnormalized import SOURCE_ROOT
from report_formal import command
from report_shape import gaps_summary, read_samples
from run_composition import AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING
from training_data import TENV_DEFAULT

STAGE = STEP.parents[1]
RENDERER = STEP / "execution/render_composition.py"
EVALUATOR = STEP / "execution/evaluate_composition.py"
FIGURE_TABLES = ("delta-metrics.csv", "fold-scores.csv", "epoch-objectives.csv",
                 "primary-contrasts.csv", "primary-contrast-summary.csv")
FIGURES = ("delta-metrics.png", "delta-metrics-paired-difference.png", "learning-curves.png",
           "primary-contrasts.png")
PROJECTIONS = {**{name: f"figure/data/{name}" for name in FIGURE_TABLES},
               "bandwidths.csv": "log/bandwidths.csv", "runs.csv": "log/runs.csv"}
FIGURE_RECEIPT = "log/figure-render.json"
CAP_GIB = 10
GIB = 1024 ** 3


def audit_resources(audit, coordinator, destination, smoke):
    matrix_samples = read_samples(audit / "resource-samples.jsonl")
    coordinator_samples = read_samples(coordinator)
    scope_files = sorted(audit.glob("scope-*.json"))
    if not scope_files or not any(row["phase"] == "train" for row in matrix_samples):
        raise ValueError("matrix audit requires scope completion and samples from training")
    start = matrix_samples[0]["monotonic_seconds"]
    last = matrix_samples[-1]
    # Convert the scope receipt's wall time onto the resource sampler's monotonic timeline.
    finished = int(scope_files[-1].stem.removeprefix("scope-")) / 1e9
    end = last["monotonic_seconds"] + max(0, finished - last["unix_seconds"])
    before = [row for row in coordinator_samples if row["monotonic_seconds"] <= start]
    between = [row for row in coordinator_samples if start < row["monotonic_seconds"] < end]
    after = [row for row in coordinator_samples if row["monotonic_seconds"] >= end]
    if not before or not after:
        raise ValueError("coordinator samples must bracket the completed matrix")
    coverage = [before[-1], *between, after[0]]
    gaps = gaps_summary(coverage)
    if gaps["gaps_over_12_seconds"]:
        raise ValueError("coordinator resource samples require ten-second coverage")
    host = [row["windows_host_total_and_available_kib"] for row in coverage]
    if any(not isinstance(value, str) for value in host):
        raise ValueError("resource receipt requires Windows host memory samples")
    host_free = min(int(value.split(",")[1]) for value in host)
    wsl_free = min(row["wsl"]["MemAvailable_kib"] for row in coverage)
    peak = max(axis_run.read_json(path)["memory_peak_bytes"] for path in scope_files)
    ram = min(row["wsl_available_kib"] for row in matrix_samples)
    gpu = min(row["gpu_total_mib"] - row["gpu_used_mib"] for row in matrix_samples)
    disk = min(row["disk_free_bytes"] for row in matrix_samples)
    if min(ram, wsl_free) < 12 * GIB // 1024 or gpu < 512 or disk < GIB or peak > CAP_GIB * GIB:
        raise ValueError("matrix exceeded its resource caps or available headroom")

    evidence = list(scope_files)
    for shape, fold, seed in fit_coordinates(smoke):
        path = audit / "runs" / shape / f"fold-{fold}" / f"seed-{seed}"
        for phase in ("register", "train", "predict", "replay"):
            resources = sorted(path.glob(f"{phase}-*.resources"))
            if not resources:
                raise ValueError(f"missing {phase} resource measurements for {(shape, fold, seed)}")
            for measurement in resources:
                evidence.extend((measurement, measurement.with_suffix(".log")))
    with destination.open("w", encoding="utf-8") as stream:
        for sample in coverage:
            print(json.dumps(sample), file=stream)
    return {
        **gaps_summary(matrix_samples),
        "scope_memory_max_bytes": CAP_GIB * GIB,
        "scope_swap_max_bytes": 8 * GIB,
        "peak_scope_memory_bytes": peak,
        "minimum_wsl_available_kib": ram,
        "minimum_gpu_available_mib": gpu,
        "minimum_disk_free_bytes": disk,
        "phase_audit_sha256": axis_run.file_hashes(audit, [path.relative_to(audit).as_posix() for path in evidence]),
        "coordinator": {**gaps, "minimum_windows_available_kib": host_free,
                        "minimum_wsl_available_kib": wsl_free,
                        "maximum_gpu_used_mib": max(int(row["gpu_total_used_mib_and_utilization_percent"]
                                                        .split(",")[1]) for row in coverage)},
    }


def render_figures(candidate, scratch, activity, smoke):
    csvs = [f"figure/data/{name}" for name in FIGURE_TABLES]
    pngs = [f"figure/{name}" for name in FIGURES]
    renders = []
    for attempt in ("first", "second"):
        draw = scratch / attempt
        for source in csvs:
            target = draw / source
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(candidate / source, target)
        command(RENDERER, ["--data", str(draw / "figure/data"), "--output", str(draw / "figure"),
                           *(["--smoke"] if smoke else [])], activity, f"render-{attempt}")
        renders.append(axis_run.file_hashes(draw, pngs))
    if renders[0] != renders[1]:
        raise ValueError("independent CSV-only Step 06.3 renders differ")
    for image in pngs:
        shutil.copyfile(scratch / "first" / image, candidate / image)
    return {
        "schema": "stage2-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": axis_run.platform.python_version(), "numpy": np.__version__,
                        "matplotlib": matplotlib.__version__, "pandas": pd.__version__},
        "independent_renders": 2,
        "csv_files_sha256": axis_run.file_hashes(candidate, csvs),
        "png_files_sha256": renders[0],
        "renderer_files_sha256": axis_run.file_hashes(STAGE, [
            RENDERER.relative_to(STAGE).as_posix(),
            *((STEP / "execution" / name).relative_to(STAGE).as_posix()
              for name in ("composition_roster.py", "composition_scores.py")),
            "util/output/api/stage2_util/palette.py", "util/output/api/stage2_util/plotting.py",
        ]),
    }


def verify_ledger(report, fits, smoke):
    receipt = axis_run.read_json(report / "manifest.json")
    roster = coordinates(smoke)
    if ((receipt["schema"], receipt["axis"], receipt["smoke"], receipt["runs"], receipt["new_runs"],
         receipt["metrics"], receipt["fold_scores"], receipt["status"])
            != (REPORT_SCHEMA, AXIS, smoke, len(roster), len(fits), list(METRICS), len(roster) * len(METRICS), "PASS")):
        raise ValueError("score receipt differs from the registered comparison")
    axis_run.check_outputs(report, receipt)
    runs = pd.read_csv(report / "runs.csv", float_precision="round_trip")
    summary = pd.read_csv(report / "delta-metrics.csv", float_precision="round_trip")
    contrasts = pd.read_csv(report / "primary-contrast-summary.csv", float_precision="round_trip")
    if list(runs[["shape", "fold", "seed"]].itertuples(index=False, name=None)) != list(roster):
        raise ValueError("scoring ledger differs from the registered fit roster")
    all_fits = dict(fits)
    if not smoke:
        for source in receipt["prior_release"].values():
            all_fits.update({(fit["shape"], fit["fold"], fit["seed"]): fit for fit in source["coordinates"]})
    for row in runs.itertuples(index=False):
        fit = all_fits[row.shape, row.fold, row.seed]
        if ((row.run_sha256, row.manifest_sha256, row.best_epoch, row.stopping_epoch)
                != (fit["registration_sha256"], fit["run_manifest_sha256"], fit["best_epoch"], fit["stopping_epoch"])):
            raise ValueError(f"scored fit identity differs at {(row.shape, row.fold, row.seed)}")
    return receipt, {
        "equal_coordinate_aggregates": len(summary),
        "primary_contrast_aggregates": len(contrasts),
        "fold_scores": receipt["fold_scores"],
        "exact_reused_scores": receipt["exact_reused_scores"],
        "primary_contrast_scores": receipt["primary_contrast_scores"],
    }


def publish(args):
    paths, fits = verify_matrix(args.matrix, args.smoke)
    expected, evidence = preflight(args.training, args.representations, args.source, args.smoke)
    if axis_run.read_json(args.matrix / "preflight.json") != evidence:
        raise ValueError("matrix preflight differs from the current frozen inputs or runtime")
    check_registrations(paths, expected, args.smoke)

    build = args.command == "build"
    receipt_name = "log/356-normalized-smoke.json" if args.smoke else "log/356-formal.json"
    matrix_samples = "log/356-normalized-smoke-resources.jsonl" if args.smoke else "log/356-resources.jsonl"
    coordinator_samples = ("log/356-normalized-smoke-coordinator-resources.jsonl" if args.smoke
                           else "log/356-coordinator-resources.jsonl")
    args.activity.mkdir(parents=True, exist_ok=True)
    evaluator_options = [
        "--matrix", str(args.matrix), "--output", str(args.report), "--scratch", str(args.scratch),
        "--training", str(args.training), "--bandwidth-training", str(args.bandwidth_training),
        "--round-one-matrix", str(args.round_one_matrix), "--round-two-matrix", str(args.round_two_matrix),
        "--prior-unnormalized", str(args.prior_unnormalized),
        *(["--smoke"] if args.smoke else []),
    ]
    with (axis_run.publication_lock(args.publication) as release,
          TemporaryDirectory(prefix="step06.3-report-", dir=args.scratch) as temporary):
        if build and release.exists():
            raise ValueError("publication build requires an unused output directory")
        if build and not args.report.exists():
            command(EVALUATOR, ["score", *evaluator_options], args.activity, "score")
        command(EVALUATOR, ["verify", *evaluator_options], args.activity, "verify-score")
        score_receipt, outcomes = verify_ledger(args.report, fits, args.smoke)
        candidate = Path(temporary) / "publication"
        sources = {destination: args.report / source for source, destination in PROJECTIONS.items()}
        sources[matrix_samples] = args.matrix_audit / "resource-samples.jsonl"
        for name, source in sources.items():
            target = candidate / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        coordinator = args.coordinator_samples if build else release / coordinator_samples
        resources = audit_resources(args.matrix_audit, coordinator, candidate / coordinator_samples, args.smoke)
        figures = render_figures(candidate, Path(temporary), args.activity, args.smoke)
        revision = args.code_revision if build else axis_run.read_json(release / receipt_name)["matrix_code_revision"]
        receipt = {
            "status": "PASS", "ticket": 356,
            "prediction_head_normalization": {shape: config["model"]["prediction_head_normalization"]
                                              for (shape, _, _), config in expected.items()
                                              if "prediction_head_normalization" in config["model"]},
            "matrix_code_revision": revision, "smoke": args.smoke,
            "new_fits": len(fits), "reused_fits": score_receipt["reused_runs"],
            "independent_exact_new_replays": len(fits), "independent_exact_score_replay": True,
            "training_seeds": score_receipt["training_seeds"], "held_out_folds": score_receipt["held_out_folds"],
            "registered_partition_sizes_and_hashes": evidence["folds"],
            "prior_release": score_receipt["prior_release"], **outcomes,
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "preflight_sha256": axis_run.digest(args.matrix / "preflight.json"),
            "report_manifest_sha256": axis_run.digest(args.report / "manifest.json"),
            "report_outputs_sha256": score_receipt["outputs"],
            "resources": resources,
            "scientific_code_sha256": axis_run.file_hashes(STEP / "execution", [
                "composition_model.py", "composition_roster.py", "composition_scores.py",
                "matrix_composition.py", "prior_composition.py", "prior_unnormalized.py", "run_composition.py",
                "verify_composition.py", "evaluate_composition.py", Path(__file__).name,
            ]),
            "outputs": axis_run.file_hashes(candidate, [*sources, coordinator_samples]),
        }
        axis_run.write_json(candidate / receipt_name, receipt)
        figures["scientific_receipt_sha256"] = axis_run.digest(candidate / receipt_name)
        axis_run.write_json(candidate / FIGURE_RECEIPT, figures)
        if build:
            candidate.rename(release)
        else:
            for name in (receipt_name, FIGURE_RECEIPT):
                if axis_run.digest(candidate / name) != axis_run.digest(release / name):
                    raise ValueError(f"Step 06.3 receipt differs from independently verified artifacts: {name}")
            axis_run.check_outputs(release, receipt)
            for name in ("csv_files_sha256", "png_files_sha256"):
                axis_run.check_outputs(release, {"outputs": figures[name]})
    return {"status": "PASS", "new_fits": len(fits), **outcomes,
            "receipt_sha256": axis_run.digest(args.publication / receipt_name),
            "figure_receipt_sha256": axis_run.digest(args.publication / FIGURE_RECEIPT)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    for name in ("matrix", "matrix-audit", "report", "publication", "activity", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--round-one-matrix", type=Path, default=ROUND_ONE_MATRIX)
    parser.add_argument("--round-two-matrix", type=Path, default=ROUND_TWO_MATRIX)
    parser.add_argument("--prior-unnormalized", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--coordinator-samples", type=Path)
    parser.add_argument("--code-revision")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.command == "build" and (args.coordinator_samples is None or not args.code_revision):
        parser.error("build requires --coordinator-samples and --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
