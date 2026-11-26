"""Publish and independently verify the complete formal Step 08 comparison."""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
from platform import python_version
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

import pandas as pd

from audit_chem_report import audit
from audit_chem_resources import audit_resources
from chem_figures import FIGURES, RENDERER, replay
from chem_matrix import preflight, run_paths, verify_registrations
from chem_reference import GATE, reference_plan
from chem_report import (METRICS, PRIMARY_METRICS, SCHEMA as REPORT_SCHEMA, report_tables, verified_runs,
                         verify_initializations, write_report)
from chem_runs import REFERENCE, ROUND_FOUR, STAGE, STEP, add_inputs, axis_run, coordinates
from training_data import TENV_DEFAULT

RECEIPT = "log/364-formal.json"
FIGURE_RECEIPT = "log/figure-render.json"
MATRIX_SAMPLES = "log/364-resources.jsonl"
COORDINATOR_SAMPLES = "log/364-coordinator-resources.jsonl"
VARIANT = "execution/selected-representation.json"
PROJECTIONS = {
    **{f"{name}.csv": f"figure/data/{name}.csv" for name in
       ("fold-scores", "paired-by-fold", "paired-summary", "epoch-objectives")},
    "runs.csv": "log/runs.csv", "bandwidths.csv": "log/bandwidths.csv",
}
REPORT_TABLES = {*PROJECTIONS, "condition-scores.csv"}
RNG_FIELDS = ("downstream_sha256", "torch_cpu_rng_sha256", "torch_cuda_rng_sha256")


def select_configuration(summary):
    """Keep KPGT unless one arm leads both metrics beyond paired repeat variability."""
    primary = summary.loc[summary.metric.isin(PRIMARY_METRICS)]
    winners, comparisons = {}, {}
    for metric, direction in zip(PRIMARY_METRICS, (1, -1), strict=True):
        rows = primary.loc[primary.metric.eq(metric)].set_index("arm")
        improvement = direction * rows.score_difference
        leading = improvement.index[improvement.eq(improvement.max())].tolist()
        winners[metric] = leading[0] if len(leading) == 1 else None
        comparisons[metric] = rows
    first, second = (winners[name] for name in PRIMARY_METRICS)
    default = REFERENCE
    reason = "primary metrics have different winners or a tied winner"
    if first == second and first is not None:
        if first == REFERENCE:
            reason = "KPGT leads both primary metrics"
        else:
            pcc, rmse = [comparisons[name].loc[first] for name in PRIMARY_METRICS]
            if (pcc.n_paired == rmse.n_paired == 9
                    and pcc.score_difference > pcc.score_difference_sd
                    and -rmse.score_difference > rmse.score_difference_sd):
                default, reason = first, "one representation leads both metrics beyond paired sample SD"
            else:
                reason = "paired primary improvements are small relative to repeat variability"
    return {"working_default": default, "metric_winners": winners, "reason": reason}


def score(args):
    matrix = args.bundle / "matrix"
    actual_plan = axis_run.read_json(args.bundle / "reuse-preflight.json")
    reused, reproduced_plan = reference_plan(args.round_four, args.training, args.representations,
                                             args.source, args.bandwidth_training, accepted_plan=actual_plan)
    if actual_plan != reproduced_plan:
        raise ValueError("fresh KPGT reference replay differs from the frozen reuse decisions")
    fresh = verified_runs(matrix, False)
    if set(fresh) != set(coordinates()) - reused.keys():
        raise ValueError("fresh matrix and KPGT decisions do not cover exactly 63 coordinates")
    paths = run_paths(matrix, fresh)
    expected, plan = preflight(args.training, args.representations, args.source, False)
    if axis_run.read_json(matrix / "preflight.json") != {**plan, "fresh_fits": [list(key) for key in fresh]}:
        raise ValueError("formal matrix preflight differs from current inputs, runtime or partitions")
    verify_registrations(paths, expected, False)
    rng_by_fold_seed = verify_initializations(fresh)
    for ref in actual_plan["decisions"]:
        arm, fold, seed = ref["coordinate"]
        if (arm, fold, seed) in reused and rng_by_fold_seed[fold, seed] != tuple(
            ref["initialization"][key] for key in RNG_FIELDS
        ):
            raise ValueError(f"fresh fits are not paired to the retained KPGT initialization: {(fold, seed)}")
    report = report_tables(fresh, reused, args.training, args.representations, args.source,
                           args.bandwidth_training, False)
    return write_report(args.output, report, False)


def checked_report(report):
    manifest = axis_run.read_json(report / "manifest.json")
    if (manifest["schema"], manifest["smoke"], manifest["scores"]) != (REPORT_SCHEMA, False, len(coordinates()) * len(METRICS)):
        raise ValueError("invalid formal Step 08 report")
    if (set(manifest["csv_sha256"]) != REPORT_TABLES
            or manifest["csv_sha256"] != axis_run.file_hashes(report, REPORT_TABLES)):
        raise ValueError("formal report tables differ from their sealed content")
    return manifest


def score_in_process(args, work, name):
    report = work / name
    resources = (work / name).with_suffix(".resources")
    command = ["/usr/bin/time", "-v", "-o", str(resources), sys.executable, __file__, "score"]
    for field in ("bundle", "round_four", "training", "representations", "source", "bandwidth_training"):
        command.extend((f"--{field.replace('_', '-')}", str(getattr(args, field))))
    log_path = (work / name).with_suffix(".log")
    try:
        with log_path.open("x", encoding="utf-8") as log:
            subprocess.run([*command, "--output", str(report)], check=True, stdout=log, stderr=log)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"{name} failed:\n{log_path.read_text()[-4000:]}") from error
    rss = next(int(line.rsplit(": ", 1)[1]) for line in resources.read_text().splitlines()
               if "Maximum resident set size (kbytes):" in line)
    return report, rss


def figure_receipt(figure, scientific_receipt_sha256):
    render_deps = [RENDERER, *(STAGE / f"util/output/api/stage2_util/{name}.py" for name in ("palette", "plotting"))]
    return {
        "schema": "stage2-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": python_version(), **{name: version(name) for name in ("matplotlib", "numpy", "pandas")}},
        "independent_renders": figure["independent_csv_only_renders"],
        "csv_files_sha256": {f"figure/data/{name}": digest for name, digest in figure["csv_sha256"].items()},
        "png_files_sha256": {f"figure/{name}": digest for name, digest in figure["png_sha256"].items()},
        "renderer_files_sha256": axis_run.file_hashes(STAGE, [path.relative_to(STAGE).as_posix()
                                                            for path in render_deps]),
        "scientific_receipt_sha256": scientific_receipt_sha256,
    }


def publish(args):
    build = args.mode == "build"
    matrix, activity, report = (args.bundle / name for name in ("matrix", "matrix-audit", "report"))
    previous = {} if build else axis_run.read_json(args.publication / RECEIPT)
    with TemporaryDirectory(prefix="step08-formal-", dir=args.scratch) as temporary:
        work = Path(temporary)
        independently_scored, peak = score_in_process(args, work, "score-replay")
        replayed = checked_report(independently_scored)
        if build:
            first, first_peak = score_in_process(args, work, "score-first")
            if replayed != checked_report(first):
                raise ValueError("independent formal scoring processes disagree")
            score_peaks = [first_peak, peak]
            if not report.exists():
                first.rename(report)
        else:
            score_peaks = previous["score_process_peak_rss_kib"]
        if replayed != checked_report(report):
            raise ValueError("independently regenerated scores differ from the formal report")

        numerical = audit(report)
        registrations = axis_run.read_json(matrix / "preflight.json")
        fresh = tuple(map(tuple, registrations["fresh_fits"]))
        reference_path = args.bundle / "reuse-preflight.json"
        reference = axis_run.read_json(reference_path)
        coordinator = args.coordinator if build else args.publication / COORDINATOR_SAMPLES
        coverage, resources = audit_resources(activity, coordinator, fresh, bracket=build)
        configuration = select_configuration(pd.read_csv(report / "paired-summary.csv", float_precision="round_trip"))

        candidate = work / "publication" if build else args.publication
        figure_path = candidate / "figure"
        if build:
            for source, relative in PROJECTIONS.items():
                target = candidate / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(report / source, target)
            shutil.copyfile(activity / "resource-samples.jsonl", candidate / MATRIX_SAMPLES)
            (candidate / COORDINATOR_SAMPLES).write_text("".join(json.dumps(row) + "\n" for row in coverage))
        picture = replay(figure_path / "data", figure_path, work, verify=not build)
        names = [*PROJECTIONS.values(), MATRIX_SAMPLES, COORDINATOR_SAMPLES, *(f"figure/{name}" for name in FIGURES)]
        if configuration["working_default"] != REFERENCE:
            names.append(VARIANT)
            if build:
                (candidate / VARIANT).parent.mkdir()
                arm = configuration["working_default"]
                axis_run.write_json(candidate / VARIANT, {"representation": arm,
                                    "alpha_release_sha256": reference["alpha_release_sha256"],
                                    "native_width": registrations["frozen_representations"][arm]["width"]})
        receipt = {
            "schema": "carps-vae-alpha-chemical-formal-v1", "status": "PASS", "ticket": 364,
            "total_coordinates": len(coordinates()), "held_out_folds": [0, 1, 2], "training_seeds": [17, 42, 2026],
            "sample_seeds": {"validation": axis_run.VALIDATION_SEED, "test": axis_run.TEST_SEED},
            "representations": registrations["frozen_representations"],
            "parent_count": registrations["parent_count"], "condition_count": registrations["condition_count"],
            "new_exact_independent_replays": len(fresh), "exact_scoring_processes": 2,
            "reused_kpgt_coordinates": len(coordinates()) - len(fresh), "reference_decisions": reference["decisions"],
            "alpha_release_sha256": reference["alpha_release_sha256"], "alpha_receipt_sha256": reference["alpha_receipt_sha256"],
            "bounded_gate_sha256": axis_run.digest(GATE), "reuse_plan_sha256": axis_run.digest(reference_path),
            "matrix_sha256": axis_run.digest(matrix / "manifest.json"),
            "preflight_sha256": axis_run.digest(matrix / "preflight.json"),
            "report_sha256": axis_run.digest(report / "manifest.json"), "configuration": configuration,
            "scores_audit": numerical, "resources": resources, "score_process_peak_rss_kib": score_peaks,
            "release_files_sha256": axis_run.file_hashes(candidate, names),
        }
        if build:
            axis_run.write_json(candidate / RECEIPT, receipt)
        elif previous != receipt:
            raise ValueError("formal acceptance evidence differs from the sealed publication")
        scientific_sha = axis_run.digest(candidate / RECEIPT)
        plots = figure_receipt(picture, scientific_sha)
        if build:
            axis_run.write_json(candidate / FIGURE_RECEIPT, plots)
            candidate.rename(args.publication)
        elif plots != axis_run.read_json(candidate / FIGURE_RECEIPT):
            raise ValueError("formal figure receipt differs from its independent CSV renders")
    return {"status": "PASS", "receipt_sha256": scientific_sha, "scores": numerical["score_rows"],
            "new_fits": len(fresh), "reused_fits": len(coordinates()) - len(fresh), **configuration}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("score", "build", "verify"))
    add_inputs(parser)
    parser.add_argument("--round-four", type=Path, default=ROUND_FOUR)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--bundle", type=Path, default=STEP / "output/formal")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--publication", type=Path, default=STEP)
    parser.add_argument("--coordinator", type=Path)
    parser.add_argument("--scratch", type=Path)
    args = parser.parse_args()
    if args.mode == "score":
        if args.output is None:
            parser.error("score requires --output")
        result = score(args)
    else:
        if args.scratch is None:
            parser.error("publication requires --scratch")
        if args.mode == "build" and (args.coordinator is None or args.publication.exists()):
            parser.error("build requires --coordinator and an unused --publication directory")
        result = publish(args)
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
