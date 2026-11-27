"""Fit, audit and verify the formal CarPS-Transcriptome prediction-axis comparison."""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
from platform import python_version
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from axis_figures import RENDERER, replay
from axis_matrix import execute, verify_bounded_gate, verify_matrix
from axis_protocol import REFERENCE_AXIS, SEEDS, STAGE, STEP, coordinates
from axis_reference import reference_plan, verify_publication
from axis_report import (PLOT_TABLES, REPORT_TABLES, SCHEMA as REPORT_SCHEMA, add_report_inputs,
                         figure_data, read_reused, report_tables, write_report)
from axis_scores import REPORT_METRICS
from axis_runs import add_inputs, axis_run
from audit_axis_report import audit_plotted_csvs, audit_report
from matrix import enforce_record
from cpst_util.derivation.release import replace
from cpst_util.publication import replace_release

SCHEMA = "cpst-unimol-axis-formal-v2"
HISTORICAL_REPORT_SHA256 = "7ed1f80d1673f7d99a129527e8c5ab07fa5a2696852c0cba3cb9461375f5f5f0"
RECEIPT = "log/371-formal.json"
FIGURE_RECEIPT = "log/figure-render.json"
MATRIX_SAMPLES = "log/368-resources.jsonl"
COORDINATOR_SAMPLES = "log/368-coordinator-resources.jsonl"
REFERENCE_PLAN = "log/368-reference.json"
PROJECTIONS = {name: f"log/{name}" for name in ("bandwidths.csv", "runs.csv", "summary.csv")}


def fresh_coordinates(plan):
    """Resolve a complete per-coordinate reuse plan to the remaining ordered fits."""
    roster = coordinates()
    references = {key for key in roster if key[0] == REFERENCE_AXIS}
    decisions = {tuple(row["coordinate"]): row for row in plan["decisions"]}
    if set(decisions) != references or len(decisions) != len(plan["decisions"]):
        raise ValueError("reuse audit must decide each of the nine Step 08 coordinates once")
    reused = set()
    for key, row in decisions.items():
        if (row["decision"] == "reuse" and row.get("validation_replay") == "exact"
                and row.get("test_replay") == "exact"):
            reused.add(key)
        elif row["decision"] != "refit" or not row.get("reason"):
            raise ValueError(f"reuse decision lacks exact replay or a refit reason: {key}")
    return tuple(key for key in roster if key not in reused)


def fit(output, step03, step08, representations, source):
    """Audit sealed references, persist the disposition, and independently replay each new fit."""
    _, gate_hash = verify_bounded_gate()
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "reuse-preflight.json"
    plan = reference_plan(step08, step03, representations, source)
    if not plan_path.exists() or "bounded_gate_sha256" in axis_run.read_json(plan_path):
        plan["bounded_gate_sha256"] = gate_hash
    fresh = fresh_coordinates(plan)
    enforce_record(plan_path, plan)
    print(json.dumps({"phase": "reuse-audit", "reused": len(coordinates()) - len(fresh), "fresh": len(fresh)}),
          flush=True)
    return execute(output / "matrix", output / "matrix-audit", step03, representations, source, fresh=fresh)


def collect_runs(bundle, step08):
    """Verify the full fresh and reused roster and pair the two sources at initialization."""
    plan = axis_run.read_json(bundle / "reuse-preflight.json")
    expected_fresh = fresh_coordinates(plan)
    if (len(expected_fresh) != len(coordinates())
            and plan["step08_acceptance_sha256"] != verify_publication(step08)[-1]):
        raise ValueError("reuse plan is not bound to the sealed Step 08 publication")
    fresh = verify_matrix(bundle / "matrix", False)
    reused = read_reused(step08, plan)
    if (tuple(fresh) != expected_fresh or set(fresh) & set(reused)
            or set(fresh) | set(reused) != set(coordinates())):
        raise ValueError("formal run evidence does not cover each registered coordinate once")
    for decision in plan["decisions"]:
        if decision["decision"] == "reuse":
            _, fold, seed = decision["coordinate"]
            common = {key: val for key, val in decision["paired_initialization"].items()
                      if key != "axis_projection_shapes"}
            for key, (_, _, _, training, _) in fresh.items():
                if key[1:] == (fold, seed):
                    fresh_common = {field: val for field, val in training["initialization"].items()
                                    if field != "axis_projection_shapes"}
                    if fresh_common != common:
                        raise ValueError(f"fresh initialization and training RNG are not reference-paired: {key}")
    return fresh, reused


def checked_report(report):
    """Verify that the sealed formal report contains the expected files and metrics."""
    record = axis_run.read_json(report / "manifest.json")
    score_count = len(coordinates()) * len(REPORT_METRICS)
    if (record["schema"] != REPORT_SCHEMA or record["bounded"] is not False or record["scores"] != score_count
            or {path.name for path in report.iterdir()} != {*REPORT_TABLES, "manifest.json"}
            or record["outputs"] != axis_run.file_hashes(report, REPORT_TABLES)):
        raise ValueError("formal report manifest differs from the registered full-axis report")
    return record


def run_report_action(args, action, work, report, data=None):
    """Run scoring and audit as separate measured processes."""
    resources = work / f"{action}.resources"
    command = ["/usr/bin/time", "-v", "-o", str(resources), sys.executable, __file__, action,
               "--report", str(report)]
    for field in ("bundle", "step03", "step04", "step08"):
        command.extend((f"--{field.replace('_', '-')}", str(getattr(args, field).resolve())))
    if data is not None:
        command.extend(("--figure-data", str(data)))
    log = work / f"{action}.log"
    try:
        with log.open("x", encoding="utf-8") as output:
            subprocess.run(command, check=True, stdout=output, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"{action} failed:\n{log.read_text()[-4000:]}") from error
    peak = next(int(line.rsplit(": ", 1)[1]) for line in resources.read_text().splitlines()
                if "Maximum resident set size (kbytes):" in line)
    return json.loads(log.read_text().splitlines()[-1]), peak


def figure_receipt(picture, scientific_hash):
    dependencies = [RENDERER, STEP / "execution/axis_protocol.py",
                    *(STAGE / f"util/output/api/cpst_util/{name}.py" for name in ("comparison", "palette", "plotting"))]
    return {
        "schema": "cpst-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": python_version(), **{name: version(name) for name in ("matplotlib", "numpy", "pandas")}},
        "independent_renders": picture["independent_csv_only_renders"],
        "csv_files_sha256": {f"figure/data/{path}": digest for path, digest in picture["csv_sha256"].items()},
        "png_files_sha256": {f"figure/{path}": digest for path, digest in picture["png_sha256"].items()},
        "renderer_files_sha256": axis_run.file_hashes(STAGE, [path.relative_to(STAGE).as_posix()
                                                             for path in dependencies]),
        "scientific_receipt_sha256": scientific_hash,
    }


def scientific_evidence(bundle, publication, reference, numerical, resources, score_peaks):
    """Bind complete run lineage and audit evidence to the public comparison tables."""
    expected = [
        *PROJECTIONS.values(), MATRIX_SAMPLES, COORDINATOR_SAMPLES, REFERENCE_PLAN,
        *(f"figure/data/{name}" for name in PLOT_TABLES),
    ]
    fresh = fresh_coordinates(reference)
    preflight = axis_run.read_json(bundle / "matrix/preflight.json")
    _, gate_hash = verify_bounded_gate(formal_preflight=preflight)
    if reference.get("bounded_gate_sha256", gate_hash) != gate_hash:
        raise ValueError("formal execution used different bounded gate evidence")
    if any(axis_run.digest(publication / relative) != axis_run.digest(bundle / "report" / name)
           for name, relative in PROJECTIONS.items()):
        raise ValueError("public comparison tables differ from the independently verified report")
    if axis_run.digest(publication / MATRIX_SAMPLES) != resources["matrix_samples_sha256"]:
        raise ValueError("published resource samples differ from measured CUDA activity")
    reference_hash = axis_run.digest(bundle / "reuse-preflight.json")
    if axis_run.digest(publication / REFERENCE_PLAN) != reference_hash:
        raise ValueError("published reference decisions differ from the audited plan")
    return {
        "schema": SCHEMA, "status": "PASS", "ticket": 371, "coordinates": len(coordinates()),
        "cold_drug_folds": list(range(3)), "training_seeds": list(SEEDS),
        "new_exact_replays": len(fresh), "step08_exact_reuses": len(coordinates()) - len(fresh),
        "step08_acceptance_sha256": reference["step08_acceptance_sha256"],
        "bounded_gate_sha256": gate_hash,
        "reuse_plan_sha256": reference_hash,
        "matrix_sha256": axis_run.digest(bundle / "matrix/manifest.json"),
        "preflight_sha256": axis_run.digest(bundle / "matrix/preflight.json"),
        "report_sha256": axis_run.digest(bundle / "report/manifest.json"),
        "training_manifest_sha256": preflight["training_manifests_sha256"], "unimol": preflight["unimol"],
        "fold_partitions": preflight["folds"],
        "numerical_audit": numerical, "resources": resources, "score_process_peak_rss_kib": score_peaks,
        "publication_sha256": axis_run.file_hashes(publication, expected),
    }


def finish(args):
    """Stage the report and figures, or independently verify a staged or installed publication."""
    from audit_chem_resources import audit_resources

    build = args.action == "build"
    bundle = args.bundle.resolve()
    report = bundle / "report"
    prior = {} if build else axis_run.read_json(args.publication / RECEIPT)
    with TemporaryDirectory(prefix="step09-formal-", dir=args.scratch) as temporary:
        work = Path(temporary)
        candidate = work / "publication" if build else args.publication
        figure = candidate / "figure"
        scoring_report = report
        if build:
            print(json.dumps({"phase": "formal-scoring"}), flush=True)
            scoring_report = work / "scored-report"
            _, score_peak = run_report_action(args, "score", work, scoring_report, figure / "data")
            if report.exists() and checked_report(report) != checked_report(scoring_report):
                raise ValueError("existing candidate report differs from fresh scoring")
        print(json.dumps({"phase": "independent-full-axis-audit"}), flush=True)
        numerical, audit_peak = run_report_action(args, "audit", work, scoring_report, figure / "data")
        score_peaks = {"scoring": score_peak, "independent_audit": audit_peak} if build else prior["score_process_peak_rss_kib"]
        reference = axis_run.read_json(bundle / "reuse-preflight.json")
        activity = bundle / "matrix-audit"
        coordinator = args.coordinator if build else candidate / COORDINATOR_SAMPLES
        coverage, resources = audit_resources(activity, coordinator, fresh_coordinates(reference), bracket=build)
        if build:
            if not report.exists():
                scoring_report.rename(report)
            (candidate / "log").mkdir(parents=True)
            for name, relative in PROJECTIONS.items():
                shutil.copyfile(report / name, candidate / relative)
            shutil.copyfile(activity / "resource-samples.jsonl", candidate / MATRIX_SAMPLES)
            shutil.copyfile(bundle / "reuse-preflight.json", candidate / REFERENCE_PLAN)
            (candidate / COORDINATOR_SAMPLES).write_text("".join(json.dumps(row) + "\n" for row in coverage))
        print(json.dumps({"phase": "independent-csv-figure-replays"}), flush=True)
        picture = replay(figure / "data", figure, work, verify=not build)
        receipt = scientific_evidence(bundle, candidate, reference, numerical, resources, score_peaks)
        if build:
            axis_run.write_json(candidate / RECEIPT, receipt)
        elif prior != receipt:
            raise ValueError("the installed scientific receipt differs from its verified evidence")
        scientific_hash = axis_run.digest(candidate / RECEIPT)
        figures = figure_receipt(picture, scientific_hash)
        if build:
            axis_run.write_json(candidate / FIGURE_RECEIPT, figures)
            with axis_run.publication_lock(args.publication) as publication:
                if publication.exists():
                    raise ValueError("build requires an unused report publication directory")
                candidate.rename(publication)
        elif figures != axis_run.read_json(candidate / FIGURE_RECEIPT):
            raise ValueError("figure receipt differs from the independent CSV-only renders")
    return {"status": "PASS", "receipt_sha256": scientific_hash, "new_fits": receipt["new_exact_replays"],
            "reused_fits": receipt["step08_exact_reuses"], **numerical}


def install(bundle, publication, destination, step08):
    """Install the verified full-axis report while preserving its registered fit bundle."""
    receipt = axis_run.read_json(publication / RECEIPT)
    plots = axis_run.read_json(publication / FIGURE_RECEIPT)
    if ((receipt.get("schema"), receipt.get("status")) != (SCHEMA, "PASS")
            or plots["scientific_receipt_sha256"] != axis_run.digest(publication / RECEIPT)):
        raise ValueError("formal publication is missing a passed, figure-bound receipt")
    files = receipt["publication_sha256"] | plots["png_files_sha256"]
    if axis_run.file_hashes(publication, files) != files:
        raise ValueError("published comparison files differ from their passed receipt")
    fit_fields = {"reuse-preflight.json": "reuse_plan_sha256", "matrix/manifest.json": "matrix_sha256",
                  "matrix/preflight.json": "preflight_sha256"}
    fit_files = {name: receipt[field] for name, field in fit_fields.items()}

    with axis_run.publication_lock(destination) as output:
        candidate = output if output.exists() else bundle
        report = bundle / "report"
        if (fit_files != axis_run.file_hashes(candidate, fit_files)
                or checked_report(report)["scores"] != receipt["numerical_audit"]["score_rows"]
                or axis_run.digest(report / "manifest.json") != receipt["report_sha256"]):
            raise ValueError("formal fit bundle or full-axis report differs from the passed receipt")
        preflight = axis_run.read_json(candidate / "matrix/preflight.json")
        _, gate_hash = verify_bounded_gate(formal_preflight=preflight)
        if receipt["bounded_gate_sha256"] != gate_hash:
            raise ValueError("passed receipt does not bind the verified Step 09 bounded gate")
        fresh, reused = collect_runs(candidate, step08)
        audit_plotted_csvs(publication / "figure/data", report, {**fresh, **reused})
        if not output.exists():
            candidate.rename(output)
        else:
            with axis_run.publication_lock(output / "report") as current:
                current_hash = axis_run.digest(current / "manifest.json")
                if current_hash == receipt["report_sha256"]:
                    checked_report(current)
                elif current_hash == HISTORICAL_REPORT_SHA256:
                    replace_historical_report(current, report)
                else:
                    raise ValueError("existing report is neither the verified publication nor its sealed predecessor")
    return {"status": "PASS", "output": str(output), "receipt_sha256": axis_run.digest(publication / RECEIPT)}


def replace_historical_report(current, report):
    def verify_predecessor(directory):
        predecessor = axis_run.read_json(directory / "manifest.json")["outputs"]
        if axis_run.file_hashes(directory, predecessor) != predecessor:
            raise ValueError("sealed predecessor report bytes differ from its manifest")

    verify_predecessor(current)
    with TemporaryDirectory(prefix=".report-install-", dir=current.parent) as temporary:
        candidate = shutil.copytree(report, Path(temporary) / "replacement")
        replace_release(candidate, current, replace, verify_existing=verify_predecessor)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    run = commands.add_parser("fit")
    add_inputs(run)
    run.add_argument("--output", type=Path, required=True)
    for action in ("score", "audit", "build", "verify", "install"):
        command = commands.add_parser(action)
        command.add_argument("--bundle", type=Path, default=STEP / "output/formal")
        if action != "install":
            command.add_argument("--step03", type=Path, default=STAGE / "03_data_train")
            add_report_inputs(command)
        if action in ("score", "audit"):
            command.add_argument("--report", type=Path, required=True)
        else:
            command.add_argument("--publication", type=Path, default=STEP)
            if action != "install":
                command.add_argument("--scratch", type=Path, required=True)
    for command in (run, *(commands.choices[name] for name in ("score", "audit", "build", "verify", "install"))):
        command.add_argument("--step08", type=Path, default=STAGE / "08_cpst_chem-emb")
    commands.choices["score"].add_argument("--figure-data", type=Path, required=True)
    commands.choices["audit"].add_argument("--figure-data", type=Path)
    commands.choices["build"].add_argument("--coordinator", type=Path, required=True)
    commands.choices["install"].add_argument("--destination", type=Path, default=STEP / "output/formal")
    args = parser.parse_args()
    if args.action == "fit":
        result = fit(args.output.resolve(), args.step03.resolve(), args.step08.resolve(),
                     args.representations.resolve(), args.source.resolve())
    elif args.action in ("score", "audit"):
        fresh, reused = collect_runs(args.bundle, args.step08)
        if args.action == "score":
            tables = report_tables(fresh, reused, args.step03, args.step04, False)
            result = write_report(args.report, tables, False)
            figure_data(tables, args.figure_data)
        else:
            checked_report(args.report)
            runs = {**fresh, **reused}
            result = audit_report(args.report, runs, args.step03, args.step04, args.step08)
            if args.figure_data is not None:
                audit_plotted_csvs(args.figure_data, args.report, runs)
    elif args.action == "install":
        result = install(args.bundle, args.publication, args.destination, args.step08)
    else:
        if args.action == "build" and args.publication.exists():
            parser.error("build requires an unused --publication directory")
        result = finish(args)
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
