"""Build, independently verify and install the accepted Step 10 comparison."""

import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
from platform import python_version
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from dose_figures import RENDERER, replay
from dose_gate import verify_resource_samples
from dose_interruption import InterruptedResources
from dose_matrix import BOUNDED_RECEIPT, verify_bounded_gate, verify_matrix
from dose_protocol import DEFAULT_STEP03, HELD_DOSES, SCENARIO, STAGE, STEP, coordinates, protocol
from dose_report import PLOT_TABLES, PRIMARY_METRICS, REPORT_METRICS, SCHEMA as REPORT_SCHEMA
from dose_runs import axis_run
from audit_chem_resources import audit_resources
from cpst_util.comparison import FIGURES

SCHEMA = "cpst-dose-formal-v2"
RECEIPT = "log/405-formal.json"
FIGURES_RECEIPT = "log/figure-render.json"
COORDINATOR = "log/405-coordinator-resources.jsonl"
MATRIX_SAMPLES = "log/405-matrix-resources.jsonl"
AUDIT_PEAKS = ("score", "audit")
RESOURCES = tuple(f"log/405-resources/{phase}.resources" for phase in AUDIT_PEAKS)
PUBLISHED_TABLES = ("fold-scores.csv", "foldwise-summary.csv", "summary.csv", "runs.csv", "bandwidths.csv")
REPORT_TABLES = (*PUBLISHED_TABLES, "condition-scores.csv", "epoch-objectives.csv")
PLOT_CSVS = tuple(f"figure/data/{filename}" for filename in PLOT_TABLES)
PLOT_PNGS = tuple(f"figure/{filename}" for filename in FIGURES)
PUBLICATION_FILES = (*PLOT_CSVS, COORDINATOR, MATRIX_SAMPLES, *RESOURCES,
                     *(f"log/{filename}" for filename in PUBLISHED_TABLES))


def checked_report(report):
    """Check the formal delivery inventory against the scored report receipt."""
    record = axis_run.read_json(report / "manifest.json")
    if (record["schema"] != REPORT_SCHEMA or record["bounded"] is not False
            or record["primary_metrics"] != list(PRIMARY_METRICS)
            or record["scores"] != len(coordinates()) * len(REPORT_METRICS)
            or {path.name for path in report.iterdir()} != {*REPORT_TABLES, "manifest.json"}
            or record["outputs"] != axis_run.file_hashes(report, REPORT_TABLES)):
        raise ValueError("formal report differs from its sealed inventory")
    return record


def measured_report_action(phase, bundle, report, data, step03, resources):
    """Run report derivation or numerical audit in its own measured process."""
    execution = STEP / "execution"
    script = execution / ("dose_report.py" if phase == "score" else "audit_dose_report.py")
    report_flag = "--output" if phase == "score" else "--report"
    command = ["/usr/bin/time", "-v", "-o", str(resources), sys.executable, str(script),
               "--matrix", str(bundle / "matrix"), report_flag, str(report),
               "--figure-data", str(data), "--step03", str(step03)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"Step 10 {phase} failed:\n{result.stderr}")
    return json.loads(result.stdout)


def process_peak(resources):
    lines = resources.read_text().splitlines()
    if "Exit status: 0" not in {line.strip() for line in lines}:
        raise ValueError(f"unsuccessful Step 10 report measurement: {resources}")
    return next(int(line.rsplit(": ", 1)[1]) for line in lines
                if "Maximum resident set size (kbytes):" in line)


def figure_receipt(picture, scientific_hash):
    dependencies = (RENDERER, STEP / "execution/dose_protocol.py",
                    *(STAGE / f"util/output/api/cpst_util/{name}.py" for name in ("comparison", "palette", "plotting")))
    return {
        "schema": "cpst-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": python_version(), **{name: version(name) for name in ("matplotlib", "numpy", "pandas")}},
        "independent_renders": picture["independent_csv_only_renders"],
        "csv_files_sha256": {f"figure/data/{name}": digest for name, digest in picture["csv_sha256"].items()},
        "png_files_sha256": {f"figure/{name}": digest for name, digest in picture["png_sha256"].items()},
        "renderer_files_sha256": axis_run.file_hashes(STAGE, [path.relative_to(STAGE).as_posix()
                                                             for path in dependencies]),
        "scientific_receipt_sha256": scientific_hash,
    }


def scientific_receipt(bundle, publication, report, numerical, *, recovery_sha256=None):
    matrix = bundle / "matrix"
    preflight = axis_run.read_json(matrix / "preflight.json")
    verify_bounded_gate(preflight)
    scored = checked_report(report)
    if (numerical["status"] != "PASS" or numerical["coordinates"] != len(coordinates())
            or numerical["score_rows"] != scored["scores"]):
        raise ValueError("formal publication requires an independent complete numerical audit")
    if axis_run.file_hashes(publication / "log", PUBLISHED_TABLES) != {
            name: scored["outputs"][name] for name in PUBLISHED_TABLES}:
        raise ValueError("public score tables differ from the independently audited report")
    activity = bundle / "matrix-audit"
    coordinator = publication / COORDINATOR
    interrupted = InterruptedResources(bundle, recovery_sha256) if recovery_sha256 else None
    _, resources = audit_resources(activity, coordinator, coordinates(), interrupted=interrupted)
    final_scope = sorted(activity.glob("scope-*.json"))[-1]
    final_measurement = int(final_scope.stem.removeprefix("scope-")) / 1e9
    cadence = (verify_resource_samples(coordinator, activity / "resource-samples.jsonl", final_measurement)
               if interrupted is None else resources["coordinator"])
    if axis_run.digest(publication / MATRIX_SAMPLES) != resources["matrix_samples_sha256"]:
        raise ValueError("public CUDA samples differ from the formal observations")
    return {
        "schema": SCHEMA, "status": "PASS", "ticket": 405, "coordinates": len(coordinates()),
        "independent_exact_replays": len(coordinates()), "scenario": SCENARIO,
        "held_dose_um": list(HELD_DOSES), "seeds": list(protocol().seeds),
        "primary_metrics": list(PRIMARY_METRICS), "folds": preflight["folds"],
        "bounded_receipt_sha256": axis_run.digest(STEP / BOUNDED_RECEIPT),
        "matrix_sha256": axis_run.digest(matrix / "manifest.json"),
        "preflight_sha256": axis_run.digest(matrix / "preflight.json"),
        "report_sha256": axis_run.digest(report / "manifest.json"),
        "numerical_audit": numerical, "resources": {**resources, "sampling": cadence},
        "report_peak_rss_kib": {phase: process_peak(publication / path)
                                for phase, path in zip(AUDIT_PEAKS, RESOURCES, strict=True)},
        "publication_sha256": axis_run.file_hashes(publication, PUBLICATION_FILES),
    }


def build(bundle, publication, scratch, step03, *, recovery_sha256=None):
    """Stage only verified public artifacts and retain the full local report."""
    report = bundle / "report"
    resolved_publication, resolved_report = publication.resolve(), report.resolve()
    if resolved_publication.is_relative_to(resolved_report) or resolved_report.is_relative_to(resolved_publication):
        raise ValueError("publication destination overlaps the local report directory")
    with axis_run.publication_lock(publication) as output:
        if output.exists():
            raise ValueError("build requires an unused publication directory")
        with TemporaryDirectory(prefix="step10-formal-", dir=scratch) as temporary:
            work = Path(temporary)
            candidate = work / "publication"
            figure = candidate / "figure"
            data = figure / "data"
            scored = work / "report"
            (candidate / RESOURCES[0]).parent.mkdir(parents=True)
            measured_report_action("score", bundle, scored, data, step03, candidate / RESOURCES[0])
            numerical = measured_report_action("audit", bundle, scored, data, step03, candidate / RESOURCES[1])
            for name in PUBLISHED_TABLES:
                shutil.copyfile(scored / name, candidate / "log" / name)
            for name in RESOURCES:
                measurement = candidate / name
                lines = measurement.read_text().splitlines(keepends=True)
                measurement.write_text("".join(
                    line for line in lines if not line.lstrip().startswith("Command being timed:")
                ))
            activity = bundle / "matrix-audit"
            interrupted = InterruptedResources(bundle, recovery_sha256) if recovery_sha256 else None
            coverage, _ = audit_resources(activity, bundle / "coordinator.jsonl", coordinates(),
                                          bracket=True, interrupted=interrupted)
            (candidate / COORDINATOR).write_text("".join(json.dumps(row) + "\n" for row in coverage))
            shutil.copyfile(activity / "resource-samples.jsonl", candidate / MATRIX_SAMPLES)
            picture = replay(data, figure, work)
            receipt = scientific_receipt(bundle, candidate, scored, numerical, recovery_sha256=recovery_sha256)
            axis_run.write_json(candidate / RECEIPT, receipt)
            scientific_hash = axis_run.digest(candidate / RECEIPT)
            axis_run.write_json(candidate / FIGURES_RECEIPT, figure_receipt(picture, scientific_hash))
            with axis_run.publication_lock(report) as sealed:
                if not sealed.exists():
                    scored.rename(sealed)
                elif checked_report(sealed) != checked_report(scored):
                    raise ValueError("existing scored report differs from the independently audited report")
            candidate.rename(output)
    return {"status": "PASS", "receipt_sha256": scientific_hash, **numerical}


def verify(bundle, publication, scratch, step03):
    """Recompute the numerical audit and regenerate all three PNGs from CSVs."""
    report = bundle / "report"
    figure = publication / "figure"
    with TemporaryDirectory(prefix="step10-verify-", dir=scratch) as temporary:
        work = Path(temporary)
        numerical = measured_report_action("audit", bundle, report, figure / "data", step03, work / "audit.resources")
        published = axis_run.read_json(publication / RECEIPT)
        current = scientific_receipt(bundle, publication, report, numerical,
                                     recovery_sha256=published_recovery(published))
        if current != published:
            raise ValueError("scientific receipt differs from the verified Step 10 evidence")
        picture = replay(figure / "data", figure, work, verify=True)
    scientific_hash = axis_run.digest(publication / RECEIPT)
    if figure_receipt(picture, scientific_hash) != axis_run.read_json(publication / FIGURES_RECEIPT):
        raise ValueError("figure receipt differs from the independent CSV-only reproductions")
    return {"status": "PASS", "receipt_sha256": scientific_hash, **numerical}


def published_recovery(receipt):
    return receipt["resources"]["coordinator"].get("interruption_receipt_sha256")


def checked_publication(bundle, publication):
    """Recheck passed receipts and their immutable fit, report and figure bytes."""
    receipt = axis_run.read_json(publication / RECEIPT)
    figures = axis_run.read_json(publication / FIGURES_RECEIPT)
    scientific_hash = axis_run.digest(publication / RECEIPT)
    if (receipt["schema"] != SCHEMA or receipt["status"] != "PASS"
            or receipt["independent_exact_replays"] != len(coordinates())
            or (figures["schema"], figures["status"], figures["independent_renders"])
            != ("cpst-figure-render-v1", "PASS", 2)
            or figures["scientific_receipt_sha256"] != scientific_hash
            or receipt["publication_sha256"] != axis_run.file_hashes(publication, PUBLICATION_FILES)
            or figures["png_files_sha256"] != axis_run.file_hashes(publication, PLOT_PNGS)
            or figures["csv_files_sha256"] != axis_run.file_hashes(publication, PLOT_CSVS)
            or figures["renderer_files_sha256"] != axis_run.file_hashes(STAGE, figures["renderer_files_sha256"])):
        raise ValueError("formal publication bytes differ from their passing receipts")
    verify_matrix(bundle / "matrix", False)
    if receipt != scientific_receipt(bundle, publication, bundle / "report", receipt["numerical_audit"],
                                     recovery_sha256=published_recovery(receipt)):
        raise ValueError("formal fit bundle differs from the published receipt")
    return scientific_hash


def install(bundle, publication, destination):
    """Install an immutable bundle atomically without duplicating checkpoint bytes."""
    with axis_run.publication_lock(destination) as output:
        if output.exists():
            scientific_hash = checked_publication(output, publication)
        else:
            scientific_hash = checked_publication(bundle, publication)
            with TemporaryDirectory(prefix=".step10-install-", dir=output.parent) as temporary:
                candidate = Path(temporary) / "formal"
                candidate.mkdir()
                shutil.copytree(bundle / "matrix", candidate / "matrix", copy_function=os.link)
                for name in ("matrix-audit", "report"):
                    shutil.copytree(bundle / name, candidate / name)
                candidate.rename(output)
    return {"status": "PASS", "output": str(output), "receipt_sha256": scientific_hash}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    for action in (build, verify, install):
        command = commands.add_parser(action.__name__)
        command.add_argument("--bundle", type=Path, default=STEP / "output/formal-crossed-v1")
        command.add_argument("--publication", type=Path, default=STEP)
        if action == build:
            command.add_argument("--recovery-sha256")
        if action == install:
            command.add_argument("--destination", type=Path, default=STEP / "output/formal-crossed-v1")
        else:
            command.add_argument("--scratch", type=Path, required=True)
            command.add_argument("--step03", type=Path, default=DEFAULT_STEP03)
        command.set_defaults(function=action)
    args = parser.parse_args()
    paths = {name: value.resolve() for name, value in vars(args).items() if isinstance(value, Path)}
    if args.action == "build":
        paths["recovery_sha256"] = args.recovery_sha256
    print(json.dumps(args.function(**paths), sort_keys=True))


if __name__ == "__main__":
    main()
