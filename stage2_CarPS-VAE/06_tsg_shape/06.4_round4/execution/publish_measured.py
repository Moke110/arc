"""Publish and independently verify the formal Step 06.4 report."""

import argparse
import csv
from importlib.metadata import version
import json
from pathlib import Path
from platform import python_version
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from audit_measured_figures import CSV_NAMES, replay as replay_figures
from audit_measured_scores import audit as audit_tables
from measured_roster import coordinates

STEP = Path(__file__).resolve().parents[1]
EXECUTION = STEP / "execution"
FORMAL_RECEIPT = "log/359-formal.json"
FIGURE_RECEIPT = "log/figure-render.json"
MATRIX_SAMPLES = "log/359-resources.jsonl"
COORDINATOR_SAMPLES = "log/359-coordinator-resources.jsonl"
PROJECTIONS = {**{name: f"figure/data/{name}" for name in CSV_NAMES},
               "runs.csv": "log/runs.csv", "bandwidths.csv": "log/bandwidths.csv"}


def coordinate(record):
    return record["shape"], int(record["fold"]), int(record["seed"])


def verify_ledger(table, fresh, retained):
    reused = {coordinate(fit): fit for fit in retained}
    roster = coordinates(False)
    if len(reused) != len(retained) or set(fresh) & set(reused) or set(fresh) | set(reused) != set(roster):
        raise ValueError("fit manifests differ from the registered fresh and retained rosters")
    with table.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if tuple(map(coordinate, rows)) != roster:
        raise ValueError("scoring ledger differs from the registered fit roster")
    all_fits = fresh | reused
    for row in rows:
        key = coordinate(row)
        fit = all_fits[key]
        expected = (fit["registration_sha256"], fit["run_manifest_sha256"],
                    fit["best_epoch"], fit["stopping_epoch"])
        if (row["run_sha256"], row["manifest_sha256"], int(row["best_epoch"]),
                int(row["stopping_epoch"])) != expected:
            raise ValueError(f"scored fit identity differs at {key}")


def score(phase, args, name):
    evidence = args.activity / f"{name}.resources"
    record = args.activity / f"{name}.json"
    with record.open("x", encoding="utf-8") as output:
        subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(evidence), sys.executable, str(EXECUTION / "evaluate_measured.py"),
             phase, "--matrix", str(args.matrix), "--training", str(args.training), "--bandwidth-training",
             str(args.bandwidth_training), "--source-release", str(args.source_release),
             "--output", str(args.report), "--scratch", str(args.scratch)],
            check=True, stdout=output,
        )
    receipt = json.loads(record.read_text(encoding="utf-8"))
    for line in evidence.read_text(encoding="utf-8").splitlines():
        if "Maximum resident set size (kbytes):" in line:
            return receipt, int(line.split(":")[-1].strip())
    raise ValueError(f"missing peak RSS for {phase}")


def figure_receipt(candidate, scratch, scientific_receipt_sha256):
    receipt_path = scratch / "render.json"
    replay_figures(candidate / "figure/data", candidate / "figure", receipt_path, scratch, False, False)
    replay = json.loads(receipt_path.read_text(encoding="utf-8"))
    return {
        "schema": "stage2-figure-render-v1", "status": "PASS", "step": STEP.name,
        "environment": {"python": python_version(), **{name: version(name) for name in ("matplotlib", "numpy", "pandas")}},
        "independent_renders": replay["independent_csv_only_renders"],
        "csv_files_sha256": {f"figure/data/{name}": digest for name, digest in replay["csv_sha256"].items()},
        "png_files_sha256": {f"figure/{name}": digest for name, digest in replay["png_sha256"].items()},
        "renderer_files_sha256": replay["code_sha256"],
        "scientific_receipt_sha256": scientific_receipt_sha256,
    }


def publish(args):
    from audit_measured import audit_objectives, audit_resources
    from matrix_measured import check_registrations, preflight, verify_matrix
    from measured_model import axis_run

    build = args.command == "build"
    gate_path = STEP / "log/358-gate.json"
    gate = axis_run.read_json(gate_path)
    if ((gate["status"], gate["ticket"]) != ("PASS", 358)
            or tuple(map(coordinate, gate["fits"])) != coordinates(True)):
        raise ValueError("publication requires the accepted bounded CUDA gate")
    paths, fits = verify_matrix(args.matrix, False)
    expected, plan = preflight(args.training, args.representations, args.source, False)
    if axis_run.read_json(args.matrix / "preflight.json") != plan:
        raise ValueError("matrix preflight differs from the current frozen inputs or runtime")
    check_registrations(paths, expected, False)
    args.activity.mkdir()
    with (axis_run.publication_lock(args.publication) as release,
          TemporaryDirectory(prefix="step064-publication-", dir=args.scratch) as temporary):
        scratch = Path(temporary)
        previous = {} if build else axis_run.read_json(release / FORMAL_RECEIPT)
        if build and release.exists():
            raise ValueError("publication build requires an unused output directory")
        if build:
            scored, score_peak = score("verify" if args.report.exists() else "score", args, "score")
        reproduced, replay_peak = score("verify", args, "verify")
        report_hash = axis_run.digest(args.report / "manifest.json")
        for result in (scored, reproduced) if build else (reproduced,):
            if result != {"status": "PASS", "report_sha256": report_hash}:
                raise ValueError("independent score process differs from the sealed report")
        report = axis_run.read_json(args.report / "manifest.json")
        axis_run.check_outputs(args.report, report)
        prior = report["prior_release"]
        verify_ledger(args.report / "runs.csv", fits, prior["coordinates"])
        audit = audit_tables(args.report, False)

        candidate = scratch / "publication"
        sources = {destination: args.report / source for source, destination in PROJECTIONS.items()}
        sources[MATRIX_SAMPLES] = args.matrix_audit / "resource-samples.jsonl"
        for name, source in sources.items():
            target = candidate / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        coordinator = args.coordinator_samples if build else release / COORDINATOR_SAMPLES
        resources = audit_resources(args.matrix_audit, coordinator, candidate / COORDINATOR_SAMPLES, False)
        projection = axis_run.file_hashes(candidate, [*sources, COORDINATOR_SAMPLES])
        receipt = {
            "schema": "stage2-tenv-transigen-measured-encoder-publication-v1", "status": "PASS", "ticket": 359,
            "bounded_gate_sha256": axis_run.digest(gate_path),
            "matrix_code_revision": args.code_revision if build else previous["matrix_code_revision"],
            "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
            "matrix_preflight_sha256": axis_run.digest(args.matrix / "preflight.json"),
            "report_manifest_sha256": report_hash, "report_outputs_sha256": report["outputs"],
            "training_seeds": report["training_seeds"], "held_out_folds": report["held_out_folds"],
            "registered_partitions": plan["folds"], "new_fits": len(fits), "reused_fits": len(prior["coordinates"]),
            "independent_exact_new_replays": len(fits), "independent_exact_scores": audit["fold_scores"],
            "exact_reused_scores": report["exact_reused_scores"], "prior_release": prior,
            "objectives": audit_objectives(args.matrix, False), "paired_table_audit": audit,
            "resources": resources, "score_process_peak_rss_kib": (
                {"score": score_peak, "verify": replay_peak} if build else previous["score_process_peak_rss_kib"]),
            "outputs": projection,
        }
        axis_run.write_json(candidate / FORMAL_RECEIPT, receipt)
        figure_record = figure_receipt(candidate, scratch, axis_run.digest(candidate / FORMAL_RECEIPT))
        axis_run.write_json(candidate / FIGURE_RECEIPT, figure_record)
        if build:
            candidate.rename(release)
        else:
            for name in (FORMAL_RECEIPT, FIGURE_RECEIPT):
                if axis_run.digest(candidate / name) != axis_run.digest(release / name):
                    raise ValueError(f"published evidence differs from independent verification: {name}")
            for outputs in (projection, figure_record["png_files_sha256"]):
                axis_run.check_outputs(release, {"outputs": outputs})
    return {"status": "PASS", "new_fits": len(fits), "reused_fits": len(prior["coordinates"]),
            "score_rows": audit["fold_scores"], "receipt_sha256": axis_run.digest(release / FORMAL_RECEIPT),
            "figure_receipt_sha256": axis_run.digest(release / FIGURE_RECEIPT)}


def main():
    from prior_measured import ROUND_THREE
    from run_measured import DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING
    from training_data import TENV_DEFAULT

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    for name in ("matrix", "matrix-audit", "report", "activity", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    for name, default in (("publication", STEP), ("training", DEFAULT_TRAINING),
                          ("representations", DEFAULT_REPRESENTATIONS), ("source", DEFAULT_SOURCE),
                          ("bandwidth-training", TENV_DEFAULT), ("source-release", ROUND_THREE)):
        parser.add_argument(f"--{name}", type=Path, default=default)
    parser.add_argument("--coordinator-samples", type=Path)
    parser.add_argument("--code-revision")
    args = parser.parse_args()
    if args.command == "build" and (args.coordinator_samples is None or not args.code_revision):
        parser.error("build requires --coordinator-samples and --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
