"""Build or verify the accepted Step 05 publication from the sealed formal matrix."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from audit_measured import audit_objectives, audit_resources
from audit_measured_scores import audit as audit_scores
from build_comparison import build as build_figures
from matrix_measured import check_registrations, preflight, verify_matrix
from run_module import DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING, STEP, axis_run
from training_data import TENV_DEFAULT

FORMAL_RECEIPT = "log/411-formal.json"
COORDINATOR_SAMPLES = "log/411-coordinator-resources.jsonl"
MATRIX_SAMPLES = "log/411-resources.jsonl"
ANALYSIS_TABLES = (
    "delta-metrics.csv",
    "fold-scores.csv",
    "primary-contrasts.csv",
    "primary-contrast-by-fold.csv",
    "primary-contrast-summary.csv",
)


def verify_scores(args):
    command = [
        sys.executable, str(STEP / "execution/evaluate_measured.py"), "verify",
        "--matrix", str(args.matrix),
        "--training", str(args.training),
        "--bandwidth-training", str(args.bandwidth_training),
        "--output", str(args.report),
        "--scratch", str(args.scratch),
    ]
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE, text=True)
    score_receipt = json.loads(result.stdout)
    if score_receipt != {"status": "PASS", "report_sha256": axis_run.digest(args.report / "manifest.json")}:
        raise ValueError("independent scoring did not reproduce the accepted report")
    return score_receipt


def build_candidate(args, candidate, coordinator, code_revision):
    paths, _ = verify_matrix(args.matrix, False)
    expected, plan = preflight(args.training, args.representations, args.source, False)
    preflight_path = args.matrix / "preflight.json"
    if axis_run.read_json(preflight_path) != plan:
        raise ValueError("the sealed matrix differs from the current registration plan")
    check_registrations(paths, expected, False)

    gate_path = STEP / "log/410-gate.json"
    gate = axis_run.read_json(gate_path)
    registration = gate["formal_registration"]
    preflight_sha = axis_run.digest(preflight_path)
    registered_files = sorted(path / "registration.json" for path in paths.values())
    registration_ledger = sha256("\n".join(map(axis_run.digest, registered_files)).encode()).hexdigest()
    if (gate["status"] != "PASS" or registration["new_fits"] != len(paths)
            or registration["preflight_sha256"] != preflight_sha
            or registration["registration_ledger_sha256"] != registration_ledger):
        raise ValueError("formal matrix differs from the preregistered 125-coordinate gate")

    scores = verify_scores(args)
    audit = audit_scores(args.report, False)
    report = axis_run.read_json(args.report / "manifest.json")
    if (report["new_runs"], report["reused_runs"], report["fold_scores"]) != (
        len(paths), 0, audit["fold_scores"]
    ):
        raise ValueError("the accepted report differs from the formal fit and score inventory")
    build_figures(args.report, candidate, False)
    projections = {name: f"figure/data/{name}" for name in ANALYSIS_TABLES}
    projections.update({"runs.csv": "log/runs.csv", "bandwidths.csv": "log/bandwidths.csv"})
    for name, destination in projections.items():
        shutil.copyfile(args.report / name, candidate / destination)

    shutil.copyfile(args.matrix_audit / "resource-samples.jsonl", candidate / MATRIX_SAMPLES)
    resources = audit_resources(args.matrix_audit, coordinator, candidate / COORDINATOR_SAMPLES, False)
    receipt = {
        "schema": "cpst-tenv-transigen-module-publication-v1",
        "status": "PASS",
        "ticket": 411,
        "matrix_code_revision": code_revision,
        "gate_receipt_sha256": axis_run.digest(gate_path),
        "matrix_manifest_sha256": axis_run.digest(args.matrix / "manifest.json"),
        "matrix_preflight_sha256": preflight_sha,
        "registration_ledger_sha256": registration_ledger,
        "training_seeds": report["training_seeds"],
        "held_out_folds": report["held_out_folds"],
        "registered_partitions": plan["folds"],
        "new_fits": len(paths),
        "independent_exact_replays": len(paths),
        "independent_exact_scores": audit["fold_scores"],
        "scoring_replay": scores,
        "report_outputs_sha256": report["outputs"],
        "objective_audit": audit_objectives(args.matrix, False),
        "score_table_audit": audit,
        "resources": resources,
        "outputs": axis_run.file_hashes(
            candidate,
            sorted(path.relative_to(candidate).as_posix() for path in candidate.rglob("*") if path.is_file()),
        ),
    }
    axis_run.write_json(candidate / FORMAL_RECEIPT, receipt)
    return receipt


def publish(args):
    with TemporaryDirectory(prefix="step05-publication-", dir=args.scratch) as temporary:
        candidate = Path(temporary) / "publication"
        if args.command == "build":
            with axis_run.publication_lock(args.publication) as target:
                if target.exists():
                    raise ValueError("publication requires an unused candidate directory")
                receipt = build_candidate(args, candidate, args.coordinator_samples, args.code_revision)
                candidate.rename(target)
        else:
            previous = axis_run.read_json(args.publication / FORMAL_RECEIPT)
            receipt = build_candidate(
                args, candidate, args.publication / COORDINATOR_SAMPLES, previous["matrix_code_revision"],
            )
            if receipt != previous:
                raise ValueError("independent verification differs from the published formal receipt")
            axis_run.check_outputs(args.publication, receipt)
    return {
        "status": "PASS",
        "new_fits": receipt["new_fits"],
        "independent_exact_scores": receipt["independent_exact_scores"],
        "formal_receipt_sha256": axis_run.digest(args.publication / FORMAL_RECEIPT),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    for name in ("matrix", "matrix-audit", "report", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    defaults = {
        "publication": STEP, "training": DEFAULT_TRAINING, "representations": DEFAULT_REPRESENTATIONS,
        "source": DEFAULT_SOURCE, "bandwidth-training": TENV_DEFAULT,
    }
    for name, default in defaults.items():
        parser.add_argument(f"--{name}", type=Path, default=default)
    parser.add_argument("--coordinator-samples", type=Path)
    parser.add_argument("--code-revision")
    args = parser.parse_args()
    if args.command == "build" and (args.coordinator_samples is None or not args.code_revision):
        parser.error("build requires --coordinator-samples and --code-revision")
    print(json.dumps(publish(args), sort_keys=True))


if __name__ == "__main__":
    main()
