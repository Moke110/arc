"""Check Step 10 evidence from retained predictions after authorized weight disposal."""

import argparse
from functools import partial
import json
from pathlib import Path

from audit_dose_report import audit_report
from dose_formal import FIGURES_RECEIPT, RECEIPT, SCHEMA, checked_publication, verify
from dose_matrix import run_paths, verify_matrix
from dose_protocol import DEFAULT_STEP03, SCENARIO, coordinates, protocol
from dose_runs import axis_run, read_run
from cpst_util.checkpoint_disposition import WEIGHT, dispose_checkpoints, verify_prediction_archive

DISPOSITION = "checkpoint-disposition-v1.json"


def receipt_checker(matrix, runs):
    """Allow an absent weight only for a pinned replay receipt in the disposition."""

    def read(path, name, registration_sha, files):
        manifest = path / "manifest.json"
        pin = runs.get(manifest.relative_to(matrix).as_posix())
        if pin is None or WEIGHT not in files:
            return axis_run.verified_receipt(path, name, registration_sha, files)
        checkpoint = path / WEIGHT
        if manifest.is_symlink() or axis_run.digest(manifest) != pin or checkpoint.exists() or checkpoint.is_symlink():
            raise ValueError(f"Step 10 archived checkpoint state differs: {path}")

        def check_retained(root, receipt):
            retained = {key: digest for key, digest in receipt["outputs"].items() if key != WEIGHT}
            axis_run.check_outputs(root, {"outputs": retained})

        return axis_run.verified_receipt(path, name, registration_sha, files, output_check=check_retained)

    return read


def acceptance(bundle, publication, pin):
    """Bind the post-acceptance action to the approved crossed-dose release."""
    path = publication / RECEIPT
    if path.is_symlink() or axis_run.digest(path) != pin:
        raise ValueError("Step 10 acceptance receipt digest differs")
    record = axis_run.read_json(path)
    matrix = bundle / "matrix"
    if (record["schema"] != SCHEMA or record["status"] != "PASS" or record["ticket"] != 405
            or record["scenario"] != SCENARIO or record["seeds"] != list(protocol().seeds)
            or record["coordinates"] != len(coordinates())
            or record["independent_exact_replays"] != len(coordinates())
            or record["matrix_sha256"] != axis_run.digest(matrix / "manifest.json")
            or record["preflight_sha256"] != axis_run.digest(matrix / "preflight.json")
            or record["report_sha256"] != axis_run.digest(bundle / "report/manifest.json")
            or axis_run.read_json(publication / FIGURES_RECEIPT)["scientific_receipt_sha256"] != pin):
        raise ValueError("Step 10 archive is limited to the accepted installed crossed-dose release")
    manifest = axis_run.read_json(matrix / "manifest.json")
    paths = run_paths(matrix, False)
    fits = manifest["fits"]
    weights = {}
    if len(fits) != len(paths):
        raise ValueError("Step 10 archive requires the complete accepted roster")
    for fit, (key, run) in zip(fits, paths.items(), strict=True):
        relative = run.relative_to(matrix)
        if tuple(fit["coordinate"]) != key or fit["path"] != relative.as_posix():
            raise ValueError("Step 10 archive fit path differs from the accepted roster")
        weights[(relative / "manifest.json").as_posix()] = fit["manifest_sha256"]
    return weights


def verify_retained(bundle, publication, scratch, step03, runs):
    """Recompute the accepted report and figures from authenticated retained runs."""
    check = receipt_checker(bundle / "matrix", runs)

    def audit(_phase, bundle, report, data, step03, _resources):
        run_reader = partial(read_run, read_receipt=check)
        verified = verify_matrix(bundle / "matrix", False, run_reader=run_reader)
        return audit_report(report, data, verified, step03)

    return verify(bundle, publication, scratch, step03, audit=audit)


def dispose(bundle, publication, scratch, step03, acceptance_sha256):
    """Complete all accepted Step 10 dispositions in a single transaction."""
    runs = acceptance(bundle, publication, acceptance_sha256)
    pin = dispose_checkpoints(
        bundle / "matrix", receipt=DISPOSITION, acceptance=publication / RECEIPT,
        acceptance_sha256=acceptance_sha256, runs=runs,
        verify_full=lambda: checked_publication(bundle, publication),
        verify_retained=lambda _: verify_retained(bundle, publication, scratch, step03, runs),
    )
    return {"status": "PASS", "verification": "predictions-and-figures",
            "disposition_sha256": pin, "disposed_weights": len(runs)}


def verify_archive(bundle, publication, scratch, step03, acceptance_sha256, disposition_sha256):
    runs = acceptance(bundle, publication, acceptance_sha256)
    verify_prediction_archive(bundle / "matrix", receipt=DISPOSITION, receipt_sha256=disposition_sha256,
                              acceptance=publication / RECEIPT, acceptance_sha256=acceptance_sha256, runs=runs)
    result = verify_retained(bundle, publication, scratch, step03, runs)
    return {**result, "verification": "predictions-and-figures",
            "disposition_sha256": disposition_sha256, "disposed_weights": len(runs)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    for name, operation in (("dispose", dispose), ("verify", verify_archive)):
        command = commands.add_parser(name)
        for option in ("bundle", "publication", "scratch"):
            command.add_argument(f"--{option}", type=Path, required=True)
        command.add_argument("--step03", type=Path, default=DEFAULT_STEP03)
        command.add_argument("--acceptance-sha256", required=True)
        if operation is verify_archive:
            command.add_argument("--disposition-sha256", required=True)
        command.set_defaults(operation=operation)
    args = vars(parser.parse_args())
    args.pop("action")
    operation = args.pop("operation")
    paths = {name: value.resolve(strict=True) if isinstance(value, Path) else value
             for name, value in args.items()}
    print(json.dumps(operation(**paths), sort_keys=True))


if __name__ == "__main__":
    main()
