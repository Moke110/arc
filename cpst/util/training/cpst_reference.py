"""Verify CarPS-Transcriptome alpha against the nine accepted Step 06.4 checkpoints."""

import argparse
import csv
from itertools import product
import json
import os
from pathlib import Path
import subprocess
import sys
from functools import lru_cache
from hashlib import sha256

import torch

STAGE = Path(__file__).resolve().parents[2]
LEGACY_RUNTIME = STAGE / "util/training/legacy_transigen/modules/measured_targets"
ROUND_FOUR = STAGE.parent / "archived/spec27-transigen/round4"
sys.path.insert(0, str(LEGACY_RUNTIME / "execution"))

from run_measured import axis_run
from .cpst_model import CPSTAdapter, KPGT_WIDTH, MEASURED_SHARED, PROJECTION_WIDTH, objective
from cpst_util.release_names import original_axis_manifest_sha256
from measured_roster import SEEDS, fit_coordinates
from run_measured import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                          bind)
from verify_measured import read_run
from util.training.provenance import git_code_hashes
from util.training.transigen import predict_test, record_predictions

ALPHA_REVISION = "a7613b74fc23f9f48d8ea0d4479542c008d701db"
ALPHA_ROOT = ("cpst", "07_cpst_alpha")
HISTORICAL_ROUND_FOUR = ("cpst", "06_tsg_shape", "06.4_round4")
RELEASE = (*ALPHA_ROOT, "output/release.json")
ACCEPTANCE = (*ALPHA_ROOT, "log/362-alpha.json")
RELEASE_ID = "cpst-alpha"
RELEASE_SCHEMA = "cpst-alpha-release-v1"
ALPHA_RUNTIME_REVISION = "825f48c96a0968c11b0bffa9a8b06686b0d16e03"
REFERENCE_RECEIPT = "log/359-formal.json"
REFERENCE_ROSTER = tuple(product(range(3), SEEDS))
REPORT_FILES = {"log/runs.csv": "runs.csv", "log/bandwidths.csv": "bandwidths.csv",
                "figure/data/fold-scores.csv": "fold-scores.csv"}
FORMAL_MANIFESTS = {"matrix/manifest.json": "matrix_manifest_sha256",
                    "matrix/preflight.json": "matrix_preflight_sha256",
                    "report/manifest.json": "report_manifest_sha256"}


@lru_cache(maxsize=None)
def historical_bytes(path):
    """Read an immutable artifact from the revision that accepted it."""
    return subprocess.check_output(
        ["git", "show", f"{ALPHA_REVISION}:{'/'.join(path)}"], cwd=axis_run.REPO,
    )


def historical_json(path):
    return json.loads(historical_bytes(path))


def historical_digest(path):
    return sha256(historical_bytes(path)).hexdigest()


def alpha_release_digest():
    return historical_digest(RELEASE)


def defaults():
    """Return the runtime values frozen in the release."""
    return {"axis": AXIS, "architecture": MEASURED_SHARED, "chemical_representation": "KPGT",
            "chemical_width": KPGT_WIDTH, "projection_width": PROJECTION_WIDTH,
            "objective": objective(), "validation_seed": axis_run.VALIDATION_SEED,
            "test_seed": axis_run.TEST_SEED, "folds": list(range(3)), "training_seeds": list(SEEDS)}


def verify_release():
    """Verify the fixed defaults and historical Git code behind the alpha release."""
    release = historical_json(RELEASE)
    if (release["schema"], release["release"], release["default"]) != (RELEASE_SCHEMA, RELEASE_ID, defaults()):
        raise ValueError("alpha defaults differ from the fixed CarPS-Transcriptome alpha release")
    pinned = release["runtime_sha256"]
    actual = git_code_hashes(axis_run.REPO, ALPHA_RUNTIME_REVISION, tuple(sorted(pinned)))
    if actual != pinned:
        changed = sorted(path for path in actual.keys() | pinned.keys() if actual.get(path) != pinned.get(path))
        raise ValueError(f"alpha runtime dependencies differ: {changed}")
    reference = release["step_06_4"]
    original = reference["code_sha256"]
    if git_code_hashes(axis_run.REPO, reference["matrix_code_revision"], tuple(sorted(original))) != original:
        raise ValueError("Step 06.4 source differs from its original registered Git revision")
    return release


def verify_reference_release(round_four, release):
    """Read the accepted Step 06.4 ledgers independently of fit payload availability."""
    round_four = Path(round_four)
    reference = release["step_06_4"]
    receipt_path = round_four / REFERENCE_RECEIPT
    if axis_run.digest(receipt_path) != reference["receipt_sha256"]:
        raise ValueError("Step 06.4 publication receipt identity differs")
    receipt = axis_run.read_json(receipt_path)
    if (receipt["status"], receipt["ticket"], receipt["matrix_code_revision"]) != (
            "PASS", 359, reference["matrix_code_revision"]):
        raise ValueError("Step 06.4 publication status or revision differs")
    formal = round_four / "output/formal"
    for path, field in FORMAL_MANIFESTS.items():
        if axis_run.digest(formal / path) != receipt[field]:
            raise ValueError(f"Step 06.4 sealed formal manifest differs: {path}")
    report = axis_run.read_json(formal / "report/manifest.json")
    if report["outputs"] != receipt["report_outputs_sha256"]:
        raise ValueError("Step 06.4 formal report differs from the accepted receipt")
    for published, recorded in REPORT_FILES.items():
        expected = receipt["outputs"][published]
        if (historical_digest((*HISTORICAL_ROUND_FOUR, published)) != expected
                or axis_run.digest(formal / "report" / recorded) != expected):
            raise ValueError(f"Step 06.4 published report differs: {published}")

    counts = axis_run.read_json(formal / "matrix/preflight.json")
    matrix = axis_run.read_json(formal / "matrix/manifest.json")
    full_roster = list(fit_coordinates(False))
    if (matrix["preflight_sha256"] != receipt["matrix_preflight_sha256"]
            or [tuple(key) for key in counts["new_fits"]] != full_roster
            or [(row["shape"], row["fold"], row["seed"]) for row in matrix["fits"]] != full_roster):
        raise ValueError("Step 06.4 matrix roster differs from the accepted plan")
    matrix_rows = {(row["fold"], row["seed"]): row for row in matrix["fits"] if row["shape"] == MEASURED_SHARED}
    path = round_four / "log/runs.csv"
    with path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["shape"] == MEASURED_SHARED]
    if [(int(row["fold"]), int(row["seed"])) for row in rows] != list(REFERENCE_ROSTER):
        raise ValueError("Step 06.4 selected run ledger does not cover the selected arm")
    return formal / "matrix", matrix_rows, dict(zip(REFERENCE_ROSTER, rows, strict=True))


def verify_reference(ledger, release, fold, seed):
    """Verify and return one selected checkpoint's path, registration and receipts."""
    matrix, registrations, published = ledger
    registered = registrations[fold, seed]
    row = published[fold, seed]
    path = matrix / "runs" / MEASURED_SHARED / f"fold-{fold}" / f"seed-{seed}"
    if path.relative_to(matrix).as_posix() != registered["path"]:
        raise ValueError(f"Step 06.4 selected coordinate has an invalid path: {(fold, seed)}")
    key, (_, config, replay, training, _) = read_run(path, False)
    expected = {"run_sha256": axis_run.digest(path / "registration.json"),
                "manifest_sha256": axis_run.digest(path / "manifest.json"),
                "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
                "trace_sha256": replay["outputs"]["learning_curve.csv"],
                "prediction_sha256": replay["outputs"]["prediction.npy"]}
    if (key != (MEASURED_SHARED, fold, seed)
            or config["code_sha256"] != release["step_06_4"]["code_sha256"]
            or expected["run_sha256"] != registered["registration_sha256"]
            or expected["manifest_sha256"] != registered["run_manifest_sha256"]
            or any(row[field] != digest for field, digest in expected.items())
            or row["source_shape"] != MEASURED_SHARED or row["source_ticket"] != "359"):
        raise ValueError(f"Step 06.4 selected coordinate differs from accepted ledgers: {(fold, seed)}")
    return path, config, replay, training


def verify_protocol(config, data, vectors, source):
    """Rebind relocatable inputs and compare the complete non-path protocol."""
    options = {**config["arguments"], "source": str(source)}
    observed = bind(options, data, vectors)
    if config["axis_family_manifest_sha256"] != observed["axis_family_manifest_sha256"]:
        observed["axis_family_manifest_sha256"] = original_axis_manifest_sha256(observed["axis_family_manifest_sha256"])
    ignored = {"arguments", "nonce", "code_sha256", "experiment"}
    if {key: value for key, value in observed.items() if key not in ignored} != {
            key: value for key, value in config.items() if key not in ignored}:
        raise ValueError("current alpha inputs or measured protocol differ from Step 06.4")


def replay_alpha(round_four, training, representations, source):
    """Reproduce validation and held-out values through the alpha entry point."""
    release = verify_release()
    ledger = verify_reference_release(round_four, release)
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    axis_run.kpgt.configure_backend("cuda")
    coordinates = []
    with axis_run.open_inputs(training, representations) as (data, vectors):
        for fold, seed in REFERENCE_ROSTER:
            path, config, replay, trained = verify_reference(ledger, release, fold, seed)
            verify_protocol(config, data, vectors, source)
            adapter = CPSTAdapter(source, config["model"]["transform"])
            model = adapter.build(seed, "cuda", n_genes=len(data.genes))
            registration_sha = replay["registration_sha256"]
            values, test = predict_test(path, adapter, model, data, config, registration_sha, trained)
            result = record_predictions(path, None, values, test, data, config, replay=True)
            coordinates.append({"fold": fold, "seed": seed, "registration_sha256": registration_sha,
                                "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
                                "validation_prediction_sha256": replay["outputs"]["validation_prediction.npy"],
                                "prediction_sha256": replay["outputs"]["prediction.npy"],
                                "validation_conditions": len(config["partitions"]["valid"]),
                                "test_replay": "exact", **result})
            del model
    return {"status": "PASS", "release": release["release"], "release_sha256": alpha_release_digest(),
            "validation_seed": axis_run.VALIDATION_SEED, "test_seed": axis_run.TEST_SEED,
            "step_06_4_receipt_sha256": release["step_06_4"]["receipt_sha256"], "coordinates": coordinates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round-four", type=Path, default=ROUND_FOUR)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    result = replay_alpha(**vars(args))
    if historical_json(ACCEPTANCE) != result:
        raise ValueError("alpha replay differs from the accepted evidence")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
