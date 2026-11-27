"""Audit the sealed Step 08 UniMol coordinates before their reuse in Step 09."""

import argparse
import json
from pathlib import Path

from axis_protocol import ARM, REFERENCE_AXIS, coordinates
from axis_model import AxisAdapter, paired_initialization
from axis_runs import add_inputs, axis_run, configure_cuda
import chem_runs
from tenv_axes import tenv_axis_output

DEFECTS = (OSError, ValueError, KeyError, RuntimeError)


def verify_publication(step08):
    """Validate the official acceptance and the hashed run ledger and matrix."""
    step08 = Path(step08)
    ledger = step08 / "log/runs.csv"
    acceptance = step08 / "log/364-formal.json"
    accepted = axis_run.read_json(acceptance)
    matrix = step08 / "output/formal/matrix"
    manifest = matrix / "manifest.json"
    preflight = matrix / "preflight.json"
    if (accepted["status"] != "PASS" or accepted["ticket"] != 364
            or accepted["alpha_release_sha256"] != axis_run.digest(chem_runs.ALPHA_MANIFEST)
            or accepted["matrix_sha256"] != axis_run.digest(manifest)
            or accepted["preflight_sha256"] != axis_run.digest(preflight)
            or accepted["release_files_sha256"]["log/runs.csv"] != axis_run.digest(ledger)):
        raise ValueError("sealed Step 08 acceptance, matrix or run ledger differs")
    runs = axis_run.read_json(manifest)
    if runs["preflight_sha256"] != accepted["preflight_sha256"]:
        raise ValueError("sealed Step 08 matrix has a different registration plan")
    rows = {tuple(fit["coordinate"]): fit for fit in runs["fits"] if fit["coordinate"][0] == ARM}
    table = chem_runs.axis_evaluate._csv(ledger)
    table = table[table.arm.eq(ARM)].set_index(["arm", "fold", "seed"], verify_integrity=True)
    expected = {(ARM, fold, seed) for axis, fold, seed in coordinates() if axis == REFERENCE_AXIS}
    if set(rows) != expected or set(table.index) != expected:
        raise ValueError("sealed Step 08 release does not have exactly nine UniMol coordinates")
    return matrix, rows, table, axis_run.digest(acceptance)


def check_reference(key, data, vectors, source, matrix, rows, table):
    axis, fold, seed = key
    historical = ARM, fold, seed
    fit = rows[historical]
    record = table.loc[historical]
    relative = Path("runs") / ARM / f"fold-{fold}" / f"seed-{seed}"
    path = matrix / relative
    observed, (_, config, replay, training, _) = chem_runs.read_run(path, False)
    expected = chem_runs.bind(historical, data, vectors, source, False)
    current = {field: value for field, value in config.items() if field not in ("arguments", "nonce")}
    manifest_sha = axis_run.digest(path / "manifest.json")
    hashes = {"registration_sha256": replay["registration_sha256"], "manifest_sha256": manifest_sha,
              "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
              "prediction_sha256": replay["outputs"]["prediction.npy"],
              "trace_sha256": replay["outputs"]["learning_curve.csv"]}
    if (observed != historical or current != expected or fit["path"] != relative.as_posix()
            or any(record[field] != value for field, value in hashes.items())
            or any(fit[field] != hashes[field] for field in ("registration_sha256", "manifest_sha256"))
            or fit["initialization"] != training["initialization"]
            or fit["best_epoch"] != record["best_epoch"] or fit["stopping_epoch"] != record["stopping_epoch"]):
        raise ValueError("Step 08 checkpoint, registration, protocol or published provenance differs")
    adapter = AxisAdapter(source, config["model"]["transform"])
    model = adapter.build(seed, "cuda")
    initialized = chem_runs.initialization(model, ARM)
    if initialized != training["initialization"]:
        raise ValueError("Step 08 native initialization or training-start RNG differs")
    pairing = paired_initialization(model)
    # Identical fit implementation, sampler, partitions, dose and fixed initial native RNG
    # are bound above through the complete Step 08 config and historical receipts.
    model_data = chem_runs.ChemicalInputs(data, vectors, ARM)
    values, test = axis_run.predict_test(path, adapter, model, model_data, config,
                                        hashes["registration_sha256"], training)
    axis_run.record_predictions(path, None, values, test, model_data, config, replay=True)
    return {
        "coordinate": list(key), "decision": "reuse", **hashes,
        "validation_prediction_sha256": replay["outputs"]["validation_prediction.npy"],
        "initialization": initialized, "paired_initialization": pairing,
        "training_stream": "same inherited fit loop, ordered rows, rng seeds and native initialization",
        "validation_replay": "exact", "test_replay": "exact",
    }


def reference_plan(step08, step03, representations, source):
    """Return one explicit reuse or refit decision per registered reference coordinate."""
    chem_runs.verify_release()
    matrix, rows, table, acceptance_sha = verify_publication(step08)
    configure_cuda()
    decisions = []
    training = tenv_axis_output(REFERENCE_AXIS, step03)
    with chem_runs.open_inputs(training, representations) as (data, vectors):
        for key in coordinates():
            if key[0] != REFERENCE_AXIS:
                continue
            try:
                result = check_reference(key, data, vectors, source, matrix, rows, table)
            except DEFECTS as error:
                result = {"coordinate": list(key), "decision": "refit", "reason": str(error)}
            decisions.append(result)
    return {"step08_acceptance_sha256": acceptance_sha, "decisions": decisions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_inputs(parser)
    parser.add_argument("--step08", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(reference_plan(args.step08, args.step03, args.representations, args.source), sort_keys=True))


if __name__ == "__main__":
    main()
