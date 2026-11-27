"""Register, fit and independently replay a Step 10 dose coordinate."""

import argparse
from functools import partial
import json
import os
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import torch

from dose_inputs import features_for, open_inputs
from dose_model import DoseAdapter, configuration, epoch_streams, initialization, objective
from dose_protocol import (ARMS, AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_STEP03,
                           SCENARIO, SEEDS, coordinates, representation)
from chem_runs import ALPHA_MANIFEST, axis_run, verify_release
from tsg_axis.learning import EarlyStopping, fixed_rng
from tsg_axis.protocol import MAX_EPOCHS, PATIENCE
from util.training import registered
from util.training.provenance import verify_replay_configuration, verify_sealed_code
from util.training.transigen import fit_epochs, replay_validation

SCHEMA = "carps-vae-dose-run-v1"
ORIGINAL_CODE = "b8be1ad9244c9a3e03eba3e0c77ec64a59ed695d"
DOSE_CODE = tuple(
    Path(__file__).resolve().with_name(name).relative_to(axis_run.REPO).as_posix()
    for name in ("dose_runs.py", "dose_model.py", "dose_inputs.py", "dose_protocol.py")
)
PHASES = {"train": "training.json", "predict": "predictions.json", "replay": "manifest.json"}
TRAINING_FILES = axis_run.TRAINING_FILES
PREDICTION_FILES = axis_run.PREDICTION_FILES
PHASE_FILES = {"train": TRAINING_FILES, "predict": PREDICTION_FILES,
               "replay": (*TRAINING_FILES, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])}
TRACES = ("epoch", "train_loss", "validation_loss", "validation_mse")


def arguments(key, bounded):
    arm, fold, seed = key
    if type(fold) is not int or type(seed) is not int or key not in coordinates(bounded):
        raise ValueError(f"unregistered Step 10 coordinate: {key}")
    return {"arm": arm, "fold": fold, "seed": seed, "bounded": bounded, "axis": AXIS, "device": "cuda"}


def input_identity(data, fold, bounded):
    release = data.release
    positions = data.positions(fold, bounded)
    return {
        "training_sha256": release.manifest_sha256,
        "raw_sha256": release.manifest["inputs"]["raw_sha256"],
        "folds": release.manifest["inputs"]["folds"],
        "axis_sha256": release.manifest["outputs"]["genes.csv"],
        "population": data.fold_population(fold, bounded),
        "partitions": {phase: rows.tolist() for phase, rows in positions.items()},
        "conditions": {phase: list(data.condition_ids(rows)) for phase, rows in positions.items()},
    }


def code_hashes():
    """Bind the dose runner, its shared runtime and the alpha input adapters."""
    extra = (
        *DOSE_CODE,
        "stage2_CarPS-VAE/07_cpsv_alpha/output/api/carps_vae_alpha.py",
        "stage2_CarPS-VAE/07_cpsv_alpha/execution/verify_alpha.py",
        "stage2_CarPS-VAE/08_cpsv_chem-emb/execution/chem_runs.py",
        "stage2_CarPS-VAE/03_data_train/output/api/fold_ledger.py",
        "stage2_CarPS-VAE/03_data_train/output/api/tenv_axes.py",
        "stage2_CarPS-VAE/util/output/api/stage2_util/derivation/release.py",
        "stage2_CarPS-VAE/util/output/api/stage2_util/publication.py",
    )
    return {**axis_run.code_hashes(), **axis_run.file_hashes(axis_run.REPO, extra)}


def bind(key, data, vectors, source, bounded):
    args = arguments(key, bounded)
    arm, fold, seed = key
    chemistry = representation(arm)
    native = axis_run.file_hashes(source, ("src/model.py", "src/train_TranSiGen.py", "src/utils.py"))
    if (native["src/model.py"] != axis_run.kpgt.SOURCE_DIGEST
            or native["src/train_TranSiGen.py"] != axis_run.kpgt.SCRIPT_HASH):
        raise ValueError("TranSiGen author source differs from the frozen model and configuration")
    config = {
        **input_identity(data, fold, bounded),
        "schema": SCHEMA, "scenario": SCENARIO, "seed": seed,
        "alpha_release_sha256": axis_run.digest(ALPHA_MANIFEST),
        "native_sha256": native,
        "code_sha256": code_hashes(),
        "representation": {"name": chemistry, "width": vectors.width(chemistry),
                           "manifest_sha256": vectors.manifest_sha256,
                           "vectors_sha256": vectors.manifest["files"][f"{chemistry}.npy"]},
        "model": configuration(arm), "objective": objective(arm),
        "validation_seed": axis_run.VALIDATION_SEED, "test_seed": axis_run.TEST_SEED,
        "max_epochs": 2 if bounded else MAX_EPOCHS, "patience": PATIENCE,
        "environment": {"python": platform.python_version(), "torch": str(torch.__version__),
                        "numpy": np.__version__, "cuda": torch.version.cuda,
                        "torch_threads": torch.get_num_threads(), "hardware": torch.cuda.get_device_name()},
        "arguments": args,
    }
    conditioning = arm if config["model"]["transform"] is not None else None
    return axis_run.transigen.bind_experiment(args, config, dose_conditioning=conditioning)


def register(output, key, bounded, step03, representations, source):
    verify_release()

    def definition():
        with open_inputs(step03, representations) as (data, vectors):
            return bind(key, data, vectors, source, bounded)

    return registered.register(output, definition)


def evaluate(model, inputs, config, phase, rows):
    """Run the model with the registered independent evaluation stream."""
    observed = inputs.data.model_input(phase)
    parameter = next(model.parameters())
    batch_size = config["model"]["batch_size"]
    objective_sum = squared_error = 0.0
    predictions = []
    model.eval()
    seed = config["validation_seed"] if phase == "valid" else config["test_seed"]
    with fixed_rng(seed, parameter.device), torch.no_grad():
        for start in range(0, len(rows), batch_size):
            batch = rows[start:start + batch_size]
            control = torch.as_tensor(observed.control[batch], device=parameter.device)
            features = torch.as_tensor(inputs.fingerprints(batch), device=parameter.device)
            result = model(control, features)
            predicted = result[3]
            if phase == "valid":
                measured = torch.as_tensor(observed.expression[batch], device=parameter.device)
                loss = batch_objective(model, control, measured, result)[0]
                objective_sum += loss.item()
                squared_error += (predicted - measured).square().sum().item()
            predictions.append(predicted.cpu().numpy())
    values = np.concatenate(predictions)
    objective_value, mse = objective_sum / len(rows), squared_error / values.size
    if (values.dtype != np.float32 or not np.isfinite(values).all()
            or not np.isfinite([objective_value, mse]).all()):
        raise ValueError(f"nonfinite or non-float32 {phase} result")
    return values, objective_value, mse


def batch_objective(model, control, measured, result):
    reconstructed, mu1, logvar1, predicted, mu_pred, logvar_pred, _ = result
    latent, mu2, logvar2 = model.encode_x2(measured)
    measured_reconstructed = model.decode_x2(latent)
    return model.loss(control, reconstructed, mu1, logvar1, measured, measured_reconstructed,
                      mu2, logvar2, predicted, mu_pred, logvar_pred)


def train_epoch(adapter, model, inputs, rows, optimizer, rng):
    observed = inputs.data.model_input("train")
    model.train()
    order = rng.permutation(rows)
    total, count = 0.0, 0
    for start in range(0, len(order), adapter.batch_size):
        batch = order[start:start + adapter.batch_size]
        if len(batch) == 1:
            continue
        control, measured, features = (
            torch.as_tensor(values, device=adapter.device)
            for values in (observed.control[batch], observed.expression[batch], inputs.fingerprints(batch))
        )
        optimizer.zero_grad()
        loss = batch_objective(model, control, measured, model(control, features))[0]
        if not torch.isfinite(loss):
            raise ValueError("nonfinite model training objective")
        loss.backward()
        optimizer.step()
        total += loss.item()
        count += len(batch)
    if not count:
        raise ValueError("no non-singleton training batch")
    return total / count


def fit(root, adapter, model, inputs, config, registration_sha):
    train, valid = (np.asarray(config["partitions"][part]) for part in ("train", "valid"))
    optimizer = adapter.optimizer(model)
    rng = np.random.default_rng(config["seed"])
    streams = []

    def train_step():
        loss = train_epoch(adapter, model, inputs, train, optimizer, rng)
        streams.append(epoch_streams(adapter.device, rng))
        return {"train_loss": loss}

    def validate_step():
        values, objective, mse = evaluate(model, inputs, config, "valid", valid)
        return values, {"validation_loss": objective, "validation_mse": mse}

    details = fit_epochs(root, model, config, registration_sha, train_step, validate_step)
    return {**details, "epoch_streams": streams}


def predict_test(root, model, inputs, config, registration_sha, training):
    valid, test = (np.asarray(config["partitions"][part]) for part in ("valid", "test"))
    checkpoint_validation = partial(evaluate, model, inputs, config, "valid", valid)
    replay_validation(root, model, registration_sha, training, checkpoint_validation)
    predictions, _, _ = evaluate(model, inputs, config, "test", test)
    return predictions, test


def execute(output, phase, step03, representations, source):
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    axis_run.kpgt.configure_backend("cuda")
    verify_release()
    with registered.phase(output, phase, PHASE_FILES, separate_registration=True) as run:
        config = {name: value for name, value in run.config.items() if name != "nonce"}
        args = config["arguments"]
        key = args["arm"], args["fold"], args["seed"]
        with open_inputs(step03, representations) as (data, vectors):
            current = bind(key, data, vectors, source, args["bounded"])
            verify_replay_configuration(config, current, replay=run.sealed is not None,
                                        repo=axis_run.REPO, revision=ORIGINAL_CODE, historical_paths=DOSE_CODE)
            adapter = DoseAdapter(source, key[0])
            model = adapter.build(key[2], "cuda", n_genes=len(data.genes))
            inputs = features_for(data, vectors, key[0])
            if phase == "train":
                initial = initialization(model, key[0])
                details = fit(run.scratch, adapter, model, inputs, config, run.registration_sha)
                details["initialization"] = initial
            else:
                predictions, rows = predict_test(run.root, model, inputs, config, run.registration_sha, run.training)
                details = axis_run.record_predictions(run.root, run.scratch, predictions, rows, data, config,
                                                      replay=phase == "replay")
        return run.publish(details, {**args, "schema": SCHEMA, "objective": config["objective"]},
                           report_replay_resources=True)


def read_run(path, bounded):
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
    verify_sealed_code(config, axis_run.REPO, ORIGINAL_CODE, historical_paths=DOSE_CODE)
    args = config["arguments"]
    key = args["arm"], args["fold"], args["seed"]
    registration_sha = axis_run.digest(registration)
    if (args != arguments(key, bounded) or config["schema"] != SCHEMA or config["scenario"] != SCENARIO
            or config["seed"] != key[2] or config["objective"] != objective(key[0])
            or config["model"] != configuration(key[0])):
        raise ValueError(f"invalid sealed Step 10 registration: {path}")
    training, prediction, replay = (
        axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase]) for phase in PHASES
    )
    receipts = (training, prediction, replay)
    if (len({config["nonce"], *(receipt["process_token"] for receipt in receipts)}) != 4
            or any(record.get(field) != expected for record in receipts
                   for field, expected in {**args, "schema": SCHEMA, "objective": config["objective"]}.items())
            or any(record["cuda_peak_allocated_bytes"] <= 0 for record in receipts)
            or any(replay.get(field) != "exact" for field in
                   ("validation_replay", "validation_objective_replay", "test_replay"))
            or any(replay["outputs"][name] != digest for record in receipts[:2]
                   for name, digest in record["outputs"].items())):
        raise ValueError(f"invalid Step 10 independent phase evidence: {path}")
    trace = pd.read_csv(path / "learning_curve.csv", float_precision="round_trip")
    if (tuple(trace.columns) != TRACES or trace.empty or not np.isfinite(trace.to_numpy()).all()
            or trace.epoch.to_list() != list(range(1, len(trace) + 1))
            or len(trace) != training["epochs"] or len(trace) != len(training["epoch_streams"])
            or int(trace.validation_loss.idxmin()) + 1 != training["best_epoch"]):
        raise ValueError(f"invalid Step 10 learning trace: {path}")
    best = trace.iloc[training["best_epoch"] - 1]
    if any(best[name] != training[name] for name in ("validation_loss", "validation_mse")):
        raise ValueError(f"Step 10 trace differs from its minimum checkpoint: {path}")
    stopper = EarlyStopping(patience=config["patience"], max_epochs=config["max_epochs"])
    for epoch, value in zip(trace.epoch, trace.validation_loss, strict=True):
        stopper.observe(value, epoch)
    if not stopper.stopped:
        raise ValueError(f"Step 10 trace ended before the registered stopping condition: {path}")
    return key, (path, config, replay, training, trace)


def add_inputs(parser):
    for name, default in (("step03", DEFAULT_STEP03), ("representations", DEFAULT_REPRESENTATIONS),
                          ("source", DEFAULT_SOURCE)):
        parser.add_argument(f"--{name}", type=Path, default=default)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--arm", choices=ARMS, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--seed", choices=SEEDS, type=int, required=True)
    registration.add_argument("--bounded", action="store_true")
    for phase in PHASES:
        commands.add_parser(phase)
    for command in commands.choices.values():
        add_inputs(command)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = tuple(path.resolve() for path in (args.step03, args.representations, args.source))
    if args.phase == "register":
        result = {"registration_sha256": register(args.output, (args.arm, args.fold, args.seed),
                                                   args.bounded, *inputs)}
    else:
        result = execute(args.output, args.phase, *inputs)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
