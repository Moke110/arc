"""Register, fit and independently replay a Step 10 dose coordinate."""

import argparse
import json
import os
from pathlib import Path
import platform
import resource
from tempfile import TemporaryDirectory
import time

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

SCHEMA = "carps-vae-dose-run-v1"
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


def bind(key, data, vectors, source, bounded):
    args = arguments(key, bounded)
    arm, fold, seed = key
    chemistry = representation(arm)
    native = axis_run.file_hashes(source, ("src/model.py", "src/train_TranSiGen.py", "src/utils.py"))
    if (native["src/model.py"] != axis_run.kpgt.SOURCE_DIGEST
            or native["src/train_TranSiGen.py"] != axis_run.kpgt.SCRIPT_HASH):
        raise ValueError("TranSiGen author source differs from the frozen model and configuration")
    step = Path(__file__).resolve().parent
    dependencies = tuple(step / name for name in ("dose_runs.py", "dose_model.py", "dose_inputs.py", "dose_protocol.py"))
    return {
        **input_identity(data, fold, bounded),
        "schema": SCHEMA, "scenario": SCENARIO, "seed": seed,
        "alpha_release_sha256": axis_run.digest(ALPHA_MANIFEST),
        "native_sha256": native,
        "code_sha256": {path.relative_to(axis_run.REPO).as_posix(): axis_run.digest(path) for path in dependencies},
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


def register(output, key, bounded, step03, representations, source):
    verify_release()
    with axis_run.publication_lock(output) as root:
        if root.exists():
            raise ValueError("run registration requires an unused output directory")
        with open_inputs(step03, representations) as (data, vectors):
            config = bind(key, data, vectors, source, bounded)
        root.mkdir()
        registration = root / "registration.json"
        axis_run.write_json(registration, {**config, "nonce": axis_run.PROCESS})
        return axis_run.digest(registration)


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
    stopper = EarlyStopping(patience=config["patience"], max_epochs=config["max_epochs"])
    optimizer = adapter.optimizer(model)
    rng = np.random.default_rng(config["seed"])
    history, streams = [], []
    for epoch in range(1, config["max_epochs"] + 1):
        loss = train_epoch(adapter, model, inputs, train, optimizer, rng)
        streams.append(epoch_streams(adapter.device, rng))
        values, validation, mse = evaluate(model, inputs, config, "valid", valid)
        if stopper.observe(validation, epoch):
            torch.save({"registration_sha256": registration_sha, "state_dict": model.state_dict(),
                        "epoch": epoch, "validation_loss": validation, "validation_mse": mse},
                       root / "checkpoint.pt")
            np.save(root / "validation_prediction.npy", values)
        history.append((epoch, loss, validation, mse))
        if stopper.stopped:
            break
    pd.DataFrame(history, columns=TRACES).to_csv(root / "learning_curve.csv", index=False)
    return {"best_epoch": stopper.best_epoch, "validation_loss": stopper.best,
            "validation_mse": history[stopper.best_epoch - 1][-1], "epochs": len(history),
            "parameters": sum(p.numel() for p in model.parameters()), "epoch_streams": streams}


def predict_test(root, model, inputs, config, registration_sha, training):
    checkpoint = torch.load(root / "checkpoint.pt", map_location="cpu", weights_only=True)
    metrics = ("validation_loss", "validation_mse")
    if (checkpoint["registration_sha256"] != registration_sha or checkpoint["epoch"] != training["best_epoch"]
            or any(checkpoint[key] != training[key] for key in metrics)):
        raise ValueError("checkpoint differs from its minimum-validation-objective receipt")
    model.load_state_dict(checkpoint["state_dict"])
    del checkpoint
    valid, test = (np.asarray(config["partitions"][part]) for part in ("valid", "test"))
    validation, value, mse = evaluate(model, inputs, config, "valid", valid)
    if (value, mse) != tuple(training[key] for key in metrics):
        raise ValueError("fresh process did not exactly replay the validation objective")
    if not np.array_equal(validation, np.load(root / "validation_prediction.npy", allow_pickle=False)):
        raise ValueError("fresh process did not exactly replay validation predictions")
    predictions, _, _ = evaluate(model, inputs, config, "test", test)
    return predictions, test


def execute(output, phase, step03, representations, source):
    started = time.monotonic()
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    axis_run.kpgt.configure_backend("cuda")
    verify_release()
    with axis_run.publication_lock(output) as root:
        registration = root / "registration.json"
        config = axis_run.read_json(registration)
        nonce = config.pop("nonce")
        registration_sha = axis_run.digest(registration)
        args = config["arguments"]
        key = args["arm"], args["fold"], args["seed"]
        if nonce == axis_run.PROCESS:
            raise ValueError("run phase must be separate from registration")
        receipt = root / PHASES[phase]
        if receipt.exists() and phase != "replay":
            raise ValueError(f"run phase already recorded: {phase}")
        sealed = (axis_run.verified_receipt(root, receipt.name, registration_sha, PHASE_FILES[phase])
                  if receipt.exists() else None)
        trained = (axis_run.verified_receipt(root, PHASES["train"], registration_sha, TRAINING_FILES)
                   if phase != "train" else None)
        if phase == "replay":
            axis_run.verified_receipt(root, PHASES["predict"], registration_sha, PREDICTION_FILES)
        with open_inputs(step03, representations) as (data, vectors):
            if config != bind(key, data, vectors, source, args["bounded"]):
                raise ValueError("registered Step 10 study, inputs, code or runtime changed")
            adapter = DoseAdapter(source, key[0])
            model = adapter.build(key[2], "cuda", n_genes=len(data.genes))
            inputs = features_for(data, vectors, key[0])
            with TemporaryDirectory(prefix=".phase-", dir=root) as temporary:
                scratch = Path(temporary)
                if phase == "train":
                    initial = initialization(model, key[0])
                    details = fit(scratch, adapter, model, inputs, config, registration_sha)
                    details["initialization"] = initial
                else:
                    predictions, test = predict_test(root, model, inputs, config, registration_sha, trained)
                    details = axis_run.record_predictions(root, scratch, predictions, test, data, config,
                                                          replay=phase == "replay")
                files = PHASE_FILES[phase]
                if phase == "replay":
                    details["test_replay"] = "exact"
                    hashes = axis_run.file_hashes(root, files)
                else:
                    hashes = axis_run.file_hashes(scratch, files)
                    for name in files:
                        (scratch / name).replace(root / name)
                record = {
                    **details, **args, "schema": SCHEMA, "objective": config["objective"], "outputs": hashes,
                    "registration_sha256": registration_sha, "process_token": axis_run.PROCESS,
                    "pid": os.getpid(), "wall_seconds": round(time.monotonic() - started, 3),
                    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    "cuda_peak_allocated_bytes": torch.cuda.max_memory_allocated(),
                }
                if sealed is not None:
                    return record
                axis_run.write_json(scratch / receipt.name, record)
                (scratch / receipt.name).replace(receipt)
    return record


def read_run(path, bounded):
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
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
