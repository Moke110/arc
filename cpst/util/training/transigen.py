"""Fit and independently replay a registered matched-control coordinate."""
from __future__ import annotations

import os
from pathlib import Path
import resource
import secrets
from tempfile import TemporaryDirectory
import time

import numpy as np
import pandas as pd
import torch

from cpst_util.derivation.release import check_outputs, read_json, sha256 as digest, write_json
from cpst_util.publication import publication_lock
from cpst_util.release_names import checkpoint_registration_sha256, original_axis_manifest_sha256
from tsg_axis import model as kpgt
from tsg_axis.inputs import open_inputs, score_fixed_prediction
from tsg_axis.learning import EarlyStopping, fixed_rng
from tsg_axis.model import DoseView

from . import experiment_configuration
from .additive_objective import ComponentTrace

PROCESS = secrets.token_hex(16)
INPUT_PATHS = ("training", "representations", "source")
TRAINING_FILES = ("checkpoint.pt", "validation_prediction.npy", "learning_curve.csv")
PREDICTION_FILES = ("prediction.npy", "conditions.csv", "genes.csv", "predicted_delta.npy", "measured_delta.npy")
PHASES = {"train": "training.json", "predict": "predictions.json", "replay": "manifest.json"}


def file_hashes(root, names):
    return {name: digest(root / name) for name in names}


def prediction(native, model, data, phase, rows, seed):
    with fixed_rng(seed, next(model.parameters()).device):
        result = native.predict(model, data.model_input(phase), rows)
    if result.dtype != np.float32 or not np.isfinite(result).all():
        raise ValueError("nonfinite or non-float32 model prediction")
    return result


def validate(model, data, rows, config):
    """Replay the author validation objective and its sampled predictions in one pass."""
    inputs = DoseView(data.model_input("valid"), config["model"]["transform"])
    batch_size = config["model"]["batch_size"]
    parameter = next(model.parameters())
    model.eval()
    total = prediction_squared_error = 0.0
    predictions = []
    with fixed_rng(config["validation_seed"], parameter.device), torch.no_grad():
        for start in range(0, len(rows), batch_size):
            batch = rows[start:start + batch_size]
            control, measured, features = (
                torch.as_tensor(values, dtype=parameter.dtype, device=parameter.device)
                for values in (inputs.control[batch], inputs.expression[batch], inputs.fingerprints(batch))
            )
            reconstructed, mu1, logvar1, predicted, mu_pred, logvar_pred, _ = model(control, features)
            latent, mu2, logvar2 = model.encode_x2(measured)
            measured_reconstructed = model.decode_x2(latent)
            loss = model.loss(control, reconstructed, mu1, logvar1, measured, measured_reconstructed,
                              mu2, logvar2, predicted, mu_pred, logvar_pred)[0]
            total += loss.item()
            prediction_squared_error += (predicted - measured).square().sum().item()
            predictions.append(predicted.cpu().numpy())
    values = np.concatenate(predictions)
    objective, mse = total / len(rows), prediction_squared_error / values.size
    if values.dtype != np.float32 or not np.isfinite(values).all() or not np.isfinite([objective, mse]).all():
        raise ValueError("nonfinite or non-float32 native validation result")
    return values, objective, mse


def fit_epochs(root, model, config, registration_sha, train_epoch, validate_epoch):
    """Fit a registered coordinate and retain its minimum-objective checkpoint."""
    stopper = EarlyStopping(patience=config["patience"], max_epochs=config["max_epochs"])
    history = []
    for epoch in range(1, config["max_epochs"] + 1):
        train = train_epoch()
        values, validation = validate_epoch()
        objective, mse = validation["validation_loss"], validation["validation_mse"]
        if stopper.observe(objective, epoch):
            torch.save({"registration_sha256": registration_sha, "state_dict": model.state_dict(),
                        "epoch": epoch, "validation_loss": objective, "validation_mse": mse}, root / "checkpoint.pt")
            np.save(root / "validation_prediction.npy", values)
        # Keep the sealed trace's base columns ahead of study-specific components.
        history.append({"epoch": epoch, "train_loss": train["train_loss"],
                        "validation_loss": objective, "validation_mse": mse, **train, **validation})
        if stopper.stopped:
            break
    pd.DataFrame(history).to_csv(root / "learning_curve.csv", index=False)
    return {"best_epoch": stopper.best_epoch, "validation_loss": stopper.best,
            "validation_mse": history[stopper.best_epoch - 1]["validation_mse"], "epochs": len(history),
            "parameters": sum(parameter.numel() for parameter in model.parameters())}


def fit(root, native, model, data, config, registration_sha, *, trace_factory=None):
    training = np.asarray(config["partitions"]["train"])
    valid = np.asarray(config["partitions"]["valid"])
    training_inputs = data.model_input("train")
    rng = np.random.default_rng(config["seed"])
    optimizer = native.optimizer(model)
    trace = trace_factory(model.loss) if trace_factory is not None else None
    if trace is not None:
        model.loss = trace

    def train_epoch():
        loss = native.train_epoch(model, training_inputs, training, optimizer, rng)
        if trace is None and not np.isfinite(loss):
            raise ValueError("nonfinite native training loss")
        return {"train_loss": loss, **(trace.per_condition("train") if trace is not None else {})}

    def validate_epoch():
        values, objective, mse = validate(model, data, valid, config)
        components = trace.per_condition("validation") if trace is not None else {}
        return values, {"validation_loss": objective, "validation_mse": mse, **components}

    return fit_epochs(root, model, config, registration_sha, train_epoch, validate_epoch)


def replay_validation(root, model, registration_sha, training, validate_checkpoint):
    """Restore the selected state and independently replay the validation pass."""
    checkpoint = torch.load(root / "checkpoint.pt", map_location="cpu", weights_only=True)
    metrics = ("validation_loss", "validation_mse")
    if (checkpoint["registration_sha256"] != checkpoint_registration_sha256(registration_sha)
            or checkpoint["epoch"] != training["best_epoch"]
            or any(checkpoint[key] != training[key] for key in metrics)):
        raise ValueError("validation checkpoint identity differs from its training receipt")
    model.load_state_dict(checkpoint["state_dict"])
    del checkpoint
    validation, objective, mse = validate_checkpoint()
    if (objective, mse) != tuple(training[key] for key in metrics):
        raise ValueError("fresh process did not exactly replay the validation objective")
    if not np.array_equal(validation, np.load(root / "validation_prediction.npy", allow_pickle=False)):
        raise ValueError("fresh process did not exactly replay validation predictions")


def kl_diagnostics(trace, latent):
    diagnostics = {"epoch": trace.epoch}
    for partition in ("train", "validation"):
        total = (trace[f"{partition}_control_kl"]
                 + trace[f"{partition}_treatment_kl"]
                 + trace[f"{partition}_prediction_kl"])
        diagnostics[f"{partition}_total_kl"] = total
        diagnostics[f"{partition}_kl_per_latent_dimension"] = total / latent
    return pd.DataFrame(diagnostics)


def fit_with_diagnostics(root, native, model, data, config, registration_sha):
    details = fit(root, native, model, data, config, registration_sha, trace_factory=ComponentTrace)
    trace = pd.read_csv(root / "learning_curve.csv", float_precision="round_trip")
    kl_diagnostics(trace, config["model"]["n_latent"]).to_csv(root / "kl-diagnostics.csv", index=False)
    return details


def predict_test(root, native, model, data, config, registration_sha, training):
    valid = np.asarray(config["partitions"]["valid"])
    test = np.asarray(config["partitions"]["test"])
    replay_validation(root, model, registration_sha, training, lambda: validate(model, data, valid, config))
    values = prediction(native, model, data, "test", test, config["test_seed"])
    return values, test


def record_predictions(root, output, values, test, data, config, *, replay):
    result = score_fixed_prediction(data, config["arguments"]["fold"], values, config["conditions"]["test"], rows=test)
    arrays = {"prediction": values, "predicted_delta": result.predicted_delta, "measured_delta": result.measured_delta}
    tables = {"conditions": data.condition_table(test), "genes": data.genes}
    for name, array in arrays.items():
        if replay:
            if not np.array_equal(array, np.load(root / f"{name}.npy", allow_pickle=False)):
                raise ValueError(f"fresh-process replay differs for {name}")
        else:
            np.save(output / f"{name}.npy", array)
    for name, table in tables.items():
        csv = table.to_csv(index=False)
        if replay:
            if (root / f"{name}.csv").read_text() != csv:
                raise ValueError(f"fresh-process replay differs for {name}")
        else:
            (output / f"{name}.csv").write_text(csv)
    return {"test_conditions": len(test), "validation_replay": "exact", "validation_objective_replay": "exact"}


def verified_receipt(root, name, registration_sha, files, *, output_check=check_outputs):
    record = read_json(root / name)
    if (record["registration_sha256"] != registration_sha or record["process_token"] == PROCESS
            or set(record["outputs"]) != set(files)):
        raise ValueError(f"run receipt changed or is not from a separate process: {name}")
    output_check(root, record)
    return record


def bind_experiment(args, config, **settings):
    config["experiment"] = experiment_configuration(args, config, **settings)
    return config


def register_run(output, args, bind):
    with publication_lock(output) as root:
        if root.exists():
            raise ValueError("run registration requires an unused output directory")
        with open_inputs(args["training"], args["representations"]) as (data, vectors):
            identity = bind(args, data, vectors)
        root.mkdir()
        record = root / "registration.json"
        write_json(record, {**identity, "arguments": args, "nonce": PROCESS})
        return digest(record)


def execute_phase(output, phase, paths, *, bind, build_model, fit,
                  training_files=TRAINING_FILES, record_cuda=False, expected_arguments=None):
    start = time.monotonic()
    with publication_lock(output) as root:
        registration = root / "registration.json"
        config = read_json(registration)
        registration_sha = digest(registration)
        if set(paths) - set(INPUT_PATHS):
            raise ValueError("only verified input locations can be overridden")
        args = {**config["arguments"], **paths}
        if expected_arguments and any(args[key] != value for key, value in expected_arguments.items()):
            raise ValueError("run registration belongs to a different similarity study")
        receipt = PHASES[phase]
        sealed_files = (*training_files, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])
        sealed = None
        if (root / receipt).exists():
            if phase != "replay":
                raise ValueError(f"run phase already recorded: {phase}")
            sealed = verified_receipt(root, receipt, registration_sha, sealed_files)
        if phase != "train":
            training = verified_receipt(root, PHASES["train"], registration_sha, training_files)
        if phase == "replay":
            verified_receipt(root, PHASES["predict"], registration_sha, PREDICTION_FILES)
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        with open_inputs(args["training"], args["representations"]) as (data, vectors):
            identity = bind(args, data, vectors)
            saved_identity = {key: value for key, value in config.items()
                              if key not in {"arguments", "nonce"}}
            axis = "axis_family_manifest_sha256"
            if sealed is not None and identity[axis] != saved_identity[axis]:
                identity[axis] = original_axis_manifest_sha256(identity[axis])
            if identity != saved_identity:
                raise ValueError("registered experiment, inputs, code or runtime changed")
            torch.use_deterministic_algorithms(True)
            kpgt.configure_backend(args["device"])
            native, model = build_model(args, config, data)
            with TemporaryDirectory(prefix=".phase-", dir=root) as scratch_path:
                scratch = Path(scratch_path)
                if phase == "train":
                    details = fit(scratch, native, model, data, config, registration_sha)
                    files = training_files
                else:
                    values, test = predict_test(root, native, model, data, config, registration_sha, training)
                    details = record_predictions(root, scratch, values, test, data, config, replay=phase == "replay")
                    files = PREDICTION_FILES
                if sealed is not None:
                    return sealed
                if phase == "replay":
                    details.update(test_replay="exact", axis=args["axis"], fold=args["fold"], smoke=args["smoke"],
                                   seed=config["seed"], schema=config["schema"], protocol=config["protocol"])
                    details.update({key: config[key] for key in ("shape", "objective") if key in config})
                    hashes = file_hashes(root, sealed_files)
                else:
                    hashes = file_hashes(scratch, files)
                    for name in files:
                        (scratch / name).replace(root / name)
                record = {
                    **details, "outputs": hashes, "registration_sha256": registration_sha,
                    "process_token": PROCESS, "pid": os.getpid(),
                    "wall_seconds": round(time.monotonic() - start, 3),
                    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                }
                if record_cuda:
                    record["cuda_peak_allocated_bytes"] = (
                        torch.cuda.max_memory_allocated() if args["device"] == "cuda" else None
                    )
                write_json(scratch / receipt, record)
                (scratch / receipt).replace(root / receipt)
    return record
