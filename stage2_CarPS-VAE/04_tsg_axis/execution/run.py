"""Register, fit and exactly replay one cold-drug TranSiGen axis coordinate."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import secrets
import sys
from tempfile import TemporaryDirectory
import time

import numpy as np
import pandas as pd
import torch

STEP = Path(__file__).resolve().parents[1]
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
STAGE = REPO / "stage2_CarPS-VAE"
AXIS_API = STAGE / "02_gene_axis/02.1_hvg-deg/output/api"
sys.path[:0] = [str(STEP / "output/api"), str(AXIS_API)]

from prediction_axis import axis_coordinates, open_axes
from tsg_axis.inputs import open_inputs, score_fixed_prediction
from tsg_axis.protocol import MAX_EPOCHS, PATIENCE, PROTOCOL, RUN_SCHEMA
from tsg_axis import model as kpgt
from tsg_axis.model import DoseView
from tsg_axis.learning import EarlyStopping, fixed_rng
from stage2_util.derivation.release import check_outputs, read_json, sha256 as digest, write_json
from stage2_util.publication import publication_lock
from training_data import TENV_AXIS_STEPS

SEED = 42
VALIDATION_SEED = 22617
TEST_SEED = 100042
PARTITIONS = ("train", "valid", "test")
SMOKE_LIMITS = dict(zip(PARTITIONS, (64, 24, 24), strict=True))
TRAINING_FILES = ("checkpoint.pt", "validation_prediction.npy", "learning_curve.csv")
PREDICTION_FILES = ("prediction.npy", "conditions.csv", "genes.csv", "predicted_delta.npy", "measured_delta.npy")
PROCESS = secrets.token_hex(16)


def file_hashes(root, names):
    return {name: digest(root / name) for name in names}


def code_hashes():
    api = STEP / "output/api/tsg_axis"
    paths = [
        Path(__file__),
        AXIS_API / "prediction_axis.py",
        *(api / name for name in ("inputs.py", "protocol.py", "model.py", "native_import.py", "learning.py")),
        STAGE / "03_data_train/output/api/training_data.py",
        STAGE / "util/output/api/stage2_util/representations.py",
    ]
    return {path.relative_to(REPO).as_posix(): digest(path) for path in paths}


def bind(args, data, vectors):
    family, width = axis_coordinates(args["axis"])
    axes = open_axes(family=family)
    axis_sha = axes.axis(width).sha256
    if (axis_sha != data.release.manifest["inputs"]["axis_sha256"]
            or axis_sha != data.release.manifest["outputs"]["genes.csv"] or len(data.genes) != width):
        raise ValueError("training release lacks the registered ordered prediction axis")
    rows = {name: data.rows(args["fold"], name) for name in PARTITIONS}
    parents = {name: set(data.parents[part]) for name, part in rows.items()}
    if (any(not part.size for part in rows.values())
            or sum(map(len, parents.values())) != len(set().union(*parents.values()))):
        raise ValueError("cold-drug partitions require disjoint nonempty parent populations")
    if args["smoke"]:
        rows = {name: part[:SMOKE_LIMITS[name]] for name, part in rows.items()}
    transform = kpgt.fit_transform(data.doses[rows["train"]])
    model = kpgt.configuration(transform)
    model.update(n_genes=width, method_identity="tenv-axis-transigen-kpgt-latent-numeric-v2",
                 control_input="matched plate DMSO for every partition",
                 validation_objective="native TranSiGen MSE + beta-weighted KLD, batch sums per condition",
                 adaptations=[f"{width} gene interface" if item == "5000 gene interface" else item
                              for item in model["adaptations"]
                              if item != "shared splits and fixed validation treated-expression MSE"]
                 + ["shared cold-drug splits and fixed validation native MSE + KLD"])
    model.update(batch_size=64, predictive_head_dim=532, dose_vector_dim=32)
    native = file_hashes(Path(args["source"]), ("src/model.py", "src/train_TranSiGen.py", "src/utils.py"))
    if (native["src/model.py"] != kpgt.SOURCE_DIGEST
            or native["src/train_TranSiGen.py"] != model["config_script_sha256"]):
        raise ValueError("registered native TranSiGen source differs from the pinned author release")
    hardware = (torch.cuda.get_device_name() if args["device"] == "cuda" else platform.machine())
    return {"schema": RUN_SCHEMA, "protocol": PROTOCOL, "axis_sha256": axis_sha,
            "axis_family_manifest_sha256": axes.manifest_sha256,
            "training_sha256": data.release.manifest_sha256,
            "raw_sha256": data.release.manifest["inputs"]["raw_sha256"],
            "folds": data.release.manifest["inputs"]["folds"],
            "kpgt_manifest_sha256": vectors.manifest_sha256, "kpgt_vectors_sha256": vectors.manifest["files"]["KPGT.npy"],
            "native_sha256": native,
            "code_sha256": code_hashes(), "model": model, "partitions": {key: part.tolist() for key, part in rows.items()},
            "conditions": {key: list(data.condition_ids(part)) for key, part in rows.items()},
            "environment": {"python": platform.python_version(), "torch": str(torch.__version__),
                            "numpy": np.__version__, "cuda": torch.version.cuda,
                            "torch_threads": torch.get_num_threads(), "hardware": hardware},
            "scenario": "cold-drug", "seed": SEED,
            "validation_seed": VALIDATION_SEED, "test_seed": TEST_SEED,
            "max_epochs": 2 if args["smoke"] else MAX_EPOCHS, "patience": PATIENCE}


def register(output, args):
    with publication_lock(output) as root:
        if root.exists():
            raise ValueError("run registration requires an unused output directory")
        with open_inputs(args["training"], args["representations"]) as (data, vectors):
            config = bind(args, data, vectors)
        root.mkdir()
        write_json(root / "registration.json", {**config, "arguments": args, "nonce": PROCESS})


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


def fit(root, native, model, data, config, registration_sha):
    training = np.asarray(config["partitions"]["train"])
    valid = np.asarray(config["partitions"]["valid"])
    training_inputs = data.model_input("train")
    rng = np.random.default_rng(config["seed"])
    optimizer = native.optimizer(model)
    stopper = EarlyStopping(patience=config["patience"], max_epochs=config["max_epochs"])
    history = []
    for epoch in range(1, config["max_epochs"] + 1):
        loss = native.train_epoch(model, training_inputs, training, optimizer, rng)
        values, objective, mse = validate(model, data, valid, config)
        if not np.isfinite(loss):
            raise ValueError("nonfinite native training loss")
        if stopper.observe(objective, epoch):
            torch.save({"registration_sha256": registration_sha, "state_dict": model.state_dict(),
                        "epoch": epoch, "validation_loss": objective, "validation_mse": mse}, root / "checkpoint.pt")
            np.save(root / "validation_prediction.npy", values)
        history.append({"epoch": epoch, "train_loss": loss, "validation_loss": objective, "validation_mse": mse})
        if stopper.stopped:
            break
    pd.DataFrame(history).to_csv(root / "learning_curve.csv", index=False)
    return {"best_epoch": stopper.best_epoch, "validation_loss": stopper.best,
            "validation_mse": history[stopper.best_epoch - 1]["validation_mse"], "epochs": len(history),
            "parameters": sum(p.numel() for p in model.parameters())}


def predict_test(root, native, model, data, config, registration_sha, training):
    checkpoint = torch.load(root / "checkpoint.pt", map_location="cpu", weights_only=True)
    metrics = ("validation_loss", "validation_mse")
    if (checkpoint["registration_sha256"] != registration_sha or checkpoint["epoch"] != training["best_epoch"]
            or any(checkpoint[key] != training[key] for key in metrics)):
        raise ValueError("validation checkpoint identity differs from its training receipt")
    model.load_state_dict(checkpoint["state_dict"])
    del checkpoint
    valid = np.asarray(config["partitions"]["valid"])
    test = np.asarray(config["partitions"]["test"])
    validation, objective, mse = validate(model, data, valid, config)
    if (objective, mse) != tuple(training[key] for key in metrics):
        raise ValueError("fresh process did not exactly replay the validation objective")
    if not np.array_equal(validation, np.load(root / "validation_prediction.npy", allow_pickle=False)):
        raise ValueError("fresh process did not exactly replay validation predictions")
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


def verified_receipt(root, name, registration_sha, files):
    record = read_json(root / name)
    if (record["registration_sha256"] != registration_sha or record["process_token"] == PROCESS
            or set(record["outputs"]) != set(files)):
        raise ValueError(f"run receipt changed or is not from a separate process: {name}")
    check_outputs(root, record)
    return record


def execute(output, phase, paths):
    start = time.monotonic()
    with publication_lock(output) as root:
        path = root / "registration.json"
        config, registration_sha = read_json(path), digest(path)
        args = {**config["arguments"], **paths}
        output_receipt = {"train": "training.json", "predict": "predictions.json", "replay": "manifest.json"}[phase]
        already_recorded = (root / output_receipt).exists()
        if already_recorded and phase != "replay":
            raise ValueError(f"run phase already recorded: {phase}")
        sealed_files = (*TRAINING_FILES, *PREDICTION_FILES, "training.json", "predictions.json")
        sealed = verified_receipt(root, output_receipt, registration_sha, sealed_files) if already_recorded else None
        with open_inputs(args["training"], args["representations"]) as (data, vectors):
            if bind(args, data, vectors) != {key: value for key, value in config.items()
                                           if key not in {"arguments", "nonce"}}:
                raise ValueError("registered experiment, provenance or runtime changed")
            if phase != "train":
                training = verified_receipt(root, "training.json", registration_sha, TRAINING_FILES)
            if phase == "replay":
                verified_receipt(root, "predictions.json", registration_sha, PREDICTION_FILES)
            os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
            torch.use_deterministic_algorithms(True)
            kpgt.configure_backend(args["device"])
            native = kpgt.adapter(args["source"], config["model"]["transform"])
            model = native.build(SEED, args["device"], n_genes=len(data.genes))
            with TemporaryDirectory(prefix=".phase-", dir=root) as scratch:
                scratch = Path(scratch)
                if phase == "train":
                    details = fit(scratch, native, model, data, config, registration_sha)
                    files = TRAINING_FILES
                else:
                    values, test = predict_test(root, native, model, data, config, registration_sha, training)
                    details = record_predictions(root, scratch, values, test, data, config, replay=phase == "replay")
                    files = PREDICTION_FILES
                if sealed is not None:
                    return sealed
                if phase == "replay":
                    details.update(test_replay="exact", axis=args["axis"], fold=args["fold"],
                                   smoke=args["smoke"], seed=SEED, schema=RUN_SCHEMA, protocol=PROTOCOL)
                    hashes = file_hashes(root, sealed_files)
                else:
                    hashes = file_hashes(scratch, files)
                    for name in files:
                        (scratch / name).replace(root / name)
                write_json(scratch / output_receipt, {
                    **details, "outputs": hashes, "registration_sha256": registration_sha,
                    "process_token": PROCESS, "pid": os.getpid(),
                    "wall_seconds": round(time.monotonic() - start, 3),
                    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                })
                (scratch / output_receipt).replace(root / output_receipt)
    return read_json(output / output_receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--axis", choices=TENV_AXIS_STEPS, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--training", type=Path, required=True)
    registration.add_argument("--representations", type=Path, default=STAGE / "util/output/representations")
    registration.add_argument("--source", type=Path, default=REPO / "pert_model/TranSiGen/source")
    registration.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    registration.add_argument("--smoke", action="store_true")
    for phase in ("train", "predict", "replay"):
        command = commands.add_parser(phase)
        for field in ("training", "representations", "source"):
            command.add_argument(f"--{field}", type=Path, help="same pinned input at a new location")
    for child in commands.choices.values():
        child.add_argument("--output", type=Path, required=True)
    args = vars(parser.parse_args())
    phase, output = args.pop("phase"), args.pop("output").resolve()
    if phase == "register":
        register(output, {key: str(value.resolve()) if isinstance(value, Path) else value
                          for key, value in args.items()})
        print(json.dumps({"registration_sha256": digest(output / "registration.json")}))
    else:
        paths = {key: str(value.resolve()) for key, value in args.items() if value is not None}
        print(json.dumps(execute(output, phase, paths), sort_keys=True))


if __name__ == "__main__":
    main()
