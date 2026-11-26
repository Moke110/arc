"""Register, fit and exactly replay one Step 05.1 TranSiGen loss coordinate."""

import argparse
from fractions import Fraction
from functools import partial
import importlib.util
import json
import os
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

import numpy as np
import pandas as pd
import torch

from objective import ALPHAS, FORMS, ComponentTrace, Objective

STEP = Path(__file__).resolve().parents[1]
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
STAGE = REPO / "stage2_CarPS-VAE"
AXIS_STEP = STAGE / "04_tsg_axis"
AXIS = "denoise-5k"
SCHEMA = "stage2-tenv-transigen-loss-run-v1"
PROTOCOL = "matched-plate-control-delta-loss-v1"
INPUT_PATHS = ("training", "representations", "source")
PHASES = {"train": "training.json", "predict": "predictions.json", "replay": "manifest.json"}


def load_step04(name):
    spec = importlib.util.spec_from_file_location(f"step04_{name}", AXIS_STEP / "execution" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


axis_run = load_step04("run")
SEALED_FILES = (*axis_run.TRAINING_FILES, *axis_run.PREDICTION_FILES, PHASES["train"], PHASES["predict"])
DEFAULT_TRAINING = STAGE / "03_data_train/03.7_tenv_denoise-5k/output"
DEFAULT_REPRESENTATIONS = STAGE / "util/output/representations"
DEFAULT_SOURCE = axis_run.REPO / "pert_model/TranSiGen/source"


def bind(args, data, vectors):
    if args["axis"] != AXIS:
        raise ValueError("the loss comparison requires the frozen denoise-5k axis")
    objective = Objective(args["loss"], Fraction(args["alpha"]))
    config = axis_run.bind(args, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, objective=objective.identity())
    config["model"]["validation_objective"] = "registered Step 05.1 objective, batch sums per condition"
    config["model"]["adaptations"][-1] = "shared cold-drug splits and fixed validation using the registered objective"
    for path in (Path(__file__), STEP / "execution/objective.py"):
        config["code_sha256"][path.relative_to(axis_run.REPO).as_posix()] = axis_run.digest(path)
    return config


def register(output, args):
    with axis_run.publication_lock(output) as root:
        if root.exists():
            raise ValueError("run registration requires an unused output directory")
        with axis_run.open_inputs(args["training"], args["representations"]) as (data, vectors):
            config = bind(args, data, vectors)
        root.mkdir()
        record = root / "registration.json"
        axis_run.write_json(record, {**config, "arguments": args, "nonce": axis_run.PROCESS})
        return axis_run.digest(record)


def fit(root, native, model, data, config, registration_sha):
    """Fit on matched controls and retain the selected objective trace."""
    training = np.asarray(config["partitions"]["train"])
    valid = np.asarray(config["partitions"]["valid"])
    training_inputs = data.model_input("train")
    rng = np.random.default_rng(config["seed"])
    optimizer = native.optimizer(model)
    stopper = axis_run.EarlyStopping(patience=config["patience"], max_epochs=config["max_epochs"])
    history = []
    trace = ComponentTrace(model.loss)
    model.loss = trace
    for epoch in range(1, config["max_epochs"] + 1):
        loss = native.train_epoch(model, training_inputs, training, optimizer, rng)
        train_components = trace.per_condition("train")
        values, objective, mse = axis_run.validate(model, data, valid, config)
        valid_components = trace.per_condition("validation")
        if stopper.observe(objective, epoch):
            torch.save({"registration_sha256": registration_sha, "state_dict": model.state_dict(),
                        "epoch": epoch, "validation_loss": objective, "validation_mse": mse},
                       root / "checkpoint.pt")
            np.save(root / "validation_prediction.npy", values)
        history.append({"epoch": epoch, "train_loss": loss, "validation_loss": objective,
                        "validation_mse": mse, **train_components, **valid_components})
        if stopper.stopped:
            break
    pd.DataFrame(history).to_csv(root / "learning_curve.csv", index=False)
    return {"best_epoch": stopper.best_epoch, "validation_loss": stopper.best,
            "validation_mse": history[stopper.best_epoch - 1]["validation_mse"],
            "epochs": len(history), "parameters": sum(p.numel() for p in model.parameters())}


def execute(output, phase, paths):
    start = time.monotonic()
    with axis_run.publication_lock(output) as root:
        registration = root / "registration.json"
        config = axis_run.read_json(registration)
        registration_sha = axis_run.digest(registration)
        if set(paths) - set(INPUT_PATHS):
            raise ValueError("only verified input locations can be overridden")
        args = {**config["arguments"], **paths}
        receipt = PHASES[phase]
        sealed = None
        if (root / receipt).exists():
            if phase != "replay":
                raise ValueError(f"run phase already recorded: {phase}")
            sealed = axis_run.verified_receipt(root, receipt, registration_sha, SEALED_FILES)
        if phase != "train":
            training = axis_run.verified_receipt(root, PHASES["train"], registration_sha, axis_run.TRAINING_FILES)
        if phase == "replay":
            axis_run.verified_receipt(root, PHASES["predict"], registration_sha, axis_run.PREDICTION_FILES)
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        with axis_run.open_inputs(args["training"], args["representations"]) as (data, vectors):
            identity = {key: value for key, value in config.items() if key not in {"arguments", "nonce"}}
            if bind(args, data, vectors) != identity:
                raise ValueError("registered experiment, inputs, code or runtime changed")
            torch.use_deterministic_algorithms(True)
            axis_run.kpgt.configure_backend(args["device"])
            native = axis_run.kpgt.adapter(args["source"], config["model"]["transform"])
            model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
            objective = Objective(config["objective"]["form"], Fraction(config["objective"]["alpha"]))
            model.loss = partial(objective, model.loss)
            with TemporaryDirectory(prefix=".phase-", dir=root) as scratch_path:
                scratch = Path(scratch_path)
                if phase == "train":
                    details = fit(scratch, native, model, data, config, registration_sha)
                    files = axis_run.TRAINING_FILES
                else:
                    values, test = axis_run.predict_test(root, native, model, data, config, registration_sha, training)
                    details = axis_run.record_predictions(root, scratch, values, test, data, config, replay=phase == "replay")
                    files = axis_run.PREDICTION_FILES
                if sealed is not None:
                    return sealed
                if phase == "replay":
                    details.update(test_replay="exact", axis=AXIS, fold=args["fold"],
                                   smoke=args["smoke"], seed=config["seed"], schema=SCHEMA,
                                   protocol=PROTOCOL, objective=config["objective"])
                    hashes = axis_run.file_hashes(root, SEALED_FILES)
                else:
                    hashes = axis_run.file_hashes(scratch, files)
                    for name in files:
                        (scratch / name).replace(root / name)
                record = {
                    **details, "outputs": hashes, "registration_sha256": registration_sha,
                    "process_token": axis_run.PROCESS, "pid": os.getpid(),
                    "wall_seconds": round(time.monotonic() - start, 3),
                    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                }
                axis_run.write_json(scratch / receipt, record)
                (scratch / receipt).replace(root / receipt)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--loss", choices=FORMS, required=True)
    registration.add_argument("--alpha", choices=ALPHAS, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    registration.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    registration.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    registration.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    registration.add_argument("--smoke", action="store_true")
    for phase in PHASES:
        command = commands.add_parser(phase)
        for field in INPUT_PATHS:
            command.add_argument(f"--{field}", type=Path, help="same pinned input at a new location")
    for command in commands.choices.values():
        command.add_argument("--output", type=Path, required=True)
    args = vars(parser.parse_args())
    phase, output = args.pop("phase"), args.pop("output").resolve()
    options = {key: str(value.resolve()) if isinstance(value, Path) else value
               for key, value in args.items() if value is not None}
    if phase == "register":
        result = {"registration_sha256": register(output, {"axis": AXIS, **options})}
    else:
        result = execute(output, phase, options)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
