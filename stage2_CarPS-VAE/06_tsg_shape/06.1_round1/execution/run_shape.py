"""Register, fit, predict or independently replay a Step 06.1 architecture."""

import argparse
from fractions import Fraction
from functools import partial
import json
import os
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

import pandas as pd
import torch

from shape_model import ShapeAdapter, axis_run
from shape_roster import SHAPES, SMOKE_SHAPES, dimensions
from run_sweep import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                       INPUT_PATHS, PHASES, AdditiveObjective, bind as bind_loss, fit as fit_loss)

SCHEMA = "stage2-tenv-transigen-shape-width-extension-run-v1"
PROTOCOL = "matched-plate-control-transigen-shape-width-extension-v1"
ALPHA = "1.75"
OBJECTIVE = AdditiveObjective(Fraction(ALPHA))
TRAINING_FILES = (*axis_run.TRAINING_FILES, "kl-diagnostics.csv")
PREDICTION_FILES = axis_run.PREDICTION_FILES
SEALED_FILES = (*TRAINING_FILES, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])
PHASE_FILES = {"train": TRAINING_FILES, "predict": PREDICTION_FILES, "replay": SEALED_FILES}


def kl_diagnostics(trace, latent):
    diagnostics = {"epoch": trace.epoch}
    for partition in ("train", "validation"):
        total = trace[f"{partition}_control_kl"] + trace[f"{partition}_treatment_kl"] + trace[f"{partition}_prediction_kl"]
        diagnostics[f"{partition}_total_kl"] = total
        diagnostics[f"{partition}_kl_per_latent_dimension"] = total / latent
    return pd.DataFrame(diagnostics)


def bind(args, data, vectors):
    if args["smoke"] and args["shape"] not in SMOKE_SHAPES:
        raise ValueError("shape is outside the bounded smoke roster")
    config = bind_loss({**args, "alpha": ALPHA}, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, shape=args["shape"])
    dims = dimensions(args["shape"])
    replacement = f"replace both predictive heads Linear{dims['predictive_head_dim']}→{dims['n_latent']} with PyTorch defaults"
    config["model"].update(dims, constructor_order=[config["model"]["constructor_order"][0], replacement])
    config["model"]["initialization"] = "; ".join(config["model"]["constructor_order"])
    for path in (Path(__file__), *(Path(__file__).with_name(name) for name in ("shape_model.py", "shape_roster.py"))):
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
            training = axis_run.verified_receipt(root, PHASES["train"], registration_sha, TRAINING_FILES)
        if phase == "replay":
            axis_run.verified_receipt(root, PHASES["predict"], registration_sha, PREDICTION_FILES)
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        with axis_run.open_inputs(args["training"], args["representations"]) as (data, vectors):
            identity = {key: value for key, value in config.items() if key not in {"arguments", "nonce"}}
            if bind(args, data, vectors) != identity:
                raise ValueError("registered experiment, inputs, code or runtime changed")
            torch.use_deterministic_algorithms(True)
            axis_run.kpgt.configure_backend(args["device"])
            native = ShapeAdapter(args["source"], config["model"]["transform"], config["shape"])
            model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
            model.loss = partial(OBJECTIVE, model.loss)
            with TemporaryDirectory(prefix=".phase-", dir=root) as scratch_path:
                scratch = Path(scratch_path)
                if phase == "train":
                    details = fit_loss(scratch, native, model, data, config, registration_sha)
                    trace = pd.read_csv(scratch / "learning_curve.csv", float_precision="round_trip")
                    kl_diagnostics(trace, config["model"]["n_latent"]).to_csv(scratch / "kl-diagnostics.csv", index=False)
                else:
                    values, test = axis_run.predict_test(root, native, model, data, config, registration_sha, training)
                    details = axis_run.record_predictions(root, scratch, values, test, data, config, replay=phase == "replay")
                if sealed is not None:
                    return sealed
                if phase == "replay":
                    details.update(test_replay="exact", axis=AXIS, fold=args["fold"],
                                   smoke=args["smoke"], seed=config["seed"], schema=SCHEMA, protocol=PROTOCOL,
                                   shape=config["shape"], objective=config["objective"])
                    hashes = axis_run.file_hashes(root, SEALED_FILES)
                else:
                    files = PHASE_FILES[phase]
                    hashes = axis_run.file_hashes(scratch, files)
                    for name in files:
                        (scratch / name).replace(root / name)
                record = {
                    **details, "outputs": hashes, "registration_sha256": registration_sha,
                    "process_token": axis_run.PROCESS, "pid": os.getpid(),
                    "wall_seconds": round(time.monotonic() - start, 3),
                    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    "cuda_peak_allocated_bytes": torch.cuda.max_memory_allocated() if args["device"] == "cuda" else None,
                }
                axis_run.write_json(scratch / receipt, record)
                (scratch / receipt).replace(root / receipt)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--shape", choices=SHAPES, required=True)
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
