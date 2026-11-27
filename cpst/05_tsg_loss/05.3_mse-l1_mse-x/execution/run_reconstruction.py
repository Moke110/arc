"""Register, fit and exactly replay one Step 05.3 TranSiGen reconstruction MSE coordinate."""

import argparse
from fractions import Fraction
from functools import partial
import json
import os
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

import torch

from reconstruction_objective import PREDECESSOR_EXECUTION, ReconstructionObjective
from sweep_protocol import ALPHAS

from additive_objective import LOSS_COMPARISON
from run_sweep import fit
from run_loss import (AXIS, AXIS_STEP, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                      INPUT_PATHS, PHASES, SEALED_FILES, axis_run)

STEP = Path(__file__).resolve().parents[1]
SCHEMA = "stage2-tenv-transigen-reconstruction-mse-run-v1"
PROTOCOL = "matched-plate-control-reconstruction-mse-v1"


def bind(args, data, vectors):
    if args["axis"] != AXIS:
        raise ValueError("the reconstruction sweep requires the frozen denoise-5k axis")
    objective = ReconstructionObjective(Fraction(args["alpha"]))
    config = axis_run.bind(args, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, objective=objective.identity())
    config["model"]["validation_objective"] = "registered Step 05.3 objective, batch sums per condition"
    config["model"]["adaptations"][-1] = "shared cold-drug splits and fixed validation using the registered objective"
    dependencies = (
        Path(__file__), STEP / "execution/reconstruction_objective.py",
        STEP / "execution/sweep_protocol.py", PREDECESSOR_EXECUTION / "additive_objective.py",
        PREDECESSOR_EXECUTION / "run_sweep.py", LOSS_COMPARISON / "objective.py",
    )
    for path in dependencies:
        config["code_sha256"][path.relative_to(axis_run.REPO).as_posix()] = axis_run.digest(path)
    return axis_run.transigen.bind_experiment(args, config)


def register(output, args):
    return axis_run.transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = axis_run.kpgt.adapter(args["source"], config["model"]["transform"])
    model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
    objective = ReconstructionObjective(Fraction(config["objective"]["alpha"]))
    model.loss = partial(objective, model.loss)
    return native, model


def execute(output, phase, paths):
    return axis_run.transigen.execute_phase(
        output, phase, paths, bind=bind, build_model=build_model,
        fit=fit,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
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
