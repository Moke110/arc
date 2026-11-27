"""Register, fit, predict or independently replay a Step 06.3 architecture and training seed."""

import argparse
import json
import os
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

import pandas as pd
import torch

from composition_model import CompositionAdapter, axis_run, configure_objective, objective_identity
from composition_roster import (COMBINATIONS, HEAD_PARAMETERS, NORMALIZED_FILM_SHAPES, ROUND_TWO_SHAPES, SEEDS, SHAPES,
                                fit_coordinates, objective_label)
from run_architecture import ARCHITECTURE_CHANGES, MODEL_DEPENDENCIES as ROUND_TWO_DEPENDENCIES
from run_architecture import architecture_config as round_two_config
from run_shape import fit_with_diagnostics
from run_sweep import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                       INPUT_PATHS, PHASES, bind as bind_loss, fit as fit_loss)

SCHEMA = "stage2-tenv-transigen-architecture-composition-run-v2"
PROTOCOL = "matched-plate-control-transigen-architecture-composition-v2"
ALPHA = "1.75"
TRAINING_FILES = (*axis_run.TRAINING_FILES, "kl-diagnostics.csv")
PREDICTION_FILES = axis_run.PREDICTION_FILES
SEALED_FILES = (*TRAINING_FILES, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])
PHASE_FILES = {"train": TRAINING_FILES, "predict": PREDICTION_FILES, "replay": SEALED_FILES}

MODEL_DEPENDENCIES = (
    *ROUND_TWO_DEPENDENCIES, Path(__file__),
    Path(__file__).with_name("composition_model.py"),
    Path(__file__).with_name("composition_roster.py"),
)


def architecture_config(model, shape):
    fields = round_two_config(model, ROUND_TWO_SHAPES[shape])
    fields["prediction_head_parameters"] = HEAD_PARAMETERS[shape]
    if shape in COMBINATIONS:
        fields["constructor_order"].append(ARCHITECTURE_CHANGES["z-interact_FiLM"])
    if shape in NORMALIZED_FILM_SHAPES:
        fields["prediction_head_normalization"] = {
            "type": "LayerNorm", "normalized_shape": [100], "elementwise_affine": False,
            "eps": 1e-5, "position": "after_film_before_projection",
        }
        fields["constructor_order"].append("apply parameter-free LayerNorm to the modulated control latent")
    fields["initialization"] = "; ".join(fields["constructor_order"])
    fields["validation_objective"] = f"{objective_label(shape)}, batch sums per condition"
    return fields


def bind(args, data, vectors):
    shape, fold, seed = (args[field] for field in ("shape", "fold", "seed"))
    if type(fold) is not int or type(seed) is not int or (shape, fold, seed) not in fit_coordinates(args["smoke"]):
        raise ValueError(f"unregistered new Step 06.3 fit: {(shape, fold, seed)}")
    config = bind_loss({**args, "alpha": ALPHA}, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, shape=shape, seed=seed, objective=objective_identity(shape))
    config["model"].update(architecture_config(config["model"], shape))
    config["code_sha256"].update({path.relative_to(axis_run.REPO).as_posix(): axis_run.digest(path)
                                  for path in MODEL_DEPENDENCIES})
    return axis_run.transigen.bind_experiment(args, config)


def register(output, args):
    return axis_run.transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = CompositionAdapter(args["source"], config["model"]["transform"], config["shape"])
    model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
    configure_objective(model, config["shape"])
    return native, model


def execute(output, phase, paths):
    return axis_run.transigen.execute_phase(
        output, phase, paths, bind=bind, build_model=build_model,
        fit=fit_with_diagnostics,
        training_files=TRAINING_FILES, record_cuda=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--shape", choices=SHAPES, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--seed", choices=SEEDS, type=int, required=True)
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
