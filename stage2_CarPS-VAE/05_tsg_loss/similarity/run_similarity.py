"""Register, fit and exactly replay one bounded TranSiGen Delta similarity coordinate."""

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

from similarity_objective import ALPHAS, METRICS, SimilarityObjective
from run_reconstruction import bind as bind_reconstruction

from run_sweep import fit
from run_loss import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                      INPUT_PATHS, PHASES, SEALED_FILES, axis_run)

EXECUTION = Path(__file__).resolve().parent
SCHEMA = "stage2-tenv-transigen-bounded-delta-similarity-run-v1"
PROTOCOL = "matched-plate-control-bounded-delta-similarity-v1"


def bind(args, data, vectors):
    objective = SimilarityObjective(args["metric"], Fraction(args["alpha"]))
    config = bind_reconstruction({**args, "alpha": "1"}, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, objective=objective.identity())
    config["model"]["validation_objective"] = "registered bounded Delta similarity objective per condition"
    for name in ("run_similarity.py", "similarity_objective.py"):
        path = EXECUTION / name
        config["code_sha256"][path.relative_to(axis_run.REPO).as_posix()] = axis_run.digest(path)
    return axis_run.transigen.bind_experiment(args, config)


def register(output, args):
    return axis_run.transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = axis_run.kpgt.adapter(args["source"], config["model"]["transform"])
    model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
    objective = SimilarityObjective(args["metric"], Fraction(config["objective"]["alpha"]))
    model.loss = partial(objective, model.loss)
    return native, model


def execute(output, phase, paths, metric):
    return axis_run.transigen.execute_phase(
        output, phase, paths, bind=bind, build_model=build_model,
        fit=fit,
        expected_arguments={"metric": metric},
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
        command.add_argument("--metric", choices=METRICS, required=True)
    args = vars(parser.parse_args())
    phase, output = args.pop("phase"), args.pop("output").resolve()
    options = {key: str(value.resolve()) if isinstance(value, Path) else value
               for key, value in args.items() if value is not None}
    if phase == "register":
        result = {"registration_sha256": register(output, {"axis": AXIS, **options})}
    else:
        metric = options.pop("metric")
        result = execute(output, phase, options, metric)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
