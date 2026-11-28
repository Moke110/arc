"""Register, fit, predict or independently replay a Step 06.1 architecture."""

import argparse
from fractions import Fraction
from functools import partial
import json
from pathlib import Path

from shape_model import ShapeAdapter, axis_run
from shape_roster import SHAPES, SMOKE_SHAPES, dimensions
from run_sweep import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                       INPUT_PATHS, PHASES, AdditiveObjective, bind as bind_loss)
from util.training.transigen import fit_with_diagnostics, kl_diagnostics

SCHEMA = "stage2-tenv-transigen-shape-width-extension-run-v1"
PROTOCOL = "matched-plate-control-transigen-shape-width-extension-v1"
ALPHA = "1.75"
OBJECTIVE = AdditiveObjective(Fraction(ALPHA))
TRAINING_FILES = (*axis_run.TRAINING_FILES, "kl-diagnostics.csv")
PREDICTION_FILES = axis_run.PREDICTION_FILES
SEALED_FILES = (*TRAINING_FILES, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])
PHASE_FILES = {"train": TRAINING_FILES, "predict": PREDICTION_FILES, "replay": SEALED_FILES}


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
    return axis_run.transigen.bind_experiment(args, config)


def register(output, args):
    return axis_run.transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = ShapeAdapter(args["source"], config["model"]["transform"], config["shape"])
    model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
    model.loss = partial(OBJECTIVE, model.loss)
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
