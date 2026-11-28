"""Register, fit, predict or independently replay a Step 05 architecture and training seed."""

import argparse
import json
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
ROUND_FOUR = STAGE / "util/training/legacy_transigen/modules/measured_targets/execution"
sys.path[:0] = [str(STAGE), str(ROUND_FOUR)]

from measured_model import MeasuredAdapter
from run_loss import load_step04
from measured_roster import CONTROL, CONTROL_ONLY, DIRECT_DELTA, MEASURED_SHARED, SHAPES, SOURCE_SHAPES
from run_composition import MODEL_DEPENDENCIES as ROUND_THREE_DEPENDENCIES
from run_composition import architecture_config as round_three_config
from run_sweep import AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING, INPUT_PATHS, PHASES
from module_objective import configure_module_objective, objective_identity
from module_roster import fit_coordinates, protocol_roster
from util.training.transigen import fit_with_diagnostics

axis_run = load_step04("run")

SCHEMA = "cpst-tenv-transigen-module-run-v1"
PROTOCOL = "matched-plate-control-transigen-module-v1"
TRAINING_FILES = (*axis_run.TRAINING_FILES, "kl-diagnostics.csv")
PREDICTION_FILES = axis_run.PREDICTION_FILES
SEALED_FILES = (*TRAINING_FILES, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])
PHASE_FILES = {"train": TRAINING_FILES, "predict": PREDICTION_FILES, "replay": SEALED_FILES}

MODEL_DEPENDENCIES = (
    *ROUND_THREE_DEPENDENCIES, Path(__file__), ROUND_FOUR / "measured_model.py",
    ROUND_FOUR / "measured_roster.py", Path(__file__).with_name("module_objective.py"),
    Path(__file__).with_name("module_roster.py"),
)


def architecture_config(model, shape):
    source = SOURCE_SHAPES.get(shape, SOURCE_SHAPES[DIRECT_DELTA])
    fields = round_three_config(model, source)
    if shape == CONTROL_ONLY:
        updates = ("remove the measured treated encoder, Gaussian heads and decoder",)
        measured, readout = None, "independent signed Delta"
    elif shape in SOURCE_SHAPES:
        updates = ()
        measured = "treated expression"
        readout = "treated expression minus control" if shape == CONTROL else "independent signed Delta"
    else:
        updates = ("encode signed measured Delta in the second Gaussian and remove its decoder output ReLU",)
        measured, readout = "signed Delta", "independent signed Delta"
        if shape == MEASURED_SHARED:
            updates += ("reuse the signed measurement decoder for prediction",)
            readout = "shared signed Delta"
    fields["measured_encoding"], fields["delta_readout"] = measured, readout
    fields["constructor_order"].extend(updates)
    fields["initialization"] = "; ".join(fields["constructor_order"])
    fields["validation_objective"] = "full registered module objective, batch sums per condition"
    return fields


def bind(args, data, vectors):
    shape, fold, seed = (args[field] for field in ("shape", "fold", "seed"))
    if type(fold) is not int or type(seed) is not int or (shape, fold, seed) not in fit_coordinates(args["smoke"]):
        raise ValueError(f"unregistered new Step 05 fit: {(shape, fold, seed)}")
    config = axis_run.bind({**args, "smoke": False}, data, vectors)
    if args["smoke"]:
        config["max_epochs"] = 2
    config.update(schema=SCHEMA, protocol=PROTOCOL, shape=shape, seed=seed, objective=objective_identity(shape))
    config["model"].update(architecture_config(config["model"], shape))
    config["model"]["adaptations"][-1] = "shared cold-drug splits and fixed validation using the full module objective"
    config["code_sha256"].update({path.relative_to(axis_run.REPO).as_posix(): axis_run.digest(path)
                                  for path in MODEL_DEPENDENCIES})
    return axis_run.transigen.bind_experiment(args, config)


def register(output, args):
    return axis_run.transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = MeasuredAdapter(args["source"], config["model"]["transform"], config["shape"])
    model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
    configure_module_objective(model, config["shape"])
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
    folds, seeds = protocol_roster()
    registration.add_argument("--fold", choices=folds, type=int, required=True)
    registration.add_argument("--seed", choices=seeds, type=int, required=True)
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
