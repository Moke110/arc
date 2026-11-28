"""Register, fit, predict or independently replay a Step 06.2 prediction architecture."""

import argparse
import json
from pathlib import Path

from architecture_model import ArchitectureAdapter, axis_run, configure_objective, objective_identity
from architecture_roster import CONTROL, HEAD_PARAMETERS, NO_X2, SHAPES, objective_label
from shape_roster import dimensions
from run_sweep import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                       INPUT_PATHS, PHASES, bind as bind_loss)
from util.training.transigen import fit_with_diagnostics

SCHEMA = "stage2-tenv-transigen-prediction-architecture-run-v1"
PROTOCOL = "matched-plate-control-transigen-prediction-architecture-v1"
ALPHA = "1.75"
TRAINING_FILES = (*axis_run.TRAINING_FILES, "kl-diagnostics.csv")
PREDICTION_FILES = axis_run.PREDICTION_FILES
SEALED_FILES = (*TRAINING_FILES, *PREDICTION_FILES, PHASES["train"], PHASES["predict"])
PHASE_FILES = {"train": TRAINING_FILES, "predict": PREDICTION_FILES, "replay": SEALED_FILES}

ARCHITECTURE_CHANGES = {
    "output-delta": "separate signed Delta decoder with native hidden structure and selective Xavier",
    NO_X2: "repurpose the treated decoder as the signed Delta decoder",
    "z-interact_MLP": "replace prediction heads with default-initialized 532→145→200 ReLU MLP",
    "z-interact_FiLM": "replace prediction heads with identity-initialized FiLM 432→136→200 and 100→200 projection",
    "shared-z": "tie the treated encoder trunk to the control trunk",
    "z-res": "add the sampled control latent to the predicted mean",
    "shared-z-res": "tie encoder trunks and add the sampled control latent to the predicted mean",
}
ROUND_TWO_FILES = (Path(__file__), Path(__file__).with_name("architecture_model.py"),
                   Path(__file__).with_name("architecture_roster.py"))
MODEL_DEPENDENCIES = (
    *(Path(__file__).parents[2] / "expression_shape/execution" / name
      for name in ("shape_model.py", "shape_roster.py", "run_shape.py")),
    *ROUND_TWO_FILES,
)


def architecture_config(model, shape):
    fields = dimensions(CONTROL)
    fields["prediction_head_parameters"] = HEAD_PARAMETERS[shape]
    fields["constructor_order"] = model["constructor_order"].copy()
    if shape in ARCHITECTURE_CHANGES:
        fields["constructor_order"].append(ARCHITECTURE_CHANGES[shape])
    fields["initialization"] = "; ".join(fields["constructor_order"])
    fields["validation_objective"] = f"{objective_label(shape)}, batch sums per condition"
    return fields


def bind(args, data, vectors):
    shape = args["shape"]
    if shape not in SHAPES or (args["smoke"] and args["fold"] != 0) or (not args["smoke"] and shape == CONTROL):
        raise ValueError("Step 06.2 smoke uses fold 0; formal original fits must be reused from Step 06.1")
    config = bind_loss({**args, "alpha": ALPHA}, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, shape=shape, objective=objective_identity(shape))
    config["model"].update(architecture_config(config["model"], shape))
    config["code_sha256"].update({path.relative_to(axis_run.REPO).as_posix(): axis_run.digest(path)
                                  for path in MODEL_DEPENDENCIES})
    return axis_run.transigen.bind_experiment(args, config)


def register(output, args):
    return axis_run.transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = ArchitectureAdapter(args["source"], config["model"]["transform"], config["shape"])
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
