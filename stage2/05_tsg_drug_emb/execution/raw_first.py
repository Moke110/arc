"""Prepare, evaluate and render the raw-first drug representation comparison."""
import argparse
from pathlib import Path
import sys

import numpy as np

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "output/api"))
from tsg_drug_emb_raw_first import ARMS, Comparison, evaluate, prepare, render
from tsg_raw_first.inputs import SOURCES, read_csv
from tsg_raw_first.preparation import prepare as prepare_inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    preparation = commands.add_parser("prepare")
    evaluation = commands.add_parser("evaluate")
    source_paths = preparation.add_mutually_exclusive_group()
    source_paths.add_argument("--inputs", type=Path)
    source_paths.add_argument("--view", type=Path)
    preparation.add_argument("--source", choices=SOURCES, default="tenv")
    evaluation.add_argument("--inputs", type=Path, required=True)
    for command in (preparation, evaluation):
        command.add_argument("--receipt", type=Path, required=True)
        command.add_argument("--representations", type=Path)
        command.add_argument("--representation-sha256")
    evaluation.add_argument("--arm", choices=(*ARMS, "Mean_cell"), required=True)
    evaluation.add_argument("--fold", type=int, choices=range(3), required=True)
    for argument in ("prediction", "prediction-ids", "output"):
        evaluation.add_argument("--" + argument, type=Path, required=True)
    plot = commands.add_parser("render")
    plot.add_argument("--summary", type=Path, required=True)
    plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.command == "render":
        render(args.summary, args.output)
        return

    options = {"representation": args.representations,
               "representation_sha256": args.representation_sha256}
    if args.command == "prepare":
        inputs = args.inputs
        if inputs is None:
            if args.receipt.exists() or args.receipt.is_symlink():
                raise FileExistsError(f"comparison receipt already exists: {args.receipt}")
            representations = ({args.source: (args.representations, args.representation_sha256)}
                               if args.representations is not None else None)
            if args.representation_sha256 is not None and representations is None:
                parser.error("--representation-sha256 requires --representations when preparing a freeze")
            freeze = prepare_inputs(
                args.receipt.with_suffix(".inputs"), view=args.view,
                representations=representations, sources=(args.source,),
            )
            inputs = freeze / args.source
        prepare(inputs, args.receipt, **options)
        return

    comparison = Comparison(args.receipt, args.inputs, **options)
    identities = read_csv(args.prediction_ids)
    if list(identities) != ["profile_id"]:
        raise ValueError("prediction identity CSV must contain only profile_id")
    predictions = np.load(args.prediction, mmap_mode="r")
    evaluate(comparison, args.arm, args.fold, predictions, identities.profile_id, args.output)


if __name__ == "__main__":
    main()
