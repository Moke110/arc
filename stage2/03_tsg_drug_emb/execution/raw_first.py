"""Prepare, evaluate and render the raw-first drug representation comparison."""
import argparse
from pathlib import Path
import sys

import numpy as np

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "output/api"))
from tsg_drug_emb_raw_first import ARMS, Comparison, evaluate, prepare, render
from tsg_raw_first.inputs import read_csv


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    preparation = commands.add_parser("prepare")
    evaluation = commands.add_parser("evaluate")
    for command in (preparation, evaluation):
        command.add_argument("--inputs", type=Path, required=True)
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
        prepare(args.inputs, args.receipt, **options)
        return

    comparison = Comparison(args.receipt, args.inputs, **options)
    identities = read_csv(args.prediction_ids)
    if list(identities) != ["profile_id"]:
        raise ValueError("prediction identity CSV must contain only profile_id")
    predictions = np.load(args.prediction, mmap_mode="r")
    evaluate(comparison, args.arm, args.fold, predictions, identities.profile_id, args.output)


if __name__ == "__main__":
    main()
