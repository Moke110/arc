"""Prepare, evaluate and render dose comparisons from the accepted raw-first view."""
import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "output/api"))
from tsg_raw_first import preparation
from tsg_raw_first.evaluation import score
from tsg_raw_first.inputs import Inputs, KPGTInputs, SOURCES, read_csv
from tsg_raw_first.reporting import render
from tsg_round1.protocol import SCENARIOS


def evaluate(args):
    reader = KPGTInputs if args.kpgt else Inputs
    data = reader(args.inputs, expected_pair_sha256=args.pair_sha256)
    predictions = np.load(args.prediction, mmap_mode="r")
    identities = read_csv(args.prediction_ids)
    if list(identities) != ["profile_id"]:
        raise ValueError("prediction identity CSV must contain only profile_id")
    samples, metrics = score(data, args.scenario, args.fold, args.partition,
                             predictions, identities.profile_id)
    if args.output.exists():
        raise FileExistsError(f"evaluation output already exists: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".evaluate-", dir=args.output.parent) as temp:
        candidate = Path(temp) / "report"
        candidate.mkdir()
        samples.to_csv(candidate / "per_sample.csv", index=False)
        generation = data.config["view_generation"]
        fields = {"source": data.source, "scenario": args.scenario, "fold": args.fold,
                  "partition": args.partition, "pair_sha256": data.pair_sha256,
                  "view_sha256": generation["view_sha256"], "axis_sha256": generation["axis_sha256"]}
        pd.DataFrame([{**fields, **metrics}]).to_csv(candidate / "summary.csv", index=False)
        candidate.rename(args.output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--view", type=Path)
    prepare.add_argument("--candidate", type=Path, required=True)
    prepare.add_argument("--source", choices=SOURCES, action="append")
    for source in SOURCES:
        prepare.add_argument(f"--representations-{source.replace('_', '-')}", type=Path)
    kpgt = commands.add_parser("prepare-kpgt")
    kpgt.add_argument("--tenv", type=Path, required=True)
    kpgt.add_argument("--candidate", type=Path, required=True)
    kpgt.add_argument("--representations", type=Path)
    evaluate_parser = commands.add_parser("evaluate")
    evaluate_parser.add_argument("--inputs", type=Path, required=True)
    evaluate_parser.add_argument("--pair-sha256")
    evaluate_parser.add_argument("--kpgt", action="store_true")
    evaluate_parser.add_argument("--scenario", choices=SCENARIOS, required=True)
    evaluate_parser.add_argument("--fold", type=int, choices=range(3), required=True)
    evaluate_parser.add_argument("--partition", choices=("validation", "test"), default="test")
    for argument in ("prediction", "prediction-ids", "output"):
        evaluate_parser.add_argument("--" + argument, type=Path, required=True)
    plot = commands.add_parser("render")
    plot.add_argument("--summary", type=Path, required=True)
    plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        representations = {source: (path, None) for source in SOURCES
                           if (path := getattr(args, f"representations_{source}")) is not None}
        preparation.prepare(args.candidate, view=args.view, representations=representations,
                            sources=args.source or SOURCES)
    elif args.command == "prepare-kpgt":
        specification = (args.representations, None) if args.representations else None
        preparation.prepare_kpgt(args.tenv, args.candidate, representation=specification)
    elif args.command == "evaluate":
        evaluate(args)
    else:
        render(args.summary, args.output)


if __name__ == "__main__":
    main()
