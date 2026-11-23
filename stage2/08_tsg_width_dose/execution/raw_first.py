"""Prepare, evaluate and plot raw-first tenv width/dose comparisons."""

import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "execution"))
from run import METHODS, SCENARIO, SEED, WIDTHS, digest

sys.path.insert(0, str(STEP.parent / "07_tsg_dose/output/api"))
from tsg_raw_first import preparation
from tsg_raw_first.evaluation import score
from tsg_raw_first.inputs import KPGTInputs, read_csv
from stage2_util.seeds import read_seeds

sys.path.insert(0, str(STEP / "figure/raw_first"))
from render_width import render


def prepare(args):
    representation = (args.representations, args.representation_sha256) if args.representations else None
    if args.representation_sha256 and not representation:
        raise ValueError("representation digest requires an explicit representation directory")
    if args.tenv:
        preparation.prepare_kpgt(args.tenv, args.output, representation=representation)
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".width-input-", dir=args.output.parent) as temp:
        pair = Path(temp) / "pair"
        specifications = {"tenv": representation} if representation else None
        preparation.prepare(pair, view=args.view, representations=specifications, sources=("tenv",))
        preparation.prepare_kpgt(pair / "tenv", args.output, representation=representation)


def evaluate(args):
    data = KPGTInputs(args.inputs, expected_pair_sha256=args.pair_sha256)
    seed = data.config["training_seed"]
    if seed != SEED or seed not in read_seeds():
        raise ValueError("width study training seed differs from the accepted Stage 2 seed")
    prediction = np.load(args.prediction, mmap_mode="r")
    identities = read_csv(args.prediction_ids)
    if list(identities) != ["profile_id"]:
        raise ValueError("prediction identity CSV must contain only profile_id")
    samples, metrics = score(data, SCENARIO, args.fold, "test", prediction, identities.profile_id)
    if args.output.exists():
        raise FileExistsError(f"evaluation output already exists: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    generation = data.config["view_generation"]
    fields = {
        "source": data.source, "scenario": SCENARIO, "width": args.width,
        "method": args.method, "fold": args.fold, "seed": seed,
        "input_manifest_sha256": digest(args.inputs / "manifest.json"),
        "pair_sha256": data.pair_sha256, "view_sha256": generation["view_sha256"],
        "axis_sha256": generation["axis_sha256"], "corpus_sha256": generation["corpus_sha256"],
        "view_payload_sha256": data.config["view_payload_sha256"],
        "representation_sha256": data.config["representation_sha256"],
    }
    with TemporaryDirectory(prefix=".width-evaluate-", dir=args.output.parent) as temp:
        candidate = Path(temp) / "report"
        candidate.mkdir()
        samples.to_csv(candidate / "per_sample.csv", index=False)
        pd.DataFrame([{**fields, **metrics}]).to_csv(candidate / "summary.csv", index=False)
        candidate.rename(args.output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    source = prepare_parser.add_mutually_exclusive_group()
    source.add_argument("--tenv", type=Path, help="reuse a verified raw-first dose freeze")
    source.add_argument("--view", type=Path, help="use a specific raw-first view release")
    prepare_parser.add_argument("--representations", type=Path)
    prepare_parser.add_argument("--representation-sha256")
    prepare_parser.add_argument("--output", type=Path, required=True)

    evaluate_parser = commands.add_parser("evaluate")
    evaluate_parser.add_argument("--inputs", type=Path, required=True)
    evaluate_parser.add_argument("--pair-sha256")
    evaluate_parser.add_argument("--width", choices=WIDTHS, required=True)
    evaluate_parser.add_argument("--method", choices=METHODS, required=True)
    evaluate_parser.add_argument("--fold", type=int, choices=range(3), required=True)
    for field in ("prediction", "prediction-ids", "output"):
        evaluate_parser.add_argument(f"--{field}", type=Path, required=True)

    plot = commands.add_parser("render")
    plot.add_argument("--summary", type=Path, required=True)
    plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        prepare(args)
    elif args.command == "evaluate":
        evaluate(args)
    else:
        render(args.summary, args.output)


if __name__ == "__main__":
    main()
