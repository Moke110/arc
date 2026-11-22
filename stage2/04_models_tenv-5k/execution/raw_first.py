"""Run the active Stage 2 single-dose preparation, scoring and rendering."""
import argparse
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "output/api"))
from tenv5k_raw_first import workflow


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--candidate", type=Path, required=True)
    prepare.add_argument("--view", type=Path)
    prepare.add_argument("--representations", type=Path)
    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("--scenario", choices=("cold-drug", "warm-start"), required=True)
    evaluate.add_argument("--fold", type=int, choices=range(workflow.FOLDS), required=True)
    evaluate.add_argument("--model", choices=workflow.MODELS, required=True)
    evaluate.add_argument("--seed", type=int, choices=workflow.read_seeds(), required=True)
    evaluate.add_argument("--run-id", required=True)
    for field in ("inputs", "prediction", "prediction-ids", "output"):
        evaluate.add_argument(f"--{field}", type=Path, required=True)
    render = commands.add_parser("render")
    render.add_argument("--summary", type=Path, required=True)
    render.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        workflow.prepare(args.candidate, view_path=args.view, representations=args.representations)
    elif args.command == "evaluate":
        workflow.evaluate(args.inputs, args.prediction, args.prediction_ids, args.output,
                          scenario=args.scenario, fold=args.fold, model=args.model,
                          seed=args.seed, run_id=args.run_id)
    else:
        workflow.render(args.summary, args.output)


if __name__ == "__main__":
    main()
