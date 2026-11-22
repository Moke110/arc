"""Prepare and assess fixed transfer predictions on the accepted raw-first view."""
import argparse
import json
from pathlib import Path
import runpy
import sys
from tempfile import TemporaryDirectory

import numpy as np

STEP = Path(__file__).resolve().parents[1]
DOSE = STEP.parent / "06_tsg_dose"
for api in (STEP / "output/api", DOSE / "output/api"):
    sys.path.insert(0, str(api))

from tsg_raw_first import preparation
from tsg_raw_first.evaluation import MAX_ROWS, score
from tsg_raw_first.inputs import Inputs, SOURCES, read_csv, read_json, sha256
from tsg_transfer.tahoe import make_tahoe_splits
from tsg_transfer.tenv import EVALUATION_PATHS, FOLDS, freeze_split, load_split

PREP_SCHEMA = "tsg-transfer-raw-first-prepared-v1"
ASSESSMENT_SCHEMA = "tsg-transfer-raw-first-assessment-v1"
ASSESSMENT_FILES = ("identity.csv", "scores.csv")
SPLIT_FILES = ("tenv/splits.csv", "tenv/manifest.json", "tahoe/splits.csv")


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def stage_destination(path):
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return TemporaryDirectory(prefix=".transfer-", dir=path.parent)


def prepare(candidate, *, pair=None, view=None):
    """Bind tenv held test roles and independent Tahoe validation to one pair."""
    if pair is not None and pair.resolve() in candidate.resolve().parents:
        raise ValueError("transfer candidate cannot reside inside the input pair")
    with stage_destination(candidate) as temporary:
        staged = Path(temporary) / "candidate"
        staged.mkdir()
        if pair is None:
            pair = preparation.prepare(staged / "pair", view=view)
        tenv, tahoe = (Inputs(pair / source) for source in SOURCES)
        if tenv.pair_sha256 != tahoe.pair_sha256:
            raise ValueError("transfer inputs must share one frozen source pair")
        freeze_split(tenv, staged / "tenv")
        (staged / "tahoe").mkdir()
        make_tahoe_splits(tahoe.obs).to_csv(staged / "tahoe/splits.csv", index=False)
        write_json(staged / "manifest.json", {
            "schema": PREP_SCHEMA,
            "pair_sha256": tenv.pair_sha256,
            "view_generation": tenv.config["view_generation"],
            "sources": tenv.pair["sources"],
            "representations": {data.source: data.config["representation_sha256"] for data in (tenv, tahoe)},
            "training_seed": tenv.config["training_seed"],
            "split_seed": tenv.config["split_seed"],
            "outputs": {name: sha256(staged / name) for name in SPLIT_FILES},
        })
        staged.rename(candidate)


def evaluate(args):
    prepared = args.prepared
    receipt = read_json(prepared / "manifest.json")
    if (receipt.get("schema") != PREP_SCHEMA
            or receipt.get("outputs") != {name: sha256(prepared / name) for name in SPLIT_FILES}):
        raise ValueError("transfer preparation or split artifacts differ")
    pair = args.pair if args.pair is not None else prepared / "pair"
    data = Inputs(pair / "tenv", expected_pair_sha256=receipt["pair_sha256"])
    if (receipt.get("view_generation") != data.config["view_generation"]
            or receipt.get("sources") != data.pair["sources"]
            or set(data.pair["sources"]) != set(SOURCES)
            or sha256(pair / "tahoe_bulk/manifest.json") != data.pair["sources"]["tahoe_bulk"]):
        raise ValueError("transfer view or paired source provenance differs")
    rows = load_split(data, prepared / "tenv", args.fold)["test"]
    if not np.array_equal(rows, data.rows("cold-drug-dose", args.fold, "test")):
        raise ValueError("transfer test identities differ from the frozen tenv split")
    prediction = np.load(args.prediction, mmap_mode="r", allow_pickle=False)
    if prediction.dtype != np.float32 or prediction.ndim != 2:
        raise ValueError("expected finite float32 test predictions")
    if any(not np.isfinite(prediction[start:start + MAX_ROWS]).all()
           for start in range(0, len(prediction), MAX_ROWS)):
        raise ValueError("expected finite float32 test predictions")
    identities = read_csv(args.prediction_ids)
    if list(identities) != ["profile_id"]:
        raise ValueError("prediction identity CSV must contain only profile_id")
    scores, metrics = score(data, "cold-drug-dose", args.fold, "test", prediction,
                            identities.profile_id)
    with stage_destination(args.output) as temporary:
        staged = Path(temporary) / "assessment"
        staged.mkdir()
        data.obs.iloc[rows][["profile_id", "drug_parent_inchikey", "dose_um", "cell_line"]].to_csv(
            staged / "identity.csv", index=False)
        scores.to_csv(staged / "scores.csv", index=False)
        write_json(staged / "manifest.json", {
            "schema": ASSESSMENT_SCHEMA,
            "path": args.path,
            "fold": args.fold,
            "metrics": metrics,
            "preparation_sha256": sha256(prepared / "manifest.json"),
            "pair_sha256": data.pair_sha256,
            "view_generation": data.config["view_generation"],
            "prediction_sha256": sha256(args.prediction),
            "files": {name: sha256(staged / name) for name in ASSESSMENT_FILES},
        })
        staged.rename(args.output)


def compare(directories, output):
    """Make the existing study's paired and equal-fold plot tables."""
    receipts = []
    for directory in directories:
        receipt = read_json(directory / "manifest.json")
        if (receipt.get("schema") != ASSESSMENT_SCHEMA
                or receipt.get("files") != {name: sha256(directory / name) for name in ASSESSMENT_FILES}):
            raise ValueError(f"transfer assessment files differ: {directory}")
        receipts.append(receipt)
    for field in ("preparation_sha256", "pair_sha256", "view_generation"):
        if any(receipt[field] != receipts[0][field] for receipt in receipts):
            raise ValueError(f"transfer assessment {field} differs across folds and paths")
    # Reuse the original study-owned aggregation without editing its sealed report driver.
    historical = runpy.run_path(str(STEP / "execution/publish.py"))
    scores, differences = historical["comparison"](receipts)
    with stage_destination(output) as temporary:
        staged = Path(temporary) / "comparison"
        staged.mkdir()
        historical["write_csv"](staged / "transfer_scores.csv", scores)
        historical["write_csv"](staged / "transfer_differences.csv", differences)
        staged.rename(output)


def render(data, output):
    """Render the transfer study's two CSV panels using the shared palette."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from stage2_util import palette, plotting

    scores = plotting.read_rows(data / "transfer_scores.csv")
    differences = plotting.read_rows(data / "transfer_differences.csv")
    folds = [row["fold"] for row in differences]
    labels = [f"Fold {row['fold']}\n{float(row['held_dose_um']):g} µM"
              if row["fold"] != "mean" else "Equal-fold\nmean" for row in differences]
    x = np.arange(len(folds))
    fig, (absolute, contrasted) = plt.subplots(ncols=2, figsize=(10, 4.3), layout="constrained")
    colors = tuple(palette.color(family, 2) for family in ("gray", "orange", "green"))
    for offset, path, label, color in zip(
            (-0.23, 0, 0.23), EVALUATION_PATHS, ("tenv-only", "Tahoe-only", "Transfer"),
            colors, strict=True):
        values = {row["fold"]: float(row["pcc"]) for row in scores if row["path"] == path}
        absolute.bar(x + offset, [values[fold] for fold in folds], width=0.20, color=color, label=label)
    for offset, key, label, color in (
        (-0.12, "transfer_minus_tenv_only", "Transfer − tenv-only", colors[0]),
        (0.12, "transfer_minus_tahoe_only", "Transfer − Tahoe-only", colors[1]),
    ):
        contrasted.bar(x + offset, [float(row[key]) for row in differences],
                       width=0.20, color=color, label=label)
    for axis in (absolute, contrasted):
        axis.set_xticks(x, labels)
        axis.axhline(0, linewidth=0.8, color=palette.color("gray", 0))
        if "mean" in folds:
            axis.axvline(folds.index("mean") - 0.5, linewidth=0.7, color=palette.color("gray", 6))
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3,
                    frameon=False, fontsize=8)
    absolute.set(title="Held-out scores", ylabel="Sample-equal Delta PCC")
    contrasted.set(title="Paired differences", ylabel="Difference in Delta PCC")
    plotting.save_figure(fig, output, dpi=180, facecolor="white")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    inputs = prep.add_mutually_exclusive_group()
    inputs.add_argument("--view", type=Path)
    inputs.add_argument("--pair", type=Path)
    prep.add_argument("--candidate", type=Path, required=True)
    assessment = commands.add_parser("evaluate")
    assessment.add_argument("--pair", type=Path)
    assessment.add_argument("--fold", type=int, choices=FOLDS, required=True)
    assessment.add_argument("--path", choices=EVALUATION_PATHS, required=True)
    for name in ("prepared", "prediction", "prediction-ids", "output"):
        assessment.add_argument("--" + name, type=Path, required=True)
    contrast = commands.add_parser("compare")
    contrast.add_argument("assessments", type=Path, nargs=9)
    contrast.add_argument("--output", type=Path, required=True)
    plot = commands.add_parser("render")
    plot.add_argument("--data", type=Path, required=True)
    plot.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        prepare(args.candidate, pair=args.pair, view=args.view)
    elif args.command == "evaluate":
        evaluate(args)
    elif args.command == "compare":
        compare(args.assessments, args.output)
    else:
        render(args.data, args.output)


if __name__ == "__main__":
    main()
