"""Build an isolated report candidate from the replay-accepted formal matrix."""
import argparse
import json
from pathlib import Path
from runpy import run_path
from statistics import mean, stdev
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
DOSE = STEP.parent / "04_tsg_dose"
sys.path[:0] = [str(STEP / "output/api"), str(DOSE / "output/api")]
from tsg_drug_emb.representations import ARMS, RepresentationInputs, sha256
from tsg_round1.inputs import Inputs
from tsg_round1.protocol import METRICS, summary

BASELINE = "Mean_cell"
METHODS = (BASELINE, *ARMS)
IDENTITY = ["fold", "profile_id", "drug_parent", "dose_um", "cell_line"]
SCORES = [*IDENTITY[1:], *METRICS, "constant_truth", "constant_prediction", "nonfinite"]
FOLDS = range(3)
ORIGIN = dict(source="tenv", scenario="cold-drug", seed=42)
PHASES = ("fit", "validation-replay", "test", "test-replay")
FIGURES = ("metrics.png", "learning_curves.png", "metrics_cold_drug_dose.png")
RENDERER = STEP / "figure/script/render.py"


def parent_pcc(scores):
    result = scores.assign(pcc=scores.pcc.where(np.isfinite(scores.pcc)))
    return result.groupby(["fold", "drug_parent"], as_index=False).pcc.agg(
        pcc="mean", defined_profiles="count", profiles="size")


def paired_pcc(arm, reference):
    pairs = arm[[*IDENTITY, "pcc"]].merge(
        reference[[*IDENTITY, "pcc"]], on=IDENTITY, how="outer",
        validate="one_to_one", indicator=True, suffixes=("", "_reference"))
    if not pairs._merge.eq("both").all():
        raise ValueError("paired scores differ in held-out profile identity")
    finite = np.isfinite(pairs.pcc) & np.isfinite(pairs.pcc_reference)
    pairs["delta_pcc"] = (pairs.pcc - pairs.pcc_reference).where(finite)
    return pairs.groupby(["fold", "drug_parent"], as_index=False).delta_pcc.agg(
        delta_pcc="mean", defined_pairs="count", profiles="size")


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def read_scores(path, data, fold):
    rows = data.rows("cold-drug", fold, "test")
    obs = data.obs.iloc[rows]
    expected = pd.DataFrame({
        "profile_id": obs.index.astype(str), "drug_parent": obs.drug_parent_inchikey.to_numpy(),
        "dose_um": obs.dose_um.to_numpy(), "cell_line": obs.cell_line.astype(str).to_numpy(),
    })
    scores = pd.read_csv(path, float_precision="round_trip",
                         dtype=dict(profile_id=str, drug_parent=str, cell_line=str))
    if (set(scores) != set(SCORES) or scores.profile_id.duplicated().any()
            or not scores[IDENTITY[1:]].equals(expected)):
        raise ValueError(f"score identities differ from frozen test fold {fold}: {path}")
    return scores.assign(fold=fold)


def read_matrix(matrix, data, representations, accepted):
    from fit import verify_phase

    input_digest = sha256(data.directory / "manifest.json")
    feature_digest = sha256(representations / "manifest.json")
    expected_fits = {(fit["fold"], fit["arm"]): fit for fit in accepted["fits"]}
    expected_baseline = {fit["fold"]: fit for fit in accepted["mean_cell"]}
    if (len(accepted["fits"]) != 21
            or set(expected_fits) != {(fold, arm) for fold in FOLDS for arm in ARMS}
            or len(accepted["mean_cell"]) != 3 or set(expected_baseline) != set(FOLDS)):
        raise ValueError("formal acceptance does not bind precisely 21 fits and 3 baseline rows")
    scores = {arm: [] for arm in METHODS}
    curves, runs = [], []
    for fold in FOLDS:
        root = matrix / "runs" / f"fold-{fold}"
        if {path.name for path in root.iterdir()} != set(ARMS):
            raise ValueError(f"formal arm inventory mismatch in fold {fold}")
        selection = sha256(matrix / f"selection-fold-{fold}.json")
        baseline = matrix / "mean-cell" / f"fold-{fold}"
        baseline_receipt = verify_phase(baseline, "test")
        if (sha256(baseline / "manifest.json") != expected_baseline[fold]["manifest_sha256"]
                or baseline_receipt["arm"] != BASELINE or baseline_receipt["fold"] != fold
                or baseline_receipt["inputs_sha256"] != input_digest):
            raise ValueError(f"Mean_cell output differs from accepted fold {fold}")
        scores[BASELINE].append(read_scores(baseline / "test_per_sample.csv", data, fold))
        for arm in ARMS:
            run = root / arm
            receipts = {phase: verify_phase(run / phase, phase) for phase in PHASES}
            digests = {phase: sha256(run / phase / "manifest.json") for phase in PHASES}
            accepted_run = expected_fits[fold, arm]
            bindings = {
                "fit_sha256": digests["fit"], "replay_sha256": digests["validation-replay"],
                "test_sha256": digests["test"], "test_replay_sha256": digests["test-replay"],
                "selection_sha256": selection,
            }
            if any(accepted_run[key] != digest for key, digest in bindings.items()):
                raise ValueError(f"run differs from formal acceptance: {run}")
            fit = receipts["fit"]
            config = fit["config"]
            expected = dict(arm=arm, fold=fold, seed=42, scenario="cold-drug",
                            inputs_sha256=input_digest, representations_sha256=feature_digest,
                            smoke_rows=None, max_epochs=1000, patience=50)
            if (any(config[key] != value for key, value in expected.items())
                    or any(config["order"][part] != data.rows("cold-drug", fold, part).tolist()
                           for part in ("train", "validation", "test"))
                    or len({receipt["process_token"] for receipt in receipts.values()}) != len(PHASES)):
                raise ValueError(f"formal fit registration mismatch: {run}")
            val, test, replay = (receipts[phase] for phase in PHASES[1:])
            if (any(receipt["fit_sha256"] != digests["fit"] for receipt in (val, test, replay))
                    or val["prediction_sha256"] != fit["outputs"]["validation_prediction.npy"]
                    or any(receipt["validation_replay_sha256"] != digests["validation-replay"]
                           for receipt in (test, replay))
                    or replay["test_sha256"] != digests["test"]
                    or replay["prediction_sha256"] != test["outputs"]["test_prediction.npy"]):
                raise ValueError(f"independent phase replay bindings differ: {run}")
            scores[arm].append(read_scores(run / "test/test_per_sample.csv", data, fold))
            learning = pd.read_csv(run / "fit/learning_curve.csv", float_precision="round_trip")
            best = learning.loc[learning.validation_mse.idxmin()]
            if (learning.epoch.tolist() != list(range(1, fit["epochs"] + 1))
                    or not np.isfinite(learning.select_dtypes("number")).all().all()
                    or fit["best_epoch"] != int(best.epoch)
                    or fit["validation_mse"] != float(best.validation_mse)
                    or fit["optimizer_steps"] != int(learning.optimizer_steps.iloc[-1])):
                raise ValueError(f"learning history does not match restored selection: {run}")
            curves.append(learning.assign(**ORIGIN, fold=fold, method=arm))
            runs.append(dict(fold=fold, arm=arm, width=config["width"], seed=42,
                             fit_epochs=fit["epochs"], optimizer_steps=fit["optimizer_steps"],
                             best_epoch=fit["best_epoch"], elapsed_seconds=fit["elapsed_seconds"],
                             validation_mse=fit["validation_mse"], phase_sha256=digests,
                             checkpoint_sha256=fit["outputs"]["checkpoint.pt"],
                             test_scores_sha256=test["outputs"]["test_per_sample.csv"],
                             selection_sha256=selection, resources=accepted_run["resources"]))
    pooled = {arm: pd.concat(parts, ignore_index=True) for arm, parts in scores.items()}
    for arm, frame in pooled.items():
        if (frame.profile_id.duplicated().any() or set(frame.profile_id) != set(data.obs.index.astype(str))
                or frame.drug_parent.nunique() != len(data.parents)
                or frame.groupby("drug_parent").fold.nunique().max() != 1
                or frame.nonfinite.any()):
            raise ValueError(f"outer replay coverage differs from the complete tenv roster: {arm}")
    return scores, pooled, pd.concat(curves, ignore_index=True), runs


def report_tables(scores, pooled, curves):
    dots, spreads, doses, parents, primary, pairs, effects = [], [], [], [], [], [], []
    for arm in METHODS:
        group = dict(**ORIGIN, method=arm)
        folds = [summary(scores[arm][fold]) for fold in FOLDS]
        dots += [dict(**group, fold=fold, **record) for fold, record in enumerate(folds)]
        for metric in METRICS:
            values = [row[metric] for row in folds if row[metric] is not None]
            spreads.append(dict(**group, metric=metric, n_defined_folds=len(values),
                                mean=mean(values) if values else None,
                                sd=stdev(values) if len(values) > 1 else None,
                                estimand="equal-fold mean and sample SD (ddof=1)"))
        for dose, result in pooled[arm].groupby("dose_um"):
            doses.append(dict(**group, dose_um=dose, **summary(result)))
        parent = parent_pcc(pooled[arm]).assign(**group)
        parents.append(parent)
        primary.append(dict(**group, pcc=parent.pcc.mean(), defined_parents=parent.pcc.count(),
                            parents=len(parent), defined_profiles=parent.defined_profiles.sum(),
                            profiles=parent.profiles.sum()))
        matched = paired_pcc(pooled[arm], pooled["ECFP4C"]).assign(**group)
        pairs.append(matched)
        for fold, result in [*matched.groupby("fold"), ("all", matched)]:
            effects.append(dict(**group, fold=fold, delta_pcc=result.delta_pcc.mean(),
                                defined_parents=result.delta_pcc.count(), parents=len(result),
                                defined_pairs=result.defined_pairs.sum(), profiles=result.profiles.sum()))
    return {"runs": pd.DataFrame(dots), "summary": pd.DataFrame(spreads),
            "dose": pd.DataFrame(doses), "curves": curves,
            "parent_pcc": pd.concat(parents, ignore_index=True),
            "pooled": pd.DataFrame(primary), "paired_parents": pd.concat(pairs, ignore_index=True),
            "paired": pd.DataFrame(effects)}


def report_readme(tables):
    primary = tables["pooled"].set_index("method")
    paired = tables["paired"].loc[lambda frame: frame.fold.eq("all")].set_index("method")
    lines = [
        "# TranSiGen drug representation comparison", "",
        "Seven representations were compared on all 459 tenv parents across three cold-drug",
        "folds at seed 42. Each arm used FiLM, the same G_5000 input, and the same fold rows.",
        "The first drug-encoder layer accepted the native feature width for each arm and",
        "produced 400 outputs. Later layers shared architecture and initialization. Each",
        "response model trained from scratch. Inner-validation MSE selected checkpoints.", "",
        "| Representation | Parent-equal OOF Δ PCC | Paired Δ PCC vs ECFP4C | Defined parents |",
        "|---|---:|---:|---:|",
    ]
    lines += [f"| {arm} | {primary.loc[arm, 'pcc']:.5f} | {paired.loc[arm, 'delta_pcc']:+.5f} | "
              f"{primary.loc[arm, 'defined_parents']}/{primary.loc[arm, 'parents']} |"
              for arm in METHODS]
    lines += [
        "", "The primary result averages finite per-profile Δ PCC within each held-out parent",
        "and gives each parent equal weight. The paired result uses profiles with defined",
        "scores for both arms before averaging within parent. Fold SD describes dispersion",
        "across the three folds for this single seed.", "",
        "- [Per-fold metrics](figure/data/runs.csv) and [equal-fold summary](figure/metrics.png)",
        "- [Native training loss and validation MSE](figure/learning_curves.png)",
        "- [Test metrics by measured dose](figure/metrics_cold_drug_dose.png)",
        "- [Primary pooled result](figure/data/pooled.csv) and [paired effects](figure/data/paired.csv)",
        "- [Run inventory](log/run_inventory.json) and [acceptance](log/acceptance.json)", "",
        "Run `python figure/script/render.py --output <new-directory>` from this step to",
        "reproduce the plots using the checked-in CSV files.", "",
    ]
    return "\n".join(lines)


def prepare(args):
    # Trust the formal gate before accepting its scientifically complete phase receipts.
    acceptance = args.matrix / "acceptance.json"
    accepted = json.loads(acceptance.read_text())
    if (accepted.get("schema") != "tsg-drug-emb-formal-acceptance-v1"
            or accepted.get("status") != "PASS" or accepted.get("formal_fits") != 21
            or accepted["plan_sha256"] != sha256(args.matrix / "formal-plan.json")):
        raise ValueError("formal fit matrix is not replay-accepted")
    inputs = Inputs(args.inputs)
    RepresentationInputs(inputs, args.representations)
    scores, pooled, curves, inventory = read_matrix(
        args.matrix, inputs, args.representations, accepted)
    tables = report_tables(scores, pooled, curves)
    target = args.destination
    data = target / "figure/data"
    data.mkdir(parents=True)
    (target / "log").mkdir()
    for name, frame in tables.items():
        frame.to_csv(data / f"{name}.csv", index=False)
    renderer = run_path(str(RENDERER))
    renderer["prepare_axes"](data)
    for output in (target / "figure", None):
        with TemporaryDirectory() as replay:
            subprocess.run([sys.executable, str(RENDERER), "--data", str(data),
                            "--output", str(output or replay)], check=True)
            if output is None:
                for name in FIGURES:
                    if sha256(Path(replay) / name) != sha256(target / "figure" / name):
                        raise ValueError(f"fresh-process renderer differs: {name}")
    (target / "README.md").write_text(report_readme(tables))
    write_json(target / "log/run_inventory.json", dict(fits=inventory, mean_cell=accepted["mean_cell"]))
    files = {str(path.relative_to(target)): sha256(path)
             for path in sorted(target.rglob("*")) if path.is_file()}
    write_json(target / "log/acceptance.json", dict(
        schema="tsg-drug-emb-report-candidate-v1", status="candidate", files=files,
        formal_acceptance_sha256=sha256(acceptance),
        inputs_sha256=sha256(args.inputs / "manifest.json"),
        representations_sha256=sha256(args.representations / "manifest.json"),
        renderer_sha256=sha256(RENDERER), fresh_csv_only_replay="PASS"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("matrix", "inputs", "representations", "destination"):
        parser.add_argument(f"--{field}", required=True, type=Path)
    prepare(parser.parse_args())
