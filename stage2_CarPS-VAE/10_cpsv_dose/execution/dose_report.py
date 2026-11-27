"""Rescore the sealed Step 10 fits against the frozen matched-control release."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from dose_inputs import open_inputs
from dose_matrix import verify_matrix
from dose_protocol import DEFAULT_STEP03, REFERENCE, STAGE, TRANSIGEN, coordinates
from dose_runs import axis_run, input_identity
from chem_runs import axis_evaluate
from stage2_util.comparison import METRIC_NAMES
from tsg_axis.inputs import score_fixed_prediction

sys.path.insert(0, str(STAGE / "09_cpsv_axis/execution"))
from axis_scores import REPORT_METRICS, score_deltas

SCHEMA = "carps-vae-dose-report-v1"
PLOT_TABLES = ("learning-curve.csv", "metrics.csv", "metrics-paired.csv")
PRIMARY_METRICS = ("Δ PCC", "Δ RMSE")


def summarize_scores(scores, roster, *, metrics=REPORT_METRICS):
    """Compute within-fold–seed CarPS-VAE pairs and equal-coordinate summaries."""
    keys = ["arm", "fold", "seed", "metric"]
    expected = pd.MultiIndex.from_tuples([(*key, metric) for key in roster for metric in metrics], names=keys)
    indexed = scores.set_index(keys)
    if not indexed.index.is_unique or len(indexed) != len(expected) or not expected.difference(indexed.index).empty:
        raise ValueError("each dose coordinate must provide every registered metric exactly once")
    counts = indexed[["n_conditions", "n_eligible", "n_undefined"]]
    if (counts.isna().any().any() or (counts < 0).any().any()
            or not counts.eq(counts.round()).all().all()
            or not counts.n_conditions.eq(counts.n_eligible + counts.n_undefined).all()
            or not indexed.value.notna().eq(counts.n_eligible.gt(0)).all()
            or not np.isfinite(indexed.value.dropna()).all()):
        raise ValueError("score values must match the registered eligibility counts")
    paired = indexed.reindex(expected).reset_index()
    pairing = ["fold", "seed", "metric"]
    reference = paired.loc[paired.arm.eq(REFERENCE)].set_index(pairing)
    pairs = pd.MultiIndex.from_frame(paired[pairing])
    if not pairs.difference(reference.index).empty:
        raise ValueError("each fold and seed needs a matching no-dose reference")
    paired["reference_value"] = reference.value.reindex(pairs).to_numpy()
    paired.loc[paired.arm.eq(TRANSIGEN), "reference_value"] = np.nan
    paired["paired_difference"] = paired.value - paired.reference_value
    aggregates = dict(mean=("value", "mean"), sample_sd=("value", "std"), n_defined=("value", "count"),
                      paired_mean=("paired_difference", "mean"), paired_sd=("paired_difference", "std"),
                      n_paired=("paired_difference", "count"), n_coordinates=("value", "size"),
                      n_eligible=("n_eligible", "sum"), n_undefined=("n_undefined", "sum"))
    return (paired, paired.groupby(["arm", "fold", "metric"], sort=False, as_index=False).agg(**aggregates),
            paired.groupby(["arm", "metric"], sort=False, as_index=False).agg(**aggregates))


def score_run(key, run, data, bounded, bandwidth):
    """Independently rebind condition, control, gene and delta values for scoring."""
    _, fold, _ = key
    path, config, replay, _, _ = run
    positions = data.positions(fold, bounded)
    if any(config[name] != expected for name, expected in input_identity(data, fold, bounded).items()):
        raise ValueError(f"{key} registration differs from the frozen input identities")
    rows = positions["test"]
    conditions = axis_evaluate._csv(path / "conditions.csv")
    if (not conditions.equals(data.condition_table(rows).reset_index(drop=True))
            or not axis_evaluate._csv(path / "genes.csv").equals(data.genes)
            or len(conditions) != replay["test_conditions"]):
        raise ValueError(f"{key} prediction row or gene order differs")
    prediction = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False)
    target = score_fixed_prediction(data, fold, prediction, conditions.condition_id, rows=rows)
    for name, expected in (("predicted_delta", target.predicted_delta), ("measured_delta", target.measured_delta)):
        value = np.load(path / f"{name}.npy", mmap_mode="r", allow_pickle=False)
        if value.dtype != expected.dtype or not np.array_equal(value, expected):
            raise ValueError(f"{key} {name} differs from the rescored values")
    observed = data.model_input("test")
    if not np.array_equal(target.measured_delta, observed.expression[rows] - observed.control[rows]):
        raise ValueError(f"{key} frozen Delta disagrees with treatment minus the matched plate control")
    genes = data.genes.gene_index.to_numpy()
    by_condition, by_fold = score_deltas(target.measured_delta, target.predicted_delta, genes, bandwidth)
    return by_condition.assign(condition_id=conditions.condition_id.to_numpy()), by_fold


def report_tables(runs, step03, bounded):
    roster = coordinates(bounded)
    if set(runs) != set(roster):
        raise ValueError("the Step 10 report requires the complete registered roster")
    scores, conditions, curves, lineage = [], [], [], []
    with open_inputs(step03) as (data, _):
        gene_indices = ";".join(map(str, data.genes.gene_index))
        bandwidths = {
            fold: {**axis_evaluate._bandwidth(data, fold), "gene_indices": gene_indices,
                   "train_rows": ";".join(map(str, data.rows(fold, "train")))}
            for fold in sorted({fold for _, fold, _ in roster})
        }
        for key in roster:
            arm, fold, seed = key
            run = runs[key]
            by_condition, by_fold = score_run(key, run, data, bounded, bandwidths[fold]["bandwidth"])
            path, config, replay, trained, trace = run
            identity = dict(arm=arm, fold=fold, seed=seed, held_dose_um=config["population"]["held_dose_um"],
                            registration_sha256=replay["registration_sha256"])
            scores.append(by_fold.assign(**identity))
            conditions.append(by_condition.assign(**identity))
            curves.append(trace.assign(**identity))
            lineage.append({
                **identity, "training_sha256": config["training_sha256"], "axis_sha256": config["axis_sha256"],
                "representation": config["representation"]["name"],
                "representation_sha256": config["representation"]["vectors_sha256"],
                "manifest_sha256": axis_run.digest(path / "manifest.json"),
                "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
                "validation_prediction_sha256": replay["outputs"]["validation_prediction.npy"],
                "prediction_sha256": replay["outputs"]["prediction.npy"],
                "trace_sha256": replay["outputs"]["learning_curve.csv"],
                "best_epoch": trained["best_epoch"], "epochs": trained["epochs"],
                "parameters": trained["parameters"],
            })
    paired, foldwise, summary = summarize_scores(pd.concat(scores, ignore_index=True), roster)
    return {
        "fold-scores.csv": paired, "foldwise-summary.csv": foldwise, "summary.csv": summary,
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "epoch-objectives.csv": pd.concat(curves, ignore_index=True),
        "runs.csv": pd.DataFrame(lineage), "bandwidths.csv": pd.DataFrame(bandwidths.values()),
    }


def write_report(output, tables, bounded):
    output.mkdir()
    for name, table in tables.items():
        table.to_csv(output / name, index=False)
    receipt = {"schema": SCHEMA, "bounded": bounded, "primary_metrics": list(PRIMARY_METRICS),
               "scores": len(tables["fold-scores.csv"]), "outputs": axis_run.file_hashes(output, tables)}
    axis_run.write_json(output / "manifest.json", receipt)
    return receipt


def figure_data(tables, output):
    """Write the complete epoch, absolute or paired inputs for every panel."""
    output.mkdir(parents=True)
    keys = ["arm", "fold", "seed", "metric"]
    scores = tables["fold-scores.csv"]
    scores = scores.loc[scores.metric.isin(METRIC_NAMES)]
    scores[[*keys, "value"]].to_csv(output / "metrics.csv", index=False)
    paired = scores.loc[~scores.arm.isin((REFERENCE, TRANSIGEN))].rename(
        columns={"value": "arm_value", "paired_difference": "value"})
    paired.assign(reference=REFERENCE)[[*keys, "reference", "arm_value", "reference_value", "value"]].to_csv(
        output / "metrics-paired.csv", index=False)
    curve = tables["epoch-objectives.csv"][["arm", "fold", "seed", "epoch", "train_loss", "validation_loss"]]
    curve.rename(columns={"train_loss": "train", "validation_loss": "validation"}).to_csv(
        output / "learning-curve.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step03", type=Path, default=DEFAULT_STEP03)
    for name in ("matrix", "output", "figure-data"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--bounded", action="store_true")
    args = parser.parse_args()
    runs = verify_matrix(args.matrix, args.bounded)
    tables = report_tables(runs, args.step03, args.bounded)
    receipt = write_report(args.output, tables, args.bounded)
    figure_data(tables, args.figure_data)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
