"""Independently derive held-out Step 09 scores and publish plotted CSV inputs."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from axis_protocol import ARM, REFERENCE_AXIS, STAGE, coordinates
from axis_matrix import verify_matrix
from axis_runs import axis_run
from axis_scores import paired_scores, score_deltas
from chem_runs import axis_evaluate, read_run as read_chemical_run
from tenv_axes import open_tenv_axis_training
from tsg_axis.inputs import AxisInputs, score_fixed_prediction
from cpst_util.comparison import METRIC_NAMES, write_paired_metrics

SCHEMA = "cpst-unimol-axis-report-v2"
STEP04 = STAGE / "04_tsg_axis"
PLOT_TABLES = ("learning-curve.csv", "metrics.csv", "metrics-paired.csv")
REPORT_TABLES = ("bandwidths.csv", "condition-scores.csv", "epoch-objectives.csv",
                 "fold-scores.csv", "runs.csv", "summary.csv")


def calibrate(step03, step04):
    """Verify and reconstruct the accepted full-axis MMD bandwidths from Step 04."""
    acceptance = axis_run.read_json(step04 / "log/319-formal.json")["report"]
    path = step04 / "output/report-matched-control/bandwidths.csv"
    if axis_run.digest(path) != acceptance["output_sha256"]["bandwidths.csv"]:
        raise ValueError("original Step 04 full-axis MMD calibration is not sealed")
    sealed = axis_evaluate._csv(path).set_index("fold", verify_integrity=True)
    records = []
    with open_tenv_axis_training("hvg-5k", step03) as release:
        base = AxisInputs(release)
        genes = ";".join(map(str, release.genes.gene_index))
        for fold in range(3):
            calibration = axis_evaluate._bandwidth(base, fold)
            if calibration != {"fold": fold, **sealed.loc[fold].to_dict()}:
                raise ValueError(f"Step 04 full-axis MMD calibration differs for fold {fold}")
            rows = ";".join(map(str, base.rows(fold, "train")))
            records.append({**calibration, "gene_indices": genes, "train_rows": rows})
    return pd.DataFrame(records)


def score_run(key, run, data, base, bounded, bandwidth):
    axis, fold, _ = key
    path, config, replay, _, _ = run
    rows, _ = axis_evaluate._verify_inputs(path, config, data, base, axis, fold, bounded)
    conditions = axis_evaluate._csv(path / "conditions.csv")
    if (not conditions.equals(data.condition_table(rows).reset_index(drop=True))
            or len(conditions) != replay["test_conditions"]):
        raise ValueError(f"independent scoring found a changed held-out roster: {key}")
    prediction = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False)
    paired = score_fixed_prediction(data, fold, prediction, conditions.condition_id, rows=rows)
    for name, expected in (("predicted_delta", paired.predicted_delta), ("measured_delta", paired.measured_delta)):
        measured = np.load(path / f"{name}.npy", mmap_mode="r", allow_pickle=False)
        if measured.dtype != expected.dtype or not np.array_equal(measured, expected):
            raise ValueError(f"{key} {name} differs from the reconstructed matched-control Delta")
    observed = data.model_input("test")
    if not np.array_equal(paired.measured_delta, observed.expression[rows] - observed.control[rows]):
        raise ValueError(f"{key} measured Delta disagrees with independent treatment-minus-plate-control")
    genes = data.genes.gene_index.to_numpy()
    by_condition, by_fold = score_deltas(paired.measured_delta, paired.predicted_delta, genes, bandwidth)
    return by_condition.assign(condition_id=conditions.condition_id.to_numpy()), by_fold


def read_reused(step08, plan):
    """Resolve reuse decisions to their sealed historical files."""
    selected = {}
    for decision in plan["decisions"]:
        if decision["decision"] != "reuse":
            continue
        key = tuple(decision["coordinate"])
        axis, fold, seed = key
        if axis != REFERENCE_AXIS or key in selected:
            raise ValueError("invalid Step 08 reuse coordinate")
        path = step08 / "output/formal/matrix/runs" / ARM / f"fold-{fold}" / f"seed-{seed}"
        historical, run = read_chemical_run(path, False)
        if (historical != (ARM, fold, seed)
                or decision["registration_sha256"] != run[2]["registration_sha256"]
                or decision["manifest_sha256"] != axis_run.digest(path / "manifest.json")
                or decision["test_replay"] != "exact" or decision["validation_replay"] != "exact"):
            raise ValueError("sealed reuse receipt differs from its historical fit")
        selected[key] = run
    return selected


def report_tables(fresh, reused, step03, step04, bounded):
    roster = coordinates(bounded)
    if (set(fresh) & set(reused) or set(fresh) | set(reused) != set(roster)
            or any(key[0] != REFERENCE_AXIS for key in reused) or bounded and reused):
        raise ValueError("fit evidence must cover every registered coordinate exactly once")
    bandwidths = calibrate(step03, step04)
    lookups = bandwidths.set_index("fold").bandwidth.to_dict()
    scores, conditions, curves, runs = [], [], [], []
    with open_tenv_axis_training("hvg-5k", step03) as baseline:
        base = AxisInputs(baseline)
        for key in roster:
            axis, fold, seed = key
            origin = "Step 08" if key in reused else "Step 09"
            fit = (reused if key in reused else fresh)[key]
            with open_tenv_axis_training(axis, step03) as release:
                by_condition, by_fold = score_run(key, fit, AxisInputs(release), base, bounded, lookups[fold])
            path, config, replay, training, trace = fit
            identity = dict(axis=axis, fold=fold, seed=seed, origin=origin,
                            registration_sha256=replay["registration_sha256"])
            scores.append(by_fold.assign(**identity))
            conditions.append(by_condition.assign(**identity))
            curves.append(trace.assign(**identity))
            runs.append({
                **identity, "training_sha256": config["training_sha256"],
                "manifest_sha256": axis_run.digest(path / "manifest.json"),
                "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
                "prediction_sha256": replay["outputs"]["prediction.npy"],
                "trace_sha256": replay["outputs"]["learning_curve.csv"],
                "best_epoch": training["best_epoch"], "epochs": training["epochs"],
                "parameters": training["parameters"],
            })
    paired, summary = paired_scores(pd.concat(scores, ignore_index=True), roster)
    return {
        "fold-scores.csv": paired, "summary.csv": summary,
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "epoch-objectives.csv": pd.concat(curves, ignore_index=True),
        "runs.csv": pd.DataFrame(runs), "bandwidths.csv": bandwidths,
    }


def write_report(output, tables, bounded):
    output.mkdir()
    for name, table in tables.items():
        table.to_csv(output / name, index=False)
    receipt = {
        "schema": SCHEMA, "bounded": bounded, "scores": len(tables["fold-scores.csv"]),
        "outputs": axis_run.file_hashes(output, REPORT_TABLES),
    }
    axis_run.write_json(output / "manifest.json", receipt)
    return receipt


def figure_data(tables, output):
    """Publish the plotted scores and complete accepted epoch traces."""
    output.mkdir(parents=True)
    scores = tables["fold-scores.csv"]
    scores = scores.loc[scores.metric.isin(METRIC_NAMES), ["axis", "fold", "seed", "metric", "value"]]
    scores.rename(columns={"axis": "arm"}).to_csv(output / "metrics.csv", index=False)
    write_paired_metrics(output, REFERENCE_AXIS, mode="matched-coordinate")
    curve = tables["epoch-objectives.csv"][["axis", "fold", "seed", "epoch", "train_loss", "validation_loss"]]
    curve.rename(columns={"axis": "arm", "train_loss": "train", "validation_loss": "validation"}).to_csv(
        output / "learning-curve.csv", index=False)


def add_report_inputs(parser):
    parser.add_argument("--step04", type=Path, default=STEP04)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step03", type=Path, default=STAGE / "03_data_train")
    add_report_inputs(parser)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure-data", type=Path, required=True)
    parser.add_argument("--bounded", action="store_true")
    args = parser.parse_args()
    fresh = verify_matrix(args.matrix, args.bounded)
    tables = report_tables(fresh, {}, args.step03, args.step04, args.bounded)
    receipt = write_report(args.output, tables, args.bounded)
    figure_data(tables, args.figure_data)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
