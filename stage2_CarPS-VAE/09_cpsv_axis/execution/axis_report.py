"""Independently derive held-out Step 09 scores and publish plotted CSV inputs."""

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from axis_protocol import ARM, DISEASE_OVERLAP, DISEASE_SIZE, REFERENCE_AXIS, STAGE, coordinates
from axis_matrix import verify_matrix
from axis_runs import axis_run
from axis_scores import COVERAGE, paired_scores, score_deltas
from chem_runs import axis_evaluate, read_run as read_chemical_run
from tenv_axes import open_tenv_axis_training
from tsg_axis.disease_reference import disease_gene_sets
from tsg_axis.inputs import AxisInputs, score_fixed_prediction
from tsg_axis.metrics import fold_bandwidth

SCHEMA = "carps-vae-unimol-axis-report-v1"
DISEASES = {"nHF", "ICM", "genetic_cardiomyopathy", "TOX"}
DISEASE_RELEASE = STAGE / "02_gene_axis/02.1_hvg-deg/output/release"
STEP04 = STAGE / "04_tsg_axis"
CSV_FIGURES = {
    "axis-metrics": ("axis", False), "axis-metrics-paired": ("axis", True),
    "disease-metrics": ("disease", False), "disease-metrics-paired": ("disease", True),
    "disease-correct-coverage": ("disease", False), "disease-correct-coverage-paired": ("disease", True),
}


def disease_union(release):
    members, sources = disease_gene_sets(release)
    union = set(np.concatenate(tuple(members.values())))
    if set(members) != DISEASES or len(union) != DISEASE_SIZE:
        raise ValueError("disease coverage requires the frozen 3,481-gene eligible union")
    return union, {"memberships_sha256": sources, "gene_indices": list(map(int, sorted(union)))}


def calibrate(step03, step04, union):
    """Independently check the Step 04 full bandwidth and calibrate the separate disease bandwidth."""
    acceptance = axis_run.read_json(step04 / "log/319-formal.json")["report"]
    path = step04 / "output/report-matched-control/bandwidths.csv"
    if axis_run.digest(path) != acceptance["output_sha256"]["bandwidths.csv"]:
        raise ValueError("original Step 04 full-axis MMD calibration is not sealed")
    sealed = axis_evaluate._csv(path).set_index("fold", verify_integrity=True)
    records = []
    with open_tenv_axis_training("hvg-5k", step03) as release:
        base = AxisInputs(release)
        genes = release.genes.gene_index.to_numpy()
        keep = np.isin(genes, list(union))
        if np.count_nonzero(keep) != DISEASE_OVERLAP["hvg-5k"]:
            raise ValueError("HVG-DEG 5k disease overlap differs from the accepted design")
        for fold in range(3):
            full = axis_evaluate._bandwidth(base, fold)
            if full != {"fold": fold, **sealed.loc[fold].to_dict()}:
                raise ValueError(f"Step 04 full-axis MMD calibration differs for fold {fold}")
            rows = base.rows(fold, "train")
            parts = []
            for start in range(0, len(base.treatments), 64):
                stop = min(start + 64, len(base.treatments))
                indexes = rows[(rows >= start) & (rows < stop)]
                if len(indexes):
                    parts.append(release.read_treatments(start, stop, "delta")[indexes - start][:, keep])
            sliced = np.concatenate(parts)
            disease = {
                **full, "width": sliced.shape[1], "bandwidth": fold_bandwidth(sliced),
                "training_delta_values_sha256": sha256(sliced.tobytes()).hexdigest(),
            }
            for scope, result, indices in (("axis", full, genes), ("disease", disease, genes[keep])):
                records.append({
                    **result, "scope": scope,
                    "gene_indices": ";".join(map(str, indices)),
                    "train_rows": ";".join(map(str, rows)),
                })
    return pd.DataFrame(records)


def score_run(key, run, data, base, bounded, union, bandwidths):
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
    if np.count_nonzero(np.isin(genes, list(union))) != DISEASE_OVERLAP[axis]:
        raise ValueError(f"{key} disease union does not have the registered overlap")
    by_condition, by_fold = score_deltas(paired.measured_delta, paired.predicted_delta, genes,
                                          union, DISEASE_SIZE, bandwidths)
    by_condition["condition_id"] = conditions.condition_id.to_numpy()[by_condition.held_out_index]
    return by_condition.drop(columns="held_out_index"), by_fold


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


def report_tables(fresh, reused, step03, step04, disease_release, bounded):
    roster = coordinates(bounded)
    if (set(fresh) & set(reused) or set(fresh) | set(reused) != set(roster)
            or any(key[0] != REFERENCE_AXIS for key in reused) or bounded and reused):
        raise ValueError("fit evidence must cover every registered coordinate exactly once")
    union, disease = disease_union(disease_release)
    bandwidths = calibrate(step03, step04, union)
    lookups = {int(fold): scope.set_index("scope").bandwidth.to_dict()
               for fold, scope in bandwidths.groupby("fold")}
    scores, conditions, curves, runs = [], [], [], []
    with open_tenv_axis_training("hvg-5k", step03) as baseline:
        base = AxisInputs(baseline)
        for key in roster:
            axis, fold, seed = key
            origin = "Step 08" if key in reused else "Step 09"
            fit = (reused if key in reused else fresh)[key]
            with open_tenv_axis_training(axis, step03) as release:
                by_condition, by_fold = score_run(key, fit, AxisInputs(release), base, bounded, union, lookups[fold])
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
    }, disease


def write_report(output, tables, disease, bounded):
    output.mkdir()
    for name, table in tables.items():
        table.to_csv(output / name, index=False)
    axis_run.write_json(output / "disease-union.json", disease)
    receipt = {
        "schema": SCHEMA, "bounded": bounded, "scores": len(tables["fold-scores.csv"]),
        "outputs": axis_run.file_hashes(output, (*tables, "disease-union.json")),
    }
    axis_run.write_json(output / "manifest.json", receipt)
    return receipt


def figure_data(tables, output):
    """Write plotted replicates, aggregate means, sample SDs and zero references."""
    output.mkdir(parents=True)
    scores = tables["fold-scores.csv"]
    summary = tables["summary.csv"]
    keys = ["axis", "scope", "metric"]
    compound = scores.merge(summary.drop(columns=["overlap", "n_eligible", "n_undefined"]),
                            on=keys, validate="many_to_one")
    for name, (scope, paired) in CSV_FIGURES.items():
        coverage = "coverage" in name
        rows = compound[compound.scope.eq(scope) & compound.metric.eq(COVERAGE).eq(coverage)].copy()
        rows["plotted_value"] = rows.paired_difference if paired else rows.value
        rows["plotted_mean"] = rows.paired_mean if paired else rows["mean"]
        rows["plotted_sd"] = rows.paired_sd if paired else rows.sample_sd
        rows["zero_reference"] = np.where(paired | rows.metric.eq("Δ R²"), 0.0, np.nan)
        rows.to_csv(output / f"{name}.csv", index=False)
    curve = tables["epoch-objectives.csv"]
    curve[["axis", "fold", "seed", "epoch", "train_loss", "validation_loss"]].to_csv(
        output / "learning-curves.csv", index=False)


def add_report_inputs(parser):
    parser.add_argument("--step04", type=Path, default=STEP04)
    parser.add_argument("--disease-release", type=Path, default=DISEASE_RELEASE)


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
    tables, disease = report_tables(fresh, {}, args.step03, args.step04, args.disease_release, args.bounded)
    receipt = write_report(args.output, tables, disease, args.bounded)
    figure_data(tables, args.figure_data)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
