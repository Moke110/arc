"""Build the Step 04 comparison publication from its accepted matched-control evidence."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
sys.path[:0] = [str(STEP / "output/api"), str(STAGE / "util/output/api")]

from stage2_util.comparison import CURVE_FIELDS, FIGURES, METRIC_FIELDS, METRIC_NAMES, write_paired_metrics
from stage2_util.comparison_publication import render_receipt
from stage2_util.derivation.release import check_outputs, read_json, sha256, write_json
from stage2_util.metrics import signed_jaccard_scores
from tsg_axis.metrics import summarize_conditions

ARMS = ("hvg-3k", "hvg-4k", "hvg-5k", "denoise-3k", "denoise-4k", "denoise-5k")
REFERENCE = "denoise-5k"
ROSTER = {(arm, fold) for arm in ARMS for fold in range(3)}
EVIDENCE = ("conditions.csv", "genes.csv", "measured_delta.npy", "predicted_delta.npy", "learning_curve.csv")


def _coordinate(row):
    return row["axis"].replace("-deg-", "-"), row["fold"]


def _table(path):
    return pd.read_csv(path, float_precision="round_trip")


def _scored_run(run, fit, sealed, reference):
    arm, fold = _coordinate(fit)
    identity = arm, fold
    manifest_path = run / "manifest.json"
    manifest = read_json(manifest_path)
    registration = run / "registration.json"
    if (sha256(manifest_path) != fit["run_manifest_sha256"]
            or sha256(registration) != fit["registration_sha256"]
            or (manifest["axis"], manifest["fold"], manifest["seed"], manifest["smoke"])
            != (fit["axis"], fold, 42, False)):
        raise ValueError(f"run differs from accepted fit: {identity}")
    for name in EVIDENCE:
        if sha256(run / name) != manifest["outputs"][name]:
            raise ValueError(f"accepted fit payload differs: {identity} {name}")

    curve = _table(run / "learning_curve.csv")
    if curve.epoch.tolist() != list(range(1, fit["epochs"] + 1)):
        raise ValueError(f"accepted epoch trace is incomplete: {identity}")
    trace = [(arm, fold, 42, row.epoch, row.train_loss, row.validation_loss)
             for row in curve.itertuples(index=False)]

    roster = _table(run / "conditions.csv").condition_id.tolist()
    config = read_json(registration)
    if roster != config["conditions"]["test"] or roster != sealed.condition_id.tolist():
        raise ValueError(f"held-out condition identity or order differs: {identity}")
    genes = _table(run / "genes.csv").gene_index.to_numpy()
    measured = np.load(run / "measured_delta.npy", mmap_mode="r", allow_pickle=False)
    predicted = np.load(run / "predicted_delta.npy", mmap_mode="r", allow_pickle=False)
    if (measured.shape != predicted.shape or measured.shape != (fit["test_conditions"], len(genes))
            or not np.isfinite(measured).all() or not np.isfinite(predicted).all()):
        raise ValueError(f"invalid sealed matched-control Delta arrays: {identity}")
    jaccard, weighted = zip(*(signed_jaccard_scores(y, yhat, genes)
                              for y, yhat in zip(measured, predicted, strict=True)), strict=True)
    if not np.array_equal(jaccard, sealed.JAC_top200, equal_nan=True):
        raise ValueError(f"held-out JAC-top200 differs from the sealed conditions: {identity}")
    top200 = summarize_conditions(jaccard)
    top200_ref = reference.loc["JAC_top200"]
    if not np.isclose(top200["value"], top200_ref.value, rtol=0, atol=1e-15, equal_nan=True):
        raise ValueError(f"JAC-top200 aggregate differs from the sealed fold: {identity}")
    for name in ("n_conditions", "n_eligible", "n_undefined"):
        if top200[name] != top200_ref[name]:
            raise ValueError(f"JAC-top200 eligibility differs for {identity}: {name}")
    return trace, summarize_conditions(weighted)["value"]


def write_plot_data(report, matrix, fits, data):
    """Publish plotted values from verified fit files and sealed scientific report tables."""
    score_manifest = read_json(report / "manifest.json")
    report_runs = {_coordinate(row): row for row in score_manifest["runs"]}
    if ({_coordinate(row) for row in fits} != ROSTER or len(fits) != len(ROSTER)
            or set(report_runs) != ROSTER or len(score_manifest["runs"]) != len(ROSTER)):
        raise ValueError("comparison requires the complete 18-fit report and matrix")
    condition = _table(report / "condition-metrics.csv")
    folds = _table(report / "fold-metrics.csv")
    curves, metrics = [], []
    for fit in fits:
        arm, fold = _coordinate(fit)
        published = report_runs[arm, fold]
        if any(published[name] != fit[name] for name in ("registration_sha256", "run_manifest_sha256")):
            raise ValueError(f"sealed report uses another fit: {(arm, fold)}")
        run = matrix / "runs" / fit["axis"] / f"fold-{fold}"
        sealed = condition.loc[condition.axis.eq(fit["axis"]) & condition.fold.eq(fold)]
        reference = folds.loc[folds.axis.eq(fit["axis"]) & folds.fold.eq(fold)].set_index("metric")
        trace, weighted = _scored_run(run, fit, sealed, reference)
        curves.extend(trace)
        metrics.extend((arm, fold, 42, metric, weighted if metric == "wJAC-2000" else reference.loc[metric, "value"])
                       for metric in METRIC_NAMES)
    data.mkdir()
    pd.DataFrame(curves, columns=CURVE_FIELDS).to_csv(data / "learning-curve.csv", index=False)
    pd.DataFrame(metrics, columns=METRIC_FIELDS).to_csv(data / "metrics.csv", index=False)
    write_paired_metrics(data, REFERENCE)


def build_accepted(matrix, report, output):
    accepted = read_json(STEP / "log/319-formal.json")
    matrix_manifest = matrix / "manifest.json"
    report_manifest = report / "manifest.json"
    if (sha256(matrix_manifest) != accepted["matrix"]["manifest_sha256"]
            or sha256(report_manifest) != accepted["report"]["report_manifest_sha256"]):
        raise ValueError("matrix or report differs from the accepted Step 04 receipt")
    matrix_records = read_json(matrix_manifest)
    if matrix_records["fits"] != accepted["matrix"]["fits"]:
        raise ValueError("matrix coordinates differ from the accepted receipt")
    if check_outputs(report, read_json(report_manifest)) != accepted["report"]["output_sha256"]:
        raise ValueError("scientific report outputs differ from the accepted receipt")

    output.mkdir()
    figure = output / "figure"
    figure.mkdir()
    data = figure / "data"
    write_plot_data(report, matrix, matrix_records["fits"], data)
    renderer = STEP / "figure/script/render.py"
    receipt = {
        **render_receipt(STAGE, STEP.name, output, renderer),
        "evidence": {"matrix_manifest_sha256": sha256(matrix_manifest),
                     "report_manifest_sha256": sha256(report_manifest),
                     "build_sha256": sha256(Path(__file__)),
                     "shared_metrics_sha256": sha256(STAGE / "util/output/api/stage2_util/metrics.py"),
                     "sealed_jac_top200_parity": True,
                     "fits": len(matrix_records["fits"])},
    }
    write_json(output / "log/figure-render.json", receipt)
    print(json.dumps({"status": "PASS", "fits": receipt["evidence"]["fits"], "figures": list(FIGURES)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=STEP / "output/formal-matched-control")
    parser.add_argument("--report", type=Path, default=STEP / "output/report-matched-control")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build_accepted(args.matrix.resolve(), args.report.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
