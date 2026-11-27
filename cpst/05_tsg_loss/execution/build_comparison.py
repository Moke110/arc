"""Build current Step 05 comparison figures from the accepted fit and report releases."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

LOSS = Path(__file__).resolve().parents[1]
STAGE = LOSS.parent
sys.path[:0] = [str(LOSS / "figure/script"), *(str(STAGE / path) for path in (
    "util/output/api", "04_tsg_axis/output/api", "03_data_train/output/api",
))]

from loss_comparison import RENDERER, STUDIES
from cpst_util.comparison import CURVE_FIELDS, METRIC_FIELDS, METRIC_NAMES, write_paired_metrics
from cpst_util.comparison_publication import render_receipt
from cpst_util.derivation.release import check_outputs, read_json, sha256, write_json
from training_data import open_tenv_training, tenv_axis_output
from tsg_axis.comparison_evidence import read_fixed_jaccard
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import summarize_conditions

TRACE = ("epoch", "train_loss", "validation_loss")
OLD_METRICS = tuple(metric for metric in METRIC_NAMES if metric != "wJAC-2000")


def _table(path):
    return pd.read_csv(path, dtype={"alpha": str}, float_precision="round_trip")


def _arm(row):
    alpha = row["alpha"]
    return f"{row['loss']}:{alpha}" if "loss" in row else alpha


def _accepted_sources(source, receipt, study):
    matrix = source / "output/formal/matrix"
    report = source / "output/formal" / (f"reports/{study.similarity}" if study.similarity else "report")
    for root, key in ((matrix, "matrix"), (report, "report")):
        if sha256(root / "manifest.json") != receipt[f"{key}_manifest_sha256"]:
            raise ValueError(f"{key} manifest differs from the accepted receipt for {source.name}")
    fits = [fit for fit in read_json(matrix / "manifest.json")["fits"]
            if fit.get("metric") == study.similarity]
    manifest = read_json(report / "manifest.json")
    expected = {(arm, fold) for arm in study.arms for fold in range(3)}
    if (receipt["status"] != "PASS" or manifest["status"] != "PASS" or manifest["seed"] != 42
            or manifest["runs"] != len(expected) or len(fits) != len(expected) or receipt["fits"] != len(expected)
            or {(_arm(fit), fit["fold"]) for fit in fits} != expected):
        raise ValueError(f"accepted fit roster differs for {source.name}")
    if manifest.get("reference_alpha", study.reference) != study.reference:
        raise ValueError(f"accepted reference differs for {source.name}")
    if check_outputs(report, manifest) != receipt["report_outputs_sha256"]:
        raise ValueError(f"accepted report table hashes differ for {source.name}")
    return matrix, report, fits


def _score_fit(run, fit, record, sealed, scores, epochs, data):
    arm, fold, seed = _arm(fit), fit["fold"], 42
    identity = arm, fold, seed
    registration_path, manifest_path = (run / name for name in ("registration.json", "manifest.json"))
    config, manifest = map(read_json, (registration_path, manifest_path))
    if (sha256(registration_path) != fit["registration_sha256"] or record.run_sha256 != fit["registration_sha256"]
            or sha256(manifest_path) != fit["run_manifest_sha256"] or record.manifest_sha256 != fit["run_manifest_sha256"]
            or record.trace_sha256 != manifest["outputs"]["learning_curve.csv"]
            or manifest["registration_sha256"] != record.run_sha256
            or (manifest["fold"], manifest["seed"], manifest["axis"], manifest["smoke"])
            != (fold, seed, "denoise-5k", False)
            or any(manifest[key] != "exact" for key in ("test_replay", "validation_replay",
                                                      "validation_objective_replay"))):
        raise ValueError(f"fit or replay identity differs for {identity}")

    curve_path = run / "learning_curve.csv"
    if sha256(curve_path) != record.trace_sha256:
        raise ValueError(f"sealed trace digest differs for {identity}")
    curve = _table(curve_path).loc[:, list(TRACE)]
    if (record.stopping_epoch != fit["stopping_epoch"] or record.best_epoch != fit["best_epoch"]
            or curve.epoch.tolist() != list(range(1, record.stopping_epoch + 1))
            or not curve.equals(epochs.loc[:, list(TRACE)].reset_index(drop=True))):
        raise ValueError(f"accepted train and validation traces differ for {identity}")

    conditions, top200, weighted = read_fixed_jaccard(run, config, manifest, data)
    if (conditions != tuple(sealed.condition_id)
            or not np.array_equal(top200, sealed.JAC_top200, equal_nan=True)):
        raise ValueError(f"held-out JAC-top200 differs from accepted conditions for {identity}")
    scores = scores.set_index("metric", verify_integrity=True)
    if not set(OLD_METRICS) <= set(scores.index):
        raise ValueError(f"accepted fold metric roster differs for {identity}")
    jacc = scores.loc["JAC_top200"]
    observed = summarize_conditions(top200)
    if (not np.isclose(jacc.value, observed["value"], atol=1e-15, rtol=0, equal_nan=True)
            or any(observed[key] != jacc[key] for key in ("n_conditions", "n_eligible", "n_undefined"))):
        raise ValueError(f"held-out JAC-top200 fold score differs for {identity}")
    weighted_mean = summarize_conditions(weighted)["value"]
    return (
        [(*identity, *row) for row in curve.itertuples(index=False, name=None)],
        [(*identity, metric, weighted_mean if metric == "wJAC-2000" else scores.at[metric, "value"])
         for metric in METRIC_NAMES],
    )


def _write_plot_data(matrix, report, fits, data, output):
    runs = _table(report / "runs.csv")
    if len(runs) != len(fits) or runs.run_sha256.duplicated().any():
        raise ValueError(f"invalid accepted fit ledger in {report}")
    records = {record.run_sha256: record for record in runs.itertuples(index=False)}
    tables = [_table(report / name).groupby("run_sha256", sort=False)
              for name in ("condition-scores.csv", "fold-scores.csv", "epoch-objectives.csv")]
    if any(set(group.groups) != records.keys() for group in tables):
        raise ValueError(f"accepted score and trace rosters differ in {report}")

    curves, metrics = [], []
    for fit in fits:
        identity = fit["registration_sha256"]
        record = records[identity]
        if (_arm(record._asdict()), record.fold, record.seed) != (_arm(fit), fit["fold"], 42):
            raise ValueError(f"accepted report coordinate differs for {identity}")
        matched = [group.get_group(identity) for group in tables]
        for table in matched:
            if (not table.fold.eq(record.fold).all() or not table.seed.eq(record.seed).all()
                    or not table.alpha.eq(record.alpha).all()
                    or ("loss" in table and not table.loss.eq(record.loss).all())):
                raise ValueError(f"accepted report table coordinate differs for {identity}")
        fit_curves, fit_metrics = _score_fit(matrix / fit["path"], fit, record, *matched, data)
        curves.extend(fit_curves)
        metrics.extend(fit_metrics)

    output.mkdir()
    pd.DataFrame(curves, columns=CURVE_FIELDS).to_csv(output / "learning-curve.csv", index=False)
    pd.DataFrame(metrics, columns=METRIC_FIELDS).to_csv(output / "metrics.csv", index=False)


def build_study(source, output, data):
    study = STUDIES[source.name]
    step = LOSS / source.name
    receipt_path = step / "log" / study.receipt
    receipt = read_json(receipt_path)
    matrix, report, fits = _accepted_sources(source, receipt, study)
    output.mkdir()
    figure = output / "figure"
    figure.mkdir()
    plotted = figure / "data"
    _write_plot_data(matrix, report, fits, data, plotted)
    write_paired_metrics(plotted, study.reference)

    renderer = step / RENDERER
    dependencies = (LOSS / "figure/script/loss_comparison.py",)
    evidence = {
        "acceptance_receipt_sha256": sha256(receipt_path),
        "matrix_manifest_sha256": receipt["matrix_manifest_sha256"],
        "report_manifest_sha256": receipt["report_manifest_sha256"],
        "build_sha256": sha256(Path(__file__)),
        "matched_controls_sha256": data.release.manifest_sha256,
        "fixed_jaccard_code_sha256": sha256(STAGE / "04_tsg_axis/output/api/tsg_axis/comparison_evidence.py"),
        "shared_metrics_sha256": sha256(STAGE / "util/output/api/cpst_util/metrics.py"),
        "sealed_jac_top200_parity": True,
        "fits": len(fits),
        "reference": study.reference,
    }
    write_json(output / "log/figure-render.json", {
        **render_receipt(STAGE, step.name, output, renderer, dependencies=dependencies),
        "evidence": evidence,
    })
    print(json.dumps({"status": "PASS", "step": step.name, "fits": len(fits)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", type=Path, default=STAGE.parent,
                        help="materialized checkout containing the accepted Step 03 and Step 05 releases")
    parser.add_argument("--step", action="append", choices=STUDIES, help="rebuild only the named Step 05 comparison")
    parser.add_argument("--output", type=Path, required=True, help="new candidate directory under tmp/")
    args = parser.parse_args()
    candidate = args.output.resolve()
    candidate.mkdir()
    evidence = args.evidence_root.resolve() / STAGE.name
    with open_tenv_training(tenv_axis_output("denoise-5k", evidence / "03_data_train")) as release:
        data = AxisInputs(release)
        for name in args.step or STUDIES:
            build_study(evidence / LOSS.name / name, candidate / name, data)


if __name__ == "__main__":
    main()
