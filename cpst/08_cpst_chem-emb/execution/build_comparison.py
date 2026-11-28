"""Build the Step 08 comparison figures from its sealed matrix and KPGT reuses."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_chem_report import audit
from chem_figures import replay
from chem_formal import FIGURE_RECEIPT, PROJECTIONS, RECEIPT, checked_report, figure_receipt
from chem_runs import DEFAULT_TRAINING, REFERENCE, ROUND_FOUR, STAGE, STEP, coordinates
from cpst_util.comparison import CURVE_FIELDS, METRIC_FIELDS, METRIC_NAMES, write_paired_metrics
from cpst_util.derivation.release import read_json, sha256, write_json
from training_data import open_tenv_training
from tsg_axis.comparison_evidence import read_fixed_jaccard
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import summarize_conditions
from util.training.cpst_reference import verify_reference, verify_reference_release, verify_release

TRACE = ["epoch", "train_loss", "validation_loss"]
IDENTITY = ["arm", "fold", "seed"]


def read_table(path):
    return pd.read_csv(path, float_precision="round_trip")


def collect_accepted(bundle, round_four, receipt):
    matrix = bundle / "matrix"
    report = bundle / "report"
    plan_path = bundle / "reuse-preflight.json"
    paths = (
        (matrix / "manifest.json", "matrix_sha256"),
        (matrix / "preflight.json", "preflight_sha256"),
        (report / "manifest.json", "report_sha256"),
        (plan_path, "reuse_plan_sha256"),
    )
    if (receipt["status"], receipt["ticket"], receipt["total_coordinates"]) != ("PASS", 364, len(coordinates())):
        raise ValueError("invalid Step 08 acceptance receipt")
    for path, name in paths:
        if sha256(path) != receipt[name]:
            raise ValueError(f"accepted digest differs for {path}")
    report_manifest = checked_report(report)
    if any(report_manifest["csv_sha256"][name] != receipt["release_files_sha256"][published]
           or sha256(STEP / published) != receipt["release_files_sha256"][published]
           for name, published in PROJECTIONS.items()):
        raise ValueError("published Step 08 scientific tables differ from their accepted report")
    audit(report)

    fits = read_json(matrix / "manifest.json")["fits"]
    plan = read_json(plan_path)
    fresh = {tuple(fit["coordinate"]): fit for fit in fits}
    reuses = {tuple(row["coordinate"]): row for row in plan["decisions"]}
    expected_reuses = {key for key in coordinates() if key[0] == REFERENCE}
    if (len(fits) != len(fresh) or len(fresh) != receipt["new_exact_independent_replays"]
            or len(reuses) != len(plan["decisions"]) or len(reuses) != receipt["reused_kpgt_coordinates"]
            or plan["decisions"] != receipt["reference_decisions"]
            or set(reuses) != expected_reuses or set(fresh) != set(coordinates()) - expected_reuses
            or any(row["decision"] != "reuse" for row in reuses.values())):
        raise ValueError("accepted Step 08 fit and KPGT reuse roster differs")
    release = verify_release()
    if (plan["alpha_release_sha256"], plan["alpha_receipt_sha256"], plan["step_06_4_receipt_sha256"]) != (
        receipt["alpha_release_sha256"], receipt["alpha_receipt_sha256"], release["step_06_4"]["receipt_sha256"],
    ):
        raise ValueError("KPGT reuse plan differs from the fixed alpha release")
    ledger = verify_reference_release(round_four, release)
    for key in coordinates():
        arm, fold, seed = key
        if key in reuses:
            decision = reuses[key]
            path, config, replayed, _ = verify_reference(ledger, release, fold, seed)
            if (replayed["registration_sha256"] != decision["registration_sha256"]
                    or any(replayed["outputs"][name] != decision[field] for name, field in (
                        ("checkpoint.pt", "checkpoint_sha256"), ("prediction.npy", "prediction_sha256"),
                        ("validation_prediction.npy", "validation_prediction_sha256"),
                    ))):
                raise ValueError(f"Step 06.4 fit differs from accepted KPGT reuse {key}")
        else:
            fit = fresh[key]
            relative = Path("runs") / arm / f"fold-{fold}" / f"seed-{seed}"
            path = matrix / relative
            if (fit["path"] != relative.as_posix()
                    or sha256(path / "registration.json") != fit["registration_sha256"]
                    or sha256(path / "manifest.json") != fit["manifest_sha256"]):
                raise ValueError(f"fresh Step 08 fit differs from the accepted matrix: {key}")
            config = read_json(path / "registration.json")
            replayed = read_json(path / "manifest.json")
        yield key, (path, config, replayed, key in reuses)


def score_fit(key, fit, record, sealed, folds, epochs, data):
    path, config, manifest, reused = fit
    outputs = manifest["outputs"]
    expected = {
        "registration_sha256": manifest["registration_sha256"],
        "manifest_sha256": sha256(path / "manifest.json"),
        "fit_origin": "Step 06.4" if reused else "Step 08",
        "checkpoint_sha256": outputs["checkpoint.pt"],
        "prediction_sha256": outputs["prediction.npy"],
        "trace_sha256": outputs["learning_curve.csv"],
    }
    if (any(record[field] != value for field, value in expected.items())
            or (manifest["fold"], manifest["seed"], manifest["smoke"]) != (*key[1:], False)
            or any(manifest[name] != "exact" for name in (
                "test_replay", "validation_replay", "validation_objective_replay",
            ))):
        raise ValueError(f"Step 08 fit or exact replay differs from its report ledger: {key}")
    for table in (sealed, folds, epochs):
        if not table.registration_sha256.eq(record.registration_sha256).all():
            raise ValueError(f"Step 08 report rows use a different fit: {key}")

    trace_path = path / "learning_curve.csv"
    trace = read_table(trace_path)[TRACE]
    if (sha256(trace_path) != record.trace_sha256
            or trace.epoch.tolist() != list(range(1, record.stopping_epoch + 1))
            or not trace.equals(epochs[TRACE].reset_index(drop=True))):
        raise ValueError(f"accepted train and validation trace differs for {key}")
    conditions, top200 = read_fixed_jaccard(path, config, manifest, data)
    if (conditions != tuple(sealed.condition_id)
            or not np.array_equal(top200, sealed.JAC_top200, equal_nan=True)):
        raise ValueError(f"held-out JAC-top200 differs from the accepted conditions: {key}")
    folds = folds.set_index("metric", verify_integrity=True)
    existing = folds.loc["JAC_top200"]
    summary = summarize_conditions(top200)
    if (not np.isclose(existing.value, summary["value"], rtol=0, atol=1e-15, equal_nan=True)
            or any(existing[name] != summary[name] for name in ("n_conditions", "n_eligible", "n_undefined"))):
        raise ValueError(f"held-out JAC-top200 fold score differs for {key}")
    return (
        [(*key, *row) for row in trace.itertuples(index=False, name=None)],
        [(*key, metric, folds.at[metric, "value"])
         for metric in METRIC_NAMES],
    )


def build(bundle, round_four, training, output, scratch):
    receipt_path = STEP / RECEIPT
    accepted = read_json(receipt_path)
    report = bundle / "report"
    records = read_table(report / "runs.csv").set_index(IDENTITY, verify_integrity=True)
    groups = [read_table(report / name).groupby(IDENTITY, sort=False)
              for name in ("condition-scores.csv", "fold-scores.csv", "epoch-objectives.csv")]
    expected = set(coordinates())
    if set(records.index) != expected or any(set(table.groups) != expected for table in groups):
        raise ValueError("Step 08 report misses accepted fit, trace or score coordinates")

    curves, metrics = [], []
    with open_tenv_training(training) as release:
        data = AxisInputs(release)
        training_sha = release.manifest_sha256
        for key, fit in collect_accepted(bundle, round_four, accepted):
            curve, scores = score_fit(key, fit, records.loc[key], *(table.get_group(key) for table in groups), data)
            curves.extend(curve)
            metrics.extend(scores)

    output.mkdir()
    figure = output / "figure"
    plotted = figure / "data"
    plotted.mkdir(parents=True)
    pd.DataFrame(curves, columns=CURVE_FIELDS).to_csv(plotted / "learning-curve.csv", index=False)
    pd.DataFrame(metrics, columns=METRIC_FIELDS).to_csv(plotted / "metrics.csv", index=False)
    write_paired_metrics(plotted, REFERENCE, mode="matched-coordinate")
    picture = replay(plotted, figure, scratch, comparison=True)
    write_json(output / FIGURE_RECEIPT, {
        **figure_receipt(picture, sha256(receipt_path), comparison=True),
        "evidence": {
            "matrix_manifest_sha256": accepted["matrix_sha256"],
            "report_manifest_sha256": accepted["report_sha256"],
            "reuse_plan_sha256": accepted["reuse_plan_sha256"],
            "step_06_4_receipt_sha256": sha256(round_four / "log/359-formal.json"),
            "build_sha256": sha256(Path(__file__)),
            "matched_controls_sha256": training_sha,
            "fixed_jaccard_code_sha256": sha256(STAGE / "04_tsg_axis/output/api/tsg_axis/comparison_evidence.py"),
            "shared_metrics_sha256": sha256(STAGE / "util/output/api/cpst_util/metrics.py"),
            "sealed_jac_top200_parity": True, "fits": len(expected), "reused_kpgt": len(accepted["reference_decisions"]),
            "epoch_rows": len(curves), "reference": REFERENCE,
        },
    })
    return {"status": "PASS", "fits": len(expected), "epoch_rows": len(curves), "figures": list(picture["png_sha256"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=STEP / "output/formal")
    parser.add_argument("--round-four", type=Path, default=ROUND_FOUR)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), sort_keys=True))


if __name__ == "__main__":
    main()
