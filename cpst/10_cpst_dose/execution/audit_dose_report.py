"""Independently rescore Step 10 and check its published tables."""

import argparse
from collections import defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial.distance import pdist

from dose_inputs import SMOKE_LIMITS, cold_dose_partitions
from dose_matrix import verify_matrix
from dose_protocol import STAGE, AXIS, DEFAULT_STEP03, REFERENCE, TRANSIGEN, coordinates, fold_identity
from dose_report import PLOT_TABLES, REPORT_METRICS
from dose_runs import axis_run
sys.path.insert(0, str(STAGE / "09_cpst_axis/execution"))

from audit_axis_report import (audit_fold, check_statistics, ordered_condition_ids, read_table, same,
                              treatment_rows, value)
from cpst_util.comparison import METRIC_NAMES
from training_data import open_tenv_axis_training

IDENTITY = ("arm", "fold", "seed")
SCORE_KEY = (*IDENTITY, "metric")
OBJECTIVES = {"train": "train_loss", "validation": "validation_loss"}


def audit_aggregates(scores, foldwise, overall):
    """Verify no-dose pairing and the equal-weight statistics at both levels."""
    by_fold, by_arm = defaultdict(list), defaultdict(list)
    for (arm, fold, seed, metric), row in scores.items():
        reference = None if arm == TRANSIGEN else value(scores[REFERENCE, fold, seed, metric])
        observed = value(row)
        difference = None if observed is None or reference is None else observed - reference
        if not same(value(row, "reference_value"), reference) or not same(value(row, "paired_difference"), difference):
            raise ValueError(f"incorrect fold and seed pairing for {(arm, fold, seed, metric)}")
        by_fold[arm, fold, metric].append((observed, difference, row))
        by_arm[arm, metric].append((observed, difference, row))

    for summaries, groups in ((foldwise, by_fold), (overall, by_arm)):
        if set(summaries) != set(groups):
            raise ValueError("summary rows differ from the registered arms, folds or metrics")
        for key, results in groups.items():
            observed, differences, rows = zip(*results, strict=True)
            summary = summaries[key]
            check_statistics(summary, observed, key, field="mean")
            check_statistics(summary, differences, key, field="paired_mean", sd="paired_sd", count="n_paired")
            if any(int(summary[field]) != sum(int(row[field]) for row in rows)
                   for field in ("n_eligible", "n_undefined")):
                raise ValueError(f"incorrect aggregate eligibility for {key}")


def matched_deltas(release, controls, rows):
    selected = release.treatments.iloc[rows]
    linked = selected.control_row.to_numpy()
    if not np.array_equal(release.controls.plate.to_numpy()[linked], selected.plate.to_numpy()):
        raise ValueError("frozen treatment and DMSO control plates differ")
    control = controls[linked]
    return control, treatment_rows(release, rows, "treatment") - control


def check_fold_identity(row, fold):
    return all(row[field] == str(expected) for field, expected in fold_identity(fold).items())


def audit_partitions(treatments, fold):
    identity = fold_identity(fold)
    old_roles = treatments[f"cold_drug_{identity['parent_group']}"]
    held_parents = old_roles.eq("test")
    held_dose = treatments.dose_um.eq(identity["held_dose_um"])
    neither = ~held_parents & ~held_dose
    masks = {"train": neither & old_roles.eq("train"), "valid": neither & old_roles.eq("valid"),
             "test": held_parents & held_dose, "excluded-held-parent": held_parents & ~held_dose,
             "excluded-held-dose": ~held_parents & held_dose}
    rows = {phase: np.flatnonzero(mask.to_numpy()) for phase, mask in masks.items()}
    declared = cold_dose_partitions(treatments, fold)
    if any(not np.array_equal(actual, declared[phase]) for phase, actual in rows.items()):
        raise ValueError(f"independent crossed partition audit disagrees for fold {fold}")
    return rows


def audit_bandwidth(release, controls, fold, record, genes, rows):
    _, delta = matched_deltas(release, controls, rows)
    bandwidth = np.median(pdist(delta)) / math.sqrt(len(genes))
    identities = release.treatments.iloc[rows][["plate", "drug_parent", "dose_um"]].copy()
    identities.insert(0, "condition_id", ordered_condition_ids(release, rows))
    if (not check_fold_identity(record, fold)
            or delta.shape != (int(record["rows"]), int(record["width"]))
            or int(record["pair_count"]) != math.comb(len(rows), 2)
            or tuple(map(int, record["train_rows"].split(";"))) != tuple(rows)
            or tuple(map(int, record["gene_indices"].split(";"))) != tuple(genes)
            or record["training_manifest_sha256"] != release.manifest_sha256
            or record["training_delta_values_sha256"] != sha256(delta.tobytes()).hexdigest()
            or record["training_conditions_sha256"] != sha256(identities.to_csv(index=False).encode()).hexdigest()
            or not same(value(record, "bandwidth"), bandwidth)):
        raise ValueError(f"incorrect frozen MMD calibration for fold {fold}")
    return bandwidth


def audit_lineage(key, run, row):
    path, config, replay, trained, _ = run
    identity = {
        **fold_identity(key[1]), "registration_sha256": replay["registration_sha256"],
        "training_sha256": config["training_sha256"], "axis_sha256": config["axis_sha256"],
        "representation": config["representation"]["name"],
        "representation_sha256": config["representation"]["vectors_sha256"],
        "manifest_sha256": axis_run.digest(path / "manifest.json"),
        "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
        "validation_prediction_sha256": replay["outputs"]["validation_prediction.npy"],
        "prediction_sha256": replay["outputs"]["prediction.npy"],
        "trace_sha256": replay["outputs"]["learning_curve.csv"],
        "best_epoch": trained["best_epoch"], "epochs": trained["epochs"], "parameters": trained["parameters"],
    }
    if any(row[field] != str(expected) for field, expected in identity.items()):
        raise ValueError(f"incorrect registered fit lineage for {key}")


def audit_plots(data, report, runs, scores):
    """Check all CSV-only inputs against audited scores and sealed learning traces."""
    if {path.name for path in data.iterdir()} != set(PLOT_TABLES):
        raise ValueError("incorrect Step 10 plot table inventory")
    scores = {key: row for key, row in scores.items() if key[-1] in METRIC_NAMES}
    raw = read_table(data / "metrics.csv", SCORE_KEY)
    paired = read_table(data / "metrics-paired.csv", SCORE_KEY)
    if set(raw) != set(scores) or any(not same(value(raw[key]), value(row)) for key, row in scores.items()):
        raise ValueError("plotted absolute metrics differ from audited scores")
    comparisons = {key: row for key, row in scores.items() if key[0] not in (REFERENCE, TRANSIGEN)}
    if set(paired) != set(comparisons) or any(
        plotted["reference"] != REFERENCE
        or not same(value(plotted, "arm_value"), value(row))
        or not same(value(plotted, "reference_value"), value(row, "reference_value"))
        or not same(value(plotted), value(row, "paired_difference"))
        for key, row in comparisons.items() for plotted in (paired[key],)
    ):
        raise ValueError("plotted paired metrics differ from audited no-dose contrasts")

    indexes = (*IDENTITY, "epoch")
    reported = read_table(report / "epoch-objectives.csv", indexes)
    curves = read_table(data / "learning-curve.csv", indexes)
    if set(curves) != set(reported):
        raise ValueError("plotted and reported curve coordinates differ")
    for key, (path, _, replay, _, _) in runs.items():
        for (epoch,), trace in read_table(path / "learning_curve.csv", ("epoch",)).items():
            coordinate = (*key, epoch)
            plotted = curves.pop(coordinate, None)
            row = reported[coordinate]
            if (plotted is None or row["registration_sha256"] != replay["registration_sha256"]
                    or not check_fold_identity(row, key[1])
                    or any(not same(value(row, source), value(trace, source)) for source in
                           (*OBJECTIVES.values(), "validation_mse"))
                    or any(not same(value(plotted, field), value(trace, source))
                           for field, source in OBJECTIVES.items())):
                raise ValueError(f"learning curve differs from the sealed trace at {coordinate}")
    if curves:
        raise ValueError("plotted learning curves include unregistered epochs")


def audit_report(report, data, runs, step03, *, bounded=False):
    """Reconstruct the complete held-out score table before acceptance."""
    roster = coordinates(bounded)
    scores = read_table(report / "fold-scores.csv", SCORE_KEY)
    conditions = read_table(report / "condition-scores.csv", (*IDENTITY, "condition_id"))
    ledger = read_table(report / "runs.csv", IDENTITY)
    calibrations = read_table(report / "bandwidths.csv", ("fold",))
    folds = {key[1] for key in roster}
    if (set(runs) != set(roster) or set(ledger) != set(roster)
            or set(scores) != {(*key, metric) for key in roster for metric in REPORT_METRICS}
            or set(calibrations) != {(fold,) for fold in folds}):
        raise ValueError("report inventory differs from the Step 10 roster")

    condition_count = len(conditions)
    with open_tenv_axis_training(AXIS, step03) as release:
        controls = release.read_controls(0, len(release.controls))
        genes = release.genes.gene_index.to_numpy()
        for fold in sorted(folds):
            partitions = audit_partitions(release.treatments, fold)
            bandwidth = audit_bandwidth(release, controls, fold, calibrations[fold,], genes, partitions["train"])
            rows = partitions["test"]
            if bounded:
                rows = rows[:SMOKE_LIMITS["test"]]
            ids = ordered_condition_ids(release, rows)
            control, truth = matched_deltas(release, controls, rows)
            for key in (key for key in roster if key[1] == fold):
                audit_lineage(key, runs[key], ledger[key])
                path, config, replay, _, _ = runs[key]
                for row in (scores[(*key, metric)] for metric in REPORT_METRICS):
                    if not check_fold_identity(row, fold) or row["registration_sha256"] != replay["registration_sha256"]:
                        raise ValueError(f"incorrect score lineage for {key}")
                for name in ids:
                    row = conditions[(*key, name)]
                    if not check_fold_identity(row, fold) or row["registration_sha256"] != replay["registration_sha256"]:
                        raise ValueError(f"incorrect condition score lineage for {key}")
                predicted = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False) - control
                expected_excluded = {reason: list(ordered_condition_ids(release, partitions[f"excluded-{reason}"]))
                                     for reason in ("held-parent", "held-dose")}
                if (config["population"]["excluded_conditions"] != expected_excluded
                        or tuple(config["partitions"]["test"]) != tuple(rows)
                        or tuple(config["conditions"]["test"]) != tuple(ids)
                        or config["population"]["control_sha256"]["test"] != sha256(control.tobytes()).hexdigest()
                        or len(ids) != replay["test_conditions"]
                        or tuple(read_table(path / "conditions.csv", ("condition_id",))) != tuple((name,) for name in ids)
                        or tuple(int(index) for index, in read_table(path / "genes.csv", ("gene_index",))) != tuple(genes)):
                    raise ValueError(f"independent ordered input audit disagrees at {key}")
                for name, actual in (("measured_delta", truth), ("predicted_delta", predicted)):
                    saved = np.load(path / f"{name}.npy", mmap_mode="r", allow_pickle=False)
                    if saved.dtype != actual.dtype or not np.array_equal(saved, actual):
                        raise ValueError(f"independent {name} calculation disagrees at {key}")
                audit_fold(key, genes, truth, predicted, ids, conditions, scores, bandwidth)
    if conditions:
        raise ValueError("Step 10 report contains unregistered held-out conditions")
    foldwise = read_table(report / "foldwise-summary.csv", ("arm", "fold", "metric"))
    overall = read_table(report / "summary.csv", ("arm", "metric"))
    if any(not check_fold_identity(row, key[1]) for key, row in foldwise.items()):
        raise ValueError("incorrect crossed fold summary lineage")
    audit_aggregates(scores, foldwise, overall)
    audit_plots(data, report, runs, scores)
    return {"status": "PASS", "coordinates": len(roster), "score_rows": len(scores), "condition_rows": condition_count,
            "foldwise_rows": len(foldwise), "summary_rows": len(overall), "fold_bandwidths": len(calibrations),
            "rtol": 1e-12, "atol": 1e-12}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step03", type=Path, default=DEFAULT_STEP03)
    parser.add_argument("--bounded", action="store_true")
    for name in ("matrix", "report", "figure-data"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    runs = verify_matrix(args.matrix, args.bounded)
    print(json.dumps(audit_report(args.report, args.figure_data, runs, args.step03, bounded=args.bounded), sort_keys=True))


if __name__ == "__main__":
    main()
