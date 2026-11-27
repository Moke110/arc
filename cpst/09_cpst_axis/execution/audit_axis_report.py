"""Independently reconstruct and check Step 09 full-axis scores and their summaries."""

from collections import defaultdict
from hashlib import sha256
import json
import math

import numpy as np
from scipy.spatial.distance import pdist

from axis_protocol import ARM, REFERENCE_AXIS, coordinates
from axis_report import PLOT_TABLES
from axis_scores import REPORT_CONDITION, REPORT_METRICS, WEIGHTED
from axis_runs import axis_run
from audit_measured_scores import check_statistics, read_table, same, value
from chem_runs import SCHEMA as STEP08_SCHEMA
from tenv_axes import open_tenv_axis_training
from cpst_util.comparison import METRIC_NAMES
from tsg_axis.metrics import DISTRIBUTION, METRICS, condition_scores, distribution_scores, signed_jaccard

IDENTITY = ("axis", "fold", "seed")
SCORE_FIELDS = (*IDENTITY, "metric")


def audit_aggregates(scores, summary):
    """Rejoin each replicate to the reference and audit equal-weight means and sample SD."""
    groups = defaultdict(list)
    for (axis, fold, seed, metric), row in scores.items():
        reference = value(scores[REFERENCE_AXIS, fold, seed, metric])
        score = value(row)
        difference = score - reference if score is not None and reference is not None else None
        if not same(value(row, "reference_value"), reference) or not same(value(row, "paired_difference"), difference):
            raise ValueError(f"incorrect fold/seed pair: {(axis, fold, seed, metric)}")
        groups[axis, metric].append((score, difference, row))

    if set(summary) != set(groups):
        raise ValueError("summary has incorrect axis or metric groups")
    for key, replicates in groups.items():
        group = summary[key]
        values, differences, rows = zip(*replicates, strict=True)
        check_statistics(group, values, key, field="mean")
        check_statistics(group, differences, key, field="paired_mean", sd="paired_sd", count="n_paired")
        if ({int(row["width"]) for row in rows} != {int(group["width"])}
                or any(int(group[field]) != sum(int(row[field]) for row in rows)
                       for field in ("n_eligible", "n_undefined"))):
            raise ValueError(f"incorrect aggregate support: {key}")


def treatment_rows(release, rows, layer):
    """Read a registered axis in small contiguous windows, maintaining its row order."""
    parts = []
    for start in range(0, len(release.treatments), 128):
        stop = min(start + 128, len(release.treatments))
        offsets = rows[(rows >= start) & (rows < stop)] - start
        if len(offsets):
            parts.append(release.read_treatments(start, stop, layer)[offsets])
    return np.concatenate(parts)


def ordered_condition_ids(release, rows):
    columns = release.treatments.iloc[rows][["plate", "drug_parent", "dose_um"]]
    return [json.dumps([plate, parent, float(dose)], ensure_ascii=False, separators=(",", ":"))
            for plate, parent, dose in columns.itertuples(index=False, name=None)]


def audit_bandwidths(report, step03, step04):
    """Check row and gene order, training Delta and the full-axis bandwidth for each fold."""
    calibrations = read_table(report / "bandwidths.csv", ("fold",))
    original_path = step04 / "output/report-matched-control/bandwidths.csv"
    sealed = axis_run.read_json(step04 / "log/319-formal.json")["report"]["output_sha256"]["bandwidths.csv"]
    if axis_run.digest(original_path) != sealed:
        raise ValueError("original Step 04 full-axis MMD calibration is not sealed")
    original = read_table(original_path, ("fold",))
    if set(calibrations) != {(fold,) for fold in range(3)}:
        raise ValueError("MMD calibrations must cover all three folds")
    result = {}
    with open_tenv_axis_training("hvg-5k", step03) as release:
        genes = release.genes.gene_index.to_numpy()
        for fold in range(3):
            rows = np.flatnonzero(release.treatments[f"cold_drug_{fold}"].to_numpy() == "train")
            matrix = treatment_rows(release, rows, "delta")
            record = calibrations[fold,]
            distance = np.median(pdist(matrix)) / math.sqrt(matrix.shape[1])
            if (matrix.shape != (int(record["rows"]), int(record["width"]))
                    or int(record["pair_count"]) != math.comb(len(rows), 2)
                    or tuple(map(int, record["train_rows"].split(";"))) != tuple(rows)
                    or tuple(map(int, record["gene_indices"].split(";"))) != tuple(genes)
                    or record["training_manifest_sha256"] != release.manifest_sha256
                    or record["training_delta_values_sha256"] != sha256(matrix.tobytes()).hexdigest()
                    or not same(value(record, "bandwidth"), distance)
                    or any(record[field] != prior for field, prior in original[fold,].items())):
                raise ValueError(f"incorrect frozen MMD calibration: fold {fold}")
            result[fold] = value(record, "bandwidth")
    return result


def audit_fold(key, genes, truth, prediction, identities, conditions, scores, bandwidth):
    """Rescore the complete axis and account for every condition and eligibility mask."""
    observed = defaultdict(list)
    for condition_id, measured, estimated in zip(identities, truth, prediction, strict=True):
        row = conditions.pop((*key, condition_id))
        recalculated = condition_scores(measured, estimated, genes)
        cutoffs = ((100, 1), (250, 2), (500, 4), (1000, 8))
        recalculated[WEIGHTED] = sum(
            weight * signed_jaccard(measured, estimated, genes, limit=cutoff)[0]
            for cutoff, weight in cutoffs
        ) / 15
        for metric in REPORT_CONDITION:
            number = recalculated[metric]
            expected = float(number) if math.isfinite(number) else None
            if not same(value(row, metric), expected):
                raise ValueError(f"incorrect held-out condition score: {(*key, condition_id, metric)}")
            if expected is not None:
                observed[metric].append(expected)

    distances = dict(zip(DISTRIBUTION, distribution_scores(truth, prediction, bandwidth), strict=True))
    for metric in REPORT_METRICS:
        record = scores[*key, metric]
        if metric in distances:
            expected, eligible = distances[metric], len(identities)
        else:
            eligible = len(observed[metric])
            expected = math.fsum(observed[metric]) / eligible if eligible else None
        if (not same(value(record), expected) or int(record["n_conditions"]) != len(identities)
                or int(record["n_eligible"]) != eligible or int(record["n_undefined"]) != len(identities) - eligible
                or int(record["width"]) != len(genes)):
            raise ValueError(f"incorrect independently rescored replicate: {(*key, metric)}")


def audit_plotted_csvs(data, report, runs):
    """Verify the plotted values against audited scores and the sealed per-fit traces."""
    if {path.name for path in data.iterdir() if path.is_file()} != set(PLOT_TABLES):
        raise ValueError("plotted CSV inventory differs from the canonical Step 09 figures")
    scores = {key: row for key, row in read_table(report / "fold-scores.csv", SCORE_FIELDS).items()
              if key[-1] in METRIC_NAMES}
    metric_keys = ("arm", "fold", "seed", "metric")
    raw = read_table(data / "metrics.csv", metric_keys)
    paired = read_table(data / "metrics-paired.csv", metric_keys)
    if set(raw) != set(scores) or any(not same(value(raw[key]), value(row)) for key, row in scores.items()):
        raise ValueError("metrics.csv plotted values differ from the independently audited scores")
    comparisons = {key: row for key, row in scores.items() if key[0] != REFERENCE_AXIS}
    if set(paired) != set(comparisons) or any(
        plotted["reference"] != REFERENCE_AXIS
        or not same(value(plotted, "arm_value"), value(row))
        or not same(value(plotted, "reference_value"), value(row, "reference_value"))
        or not same(value(plotted), value(row, "paired_difference"))
        for key, row in comparisons.items() for plotted in (paired[key],)
    ):
        raise ValueError("metrics-paired.csv differs from the audited fold and seed pairs")

    objectives = {"train": "train_loss", "validation": "validation_loss"}
    reported = read_table(report / "epoch-objectives.csv", (*IDENTITY, "epoch"))
    plotted = read_table(data / "learning-curve.csv", ("arm", "fold", "seed", "epoch"))
    if set(plotted) != set(reported) or any(
        not same(value(plotted[key], field), value(row, source))
        for key, row in reported.items() for field, source in objectives.items()
    ):
        raise ValueError("learning-curve plotted values differ from the reported objectives")
    for key, (path, *_) in runs.items():
        for (epoch,), row in read_table(path / "learning_curve.csv", ("epoch",)).items():
            coordinate = (*key, epoch)
            selected = plotted.pop(coordinate, None)
            if selected is None or any(not same(value(selected, field), value(row, source))
                                       for field, source in objectives.items()):
                raise ValueError(f"learning-curve objectives differ from the sealed fit trace: {coordinate}")
    if plotted:
        raise ValueError("learning-curve has rows beyond the sealed fit traces")


def audit_report(report, runs, step03, step04, step08, *, bounded=False):
    """Check the report against immutable fits, ordered axes and independent matched controls."""
    expected_roster = coordinates(bounded)
    scores = read_table(report / "fold-scores.csv", SCORE_FIELDS)
    summary = read_table(report / "summary.csv", ("axis", "metric"))
    conditions = read_table(report / "condition-scores.csv", (*IDENTITY, "condition_id"))
    ledger = read_table(report / "runs.csv", IDENTITY)
    keys = {(*key, metric) for key in expected_roster for metric in REPORT_METRICS}
    if set(runs) != set(expected_roster) or set(ledger) != set(expected_roster) or set(scores) != keys:
        raise ValueError("independent scoring requires the complete coordinate and metric roster")
    bandwidths = audit_bandwidths(report, step03, step04)
    original_scores = {}
    if any(run[1]["schema"] == STEP08_SCHEMA for run in runs.values()):
        original_path = "figure/data/fold-scores.csv"
        original_acceptance = axis_run.read_json(step08 / "log/364-formal.json")
        if axis_run.digest(step08 / original_path) != original_acceptance["release_files_sha256"][original_path]:
            raise ValueError("the historical Step 08 score table is not sealed")
        original_scores = read_table(step08 / original_path, ("arm", "fold", "seed", "metric"))
    initial_conditions = len(conditions)
    for key in expected_roster:
        axis, fold, seed = key
        path, config, replay, _, _ = runs[key]
        origin = "Step 08" if config["schema"] == STEP08_SCHEMA else "Step 09"
        if (ledger[key]["origin"] != origin or ledger[key]["registration_sha256"] != replay["registration_sha256"]):
            raise ValueError(f"incorrect scored run lineage: {key}")
        with open_tenv_axis_training(axis, step03) as release:
            gene_indices = release.genes.gene_index.to_numpy()
            rows = np.flatnonzero(release.treatments[f"cold_drug_{fold}"].to_numpy() == "test")
            if bounded:
                rows = rows[:replay["test_conditions"]]
            ids = ordered_condition_ids(release, rows)
            linked = release.treatments.iloc[rows].control_row.to_numpy()
            controls = release.read_controls(0, len(release.controls))[linked]
            truth = treatment_rows(release, rows, "treatment") - controls
            prediction = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False) - controls
            if (len(ids) != replay["test_conditions"]
                    or tuple(read_table(path / "conditions.csv", ("condition_id",))) != tuple((name,) for name in ids)
                    or not np.array_equal(truth, np.load(path / "measured_delta.npy", mmap_mode="r", allow_pickle=False))
                    or not np.array_equal(prediction, np.load(path / "predicted_delta.npy", mmap_mode="r", allow_pickle=False))):
                raise ValueError(f"independent held-out Delta differs: {key}")
        audit_fold(key, gene_indices, truth, prediction, ids, conditions, scores, bandwidths[fold])
        if origin == "Step 08":
            for metric in METRICS:
                original = original_scores[ARM, fold, seed, metric]
                current = scores[*key, metric]
                if (value(original) != value(current)
                        or any(original[field] != current[field] for field in ("n_conditions", "n_eligible", "n_undefined"))):
                    raise ValueError(f"reused full-axis scores differ from Step 08: {(*key, metric)}")
    if conditions:
        raise ValueError("unexpected held-out conditions in the formal report")
    audit_aggregates(scores, summary)
    return {"status": "PASS", "coordinates": len(expected_roster), "score_rows": len(scores),
            "condition_rows": initial_conditions, "summary_rows": len(summary), "fold_bandwidths": len(bandwidths)}
