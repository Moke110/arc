"""Independently reconstruct and check both Step 09 scoring scopes and their summaries."""

from collections import defaultdict
from hashlib import sha256
import json
import math

import numpy as np
from scipy.spatial.distance import pdist

from axis_protocol import ARM, DISEASE_OVERLAP, DISEASE_SIZE, REFERENCE_AXIS, SCOPES, coordinates
from axis_scores import COVERAGE
from axis_runs import axis_run
from audit_measured_scores import check_statistics, read_table, same, value
from chem_runs import SCHEMA as STEP08_SCHEMA
from tenv_axes import open_tenv_axis_training
from tsg_axis.disease_reference import disease_gene_sets
from tsg_axis.metrics import CONDITION, DISTRIBUTION, METRICS, condition_scores, distribution_scores

IDENTITY = ("axis", "fold", "seed")
SCORE_FIELDS = (*IDENTITY, "scope", "metric")


def audit_aggregates(scores, summary):
    """Rejoin each replicate to the reference and audit equal-weight means and sample SD."""
    groups = defaultdict(list)
    for (axis, fold, seed, scope, metric), row in scores.items():
        reference = value(scores[REFERENCE_AXIS, fold, seed, scope, metric])
        score = value(row)
        difference = score - reference if score is not None and reference is not None else None
        if not same(value(row, "reference_value"), reference) or not same(value(row, "paired_difference"), difference):
            raise ValueError(f"incorrect fold/seed pair: {(axis, fold, seed, scope, metric)}")
        groups[axis, scope, metric].append((score, difference, row))

    if set(summary) != set(groups):
        raise ValueError("summary has incorrect axis, scope or metric groups")
    for key, replicates in groups.items():
        group = summary[key]
        values, differences, rows = zip(*replicates, strict=True)
        check_statistics(group, values, key, field="mean")
        check_statistics(group, differences, key, field="paired_mean", sd="paired_sd", count="n_paired")
        if ({int(row["overlap"]) for row in rows} != {int(group["overlap"])}
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


def audit_bandwidths(report, step03, step04, union):
    """Check stored row and gene order, raw training Δ bytes, and each fold/scope bandwidth."""
    calibrations = read_table(report / "bandwidths.csv", ("fold", "scope"))
    original = read_table(step04 / "output/report-matched-control/bandwidths.csv", ("fold",))
    expected = {(fold, scope) for fold in range(3) for scope in SCOPES}
    if set(calibrations) != expected:
        raise ValueError("MMD calibrations must cover three folds and two scopes")
    result = {}
    with open_tenv_axis_training("hvg-5k", step03) as release:
        genes = release.genes.gene_index.to_numpy()
        mask = np.asarray([gene in union for gene in genes])
        for fold in range(3):
            rows = np.flatnonzero(release.treatments[f"cold_drug_{fold}"].to_numpy() == "train")
            full = treatment_rows(release, rows, "delta")
            for scope in SCOPES:
                record = calibrations[fold, scope]
                matrix, scoped_genes = (full, genes) if scope == "axis" else (full[:, mask], genes[mask])
                distance = np.median(pdist(matrix)) / math.sqrt(matrix.shape[1])
                if (matrix.shape != (int(record["rows"]), int(record["width"]))
                        or int(record["pair_count"]) != math.comb(len(rows), 2)
                        or tuple(map(int, record["train_rows"].split(";"))) != tuple(rows)
                        or tuple(map(int, record["gene_indices"].split(";"))) != tuple(scoped_genes)
                        or record["training_manifest_sha256"] != release.manifest_sha256
                        or record["training_conditions_sha256"] != original[fold,]["training_conditions_sha256"]
                        or record["training_delta_values_sha256"] != sha256(matrix.tobytes()).hexdigest()
                        or not same(value(record, "bandwidth"), distance)
                        or scope == "axis" and any(record[field] != prior for field, prior in original[fold,].items())
                        or scope == "disease" and len(scoped_genes) != DISEASE_OVERLAP["hvg-5k"]):
                    raise ValueError(f"incorrect frozen MMD calibration: {(fold, scope)}")
                result[fold, scope] = value(record, "bandwidth")
    return result


def audit_fold(key, scope, genes, truth, prediction, identities, conditions, scores, bandwidth):
    """Rescore the selected Δ matrices and account for every condition and eligibility mask."""
    disease = scope == "disease"
    columns = (*CONDITION, COVERAGE) if disease else CONDITION
    observed = defaultdict(list)
    for condition_id, measured, estimated in zip(identities, truth, prediction, strict=True):
        row = conditions.pop((*key, scope, condition_id))
        recalculated = condition_scores(measured, estimated, genes)
        if disease:
            both_positive = (measured > 0) & (estimated > 0)
            both_negative = (measured < 0) & (estimated < 0)
            recalculated[COVERAGE] = 100 * int(np.count_nonzero(both_positive | both_negative)) / DISEASE_SIZE
        elif value(row, COVERAGE) is not None:
            raise ValueError("disease coverage appeared in the full-axis condition table")
        for metric in columns:
            number = recalculated[metric]
            expected = float(number) if math.isfinite(number) else None
            if not same(value(row, metric), expected):
                raise ValueError(f"incorrect held-out condition score: {(*key, scope, condition_id, metric)}")
            if expected is not None:
                observed[metric].append(expected)

    distances = dict(zip(DISTRIBUTION, distribution_scores(truth, prediction, bandwidth), strict=True))
    for metric in (*columns, *DISTRIBUTION):
        record = scores[*key, scope, metric]
        if metric in distances:
            expected, eligible = distances[metric], len(identities)
        else:
            eligible = len(observed[metric])
            expected = math.fsum(observed[metric]) / eligible if eligible else None
        if (not same(value(record), expected) or int(record["n_conditions"]) != len(identities)
                or int(record["n_eligible"]) != eligible or int(record["n_undefined"]) != len(identities) - eligible
                or int(record["overlap"]) != len(genes)):
            raise ValueError(f"incorrect independently rescored replicate: {(*key, scope, metric)}")


def audit_report(report, runs, step03, step04, disease_release, step08, *, bounded=False):
    """Check the report against immutable fits, ordered axes and independent matched controls."""
    expected_roster = coordinates(bounded)
    scores = read_table(report / "fold-scores.csv", SCORE_FIELDS)
    summary = read_table(report / "summary.csv", ("axis", "scope", "metric"))
    conditions = read_table(report / "condition-scores.csv", (*IDENTITY, "scope", "condition_id"))
    ledger = read_table(report / "runs.csv", IDENTITY)
    keys = {(*key, scope, metric) for key in expected_roster for scope in SCOPES
            for metric in ((*METRICS, COVERAGE) if scope == "disease" else METRICS)}
    if set(runs) != set(expected_roster) or set(ledger) != set(expected_roster) or set(scores) != keys:
        raise ValueError("independent scoring requires the complete coordinate and metric roster")
    memberships, source_hashes = disease_gene_sets(disease_release)
    union = {int(gene) for group in memberships.values() for gene in group}
    disease = axis_run.read_json(report / "disease-union.json")
    if (len(union) != DISEASE_SIZE or set(memberships) != {"nHF", "ICM", "genetic_cardiomyopathy", "TOX"}
            or disease != {"gene_indices": sorted(union), "memberships_sha256": source_hashes}):
        raise ValueError("formal report does not use the complete frozen disease union")
    bandwidths = audit_bandwidths(report, step03, step04, union)
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
            included = np.asarray([gene in union for gene in gene_indices])
            rows = np.flatnonzero(release.treatments[f"cold_drug_{fold}"].to_numpy() == "test")
            if bounded:
                rows = rows[:replay["test_conditions"]]
            ids = ordered_condition_ids(release, rows)
            linked = release.treatments.iloc[rows].control_row.to_numpy()
            controls = release.read_controls(0, len(release.controls))[linked]
            truth = treatment_rows(release, rows, "treatment") - controls
            prediction = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False) - controls
            if (np.count_nonzero(included) != DISEASE_OVERLAP[axis] or len(ids) != replay["test_conditions"]
                    or tuple(read_table(path / "conditions.csv", ("condition_id",))) != tuple((name,) for name in ids)
                    or not np.array_equal(truth, np.load(path / "measured_delta.npy", mmap_mode="r", allow_pickle=False))
                    or not np.array_equal(prediction, np.load(path / "predicted_delta.npy", mmap_mode="r", allow_pickle=False))):
                raise ValueError(f"independent held-out Delta or disease overlap differs: {key}")
        audit_fold(key, "axis", gene_indices, truth, prediction, ids, conditions, scores, bandwidths[fold, "axis"])
        audit_fold(key, "disease", gene_indices[included], truth[:, included], prediction[:, included],
                   ids, conditions, scores, bandwidths[fold, "disease"])
        if origin == "Step 08":
            for metric in METRICS:
                original = original_scores[ARM, fold, seed, metric]
                current = scores[*key, "axis", metric]
                if (value(original) != value(current)
                        or any(original[field] != current[field] for field in ("n_conditions", "n_eligible", "n_undefined"))):
                    raise ValueError(f"reused full-axis scores differ from Step 08: {(*key, metric)}")
    if conditions:
        raise ValueError("unexpected held-out conditions in the formal report")
    audit_aggregates(scores, summary)
    return {"status": "PASS", "coordinates": len(expected_roster), "score_rows": len(scores),
            "scoped_condition_rows": initial_conditions, "summary_rows": len(summary), "fold_scope_bandwidths": len(bandwidths)}
