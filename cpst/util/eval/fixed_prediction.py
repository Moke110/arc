"""Score the frozen condition population and ordered prediction axis."""

from hashlib import sha256 as hash_bytes
import math

import numpy as np
import pandas as pd

from cpst_util.release_names import original_axis_manifest_sha256
from prediction_axis import axis_coordinates, open_axes
from tsg_axis.inputs import score_fixed_prediction
from .metrics import (CONDITION, DISTRIBUTION, METRICS, condition_scores, distribution_scores,
                      fold_bandwidth, summarize_conditions)

PARTITIONS = {"train": 64, "valid": 24, "test": 24}


def read_table(path):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip")


def verify_inputs(path, config, data, base, axis, fold, smoke):
    family, width = axis_coordinates(axis)
    axes = open_axes(family=family)
    released_axis = axes.axis(width)
    stage03 = data.release.manifest["inputs"]
    baseline = base.release.manifest["inputs"]
    if (config["axis_sha256"] != stage03["axis_sha256"]
            or config["axis_sha256"] != released_axis.sha256
            or config["axis_family_manifest_sha256"] not in {
                axes.manifest_sha256, original_axis_manifest_sha256(axes.manifest_sha256)}
            or config["training_sha256"] != data.release.manifest_sha256
            or config["raw_sha256"] != stage03["raw_sha256"]
            or config["raw_sha256"] != baseline["raw_sha256"]
            or config["folds"] != stage03["folds"] or config["folds"] != baseline["folds"]
            or config.get("fold_release_sha256") != data.fold_release_sha256
            or config.get("fold_release_sha256") != base.fold_release_sha256
            or not read_table(path / "genes.csv").equals(data.genes)):
        raise ValueError(f"axis, training, raw or frozen fold identity differs: {path}")
    for phase, limit in PARTITIONS.items():
        rows = data.rows(fold, phase)
        reference_rows = base.rows(fold, phase)
        if smoke:
            rows, reference_rows = rows[:limit], reference_rows[:limit]
        ids = data.condition_ids(rows)
        if (rows.tolist() != config["partitions"][phase] or list(ids) != config["conditions"][phase]
                or ids != base.condition_ids(reference_rows)):
            raise ValueError(f"registered {phase} condition IDs differ from the shared fold: {path}")
    test = np.asarray(config["partitions"]["test"], dtype=np.int64)
    reference_test = base.rows(fold, "test")
    return test, reference_test[:len(test)]


def fixed_bandwidth(base, fold):
    data = base.release
    rows = base.rows(fold, "train")
    parts = []
    for start in range(0, len(data.treatments), 64):
        stop = min(start + 64, len(data.treatments))
        indexes = rows[(rows >= start) & (rows < stop)]
        if len(indexes):
            parts.append(data.read_treatments(start, stop, "delta")[indexes - start])
    measured = np.concatenate(parts)
    record = {"fold": fold, "bandwidth": fold_bandwidth(measured), "rows": len(measured),
              "width": len(data.genes), "pair_count": math.comb(len(measured), 2),
              "training_manifest_sha256": data.manifest_sha256,
              "training_delta_values_sha256": hash_bytes(measured.tobytes()).hexdigest(),
              "training_conditions_sha256": hash_bytes(base.condition_table(rows).to_csv(index=False).encode()).hexdigest()}
    if base.fold_release_sha256 is not None:
        record["fold_release_sha256"] = base.fold_release_sha256
    return record


def read_paired_prediction(path, config, manifest, data, base, axis, fold, smoke):
    rows, reference_rows = verify_inputs(path, config, data, base, axis, fold, smoke)
    conditions = read_table(path / "conditions.csv")
    if (not conditions.equals(data.condition_table(rows).reset_index(drop=True))
            or len(conditions) != manifest["test_conditions"]):
        raise ValueError(f"verified run has different held-out conditions: {path}")
    prediction = np.load(path / "prediction.npy", mmap_mode="r", allow_pickle=False)
    target = score_fixed_prediction(data, fold, prediction, conditions.condition_id, rows=rows)
    for name, expected in (("predicted_delta", target.predicted_delta),
                           ("measured_delta", target.measured_delta)):
        array = np.load(path / f"{name}.npy", mmap_mode="r", allow_pickle=False)
        if array.dtype != expected.dtype or not np.array_equal(array, expected):
            raise ValueError(f"fixed prediction {name} differs from paired source arrays: {path}")

    return conditions, reference_rows, target


def score_paired_prediction(conditions, paired, genes, bandwidth):
    by_condition = pd.DataFrame(
        condition_scores(truth, prediction, genes)
        for truth, prediction in zip(paired.measured_delta, paired.predicted_delta, strict=True)
    )
    by_condition.insert(0, "condition_id", conditions.condition_id)
    distribution = dict(zip(DISTRIBUTION, distribution_scores(paired.measured_delta, paired.predicted_delta,
                                                                bandwidth), strict=True))
    by_fold = pd.DataFrame([
        {"metric": metric, **(summarize_conditions(by_condition[metric]) if metric in CONDITION else
                             {"value": distribution[metric], "n_conditions": len(conditions),
                              "n_eligible": len(conditions), "n_undefined": 0})}
        for metric in METRICS
    ])
    return by_condition, by_fold
