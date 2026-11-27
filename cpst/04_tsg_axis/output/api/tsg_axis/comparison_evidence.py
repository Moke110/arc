"""Read signed Jaccard evidence at the verified matched-control prediction boundary."""

import numpy as np
import pandas as pd

from .inputs import score_fixed_prediction
from cpst_util.derivation.release import sha256
from cpst_util.metrics import signed_jaccard

EVIDENCE_FILES = ("conditions.csv", "genes.csv", "prediction.npy", "measured_delta.npy", "predicted_delta.npy")


def read_fixed_jaccard(run, config, manifest, data):
    """Return ordered held-out IDs and JAC-top200 for a verified fit."""
    for name in EVIDENCE_FILES:
        if sha256(run / name) != manifest["outputs"][name]:
            raise ValueError(f"sealed fit payload differs for {run / name}")
    genes = pd.read_csv(run / "genes.csv", float_precision="round_trip")
    if (config["training_sha256"] != data.release.manifest_sha256
            or config["axis_sha256"] != data.release.manifest["inputs"]["axis_sha256"]
            or not genes.equals(data.genes)):
        raise ValueError(f"registered training release or prediction axis differs for {run}")
    fold = manifest["fold"]
    rows = data.rows(fold, "test")
    conditions = pd.read_csv(run / "conditions.csv", float_precision="round_trip")
    if (rows.tolist() != config["partitions"]["test"]
            or conditions.condition_id.tolist() != config["conditions"]["test"]
            or len(conditions) != manifest["test_conditions"]
            or not conditions.equals(data.condition_table(rows).reset_index(drop=True))):
        raise ValueError(f"registered held-out conditions differ for {run}")

    prediction = np.load(run / "prediction.npy", mmap_mode="r", allow_pickle=False)
    targets = score_fixed_prediction(data, fold, prediction, conditions.condition_id)
    for name in ("measured_delta", "predicted_delta"):
        recorded = np.load(run / f"{name}.npy", mmap_mode="r", allow_pickle=False)
        expected = getattr(targets, name)
        if recorded.dtype != expected.dtype or not np.array_equal(recorded, expected):
            raise ValueError(f"matched-control {name} differs from sealed fit for {run}")
    gene_indices = genes.gene_index.to_numpy()
    scores = np.asarray([
        signed_jaccard(truth, predicted, gene_indices)[0]
        for truth, predicted in zip(targets.measured_delta, targets.predicted_delta, strict=True)
    ])
    return targets.condition_ids, scores
