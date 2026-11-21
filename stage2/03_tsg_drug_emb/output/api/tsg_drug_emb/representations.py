"""Read the seven parent-aligned representation matrices for frozen tenv folds."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


WIDTHS = {
    "ECFP4C": 2048, "ECFP6C": 2048, "FCFP4C": 2048, "FCFP6C": 2048,
    "E3FP": 1024, "KPGT": 2304, "UniMol": 512,
}
ARMS = tuple(WIDTHS)
SCHEMA = "tsg-drug-emb-representations-v1"
PARENT_COLUMNS = ["drug_parent_inchikey", "drug_parent_smiles"]


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def checked_vectors(arm, values, rows):
    if (values.shape != (rows, WIDTHS[arm]) or values.dtype != np.float32
            or not np.isfinite(values).all()):
        raise ValueError(f"{arm}: expected {rows} finite float32 vectors of width {WIDTHS[arm]}")
    if arm not in ("KPGT", "UniMol") and not ((values == 0) | (values == 1)).all():
        raise ValueError(f"{arm}: fingerprint vectors must be binary")
    return values


def fold_coverage(inputs):
    coverage = []
    for fold in range(3):
        partitions = {}
        for partition in ("train", "validation", "test"):
            rows = inputs.rows("cold-drug", fold, partition)
            keys = set(inputs.obs.iloc[rows].drug_parent_inchikey)
            partitions[partition] = keys
            coverage.append(dict(fold=fold, partition=partition, rows=len(rows),
                                 parents=len(keys), covered=len(keys)))
        if partitions["test"] & (partitions["train"] | partitions["validation"]):
            raise ValueError(f"fold {fold}: held-out parents leak into fit partitions")
        if set.union(*partitions.values()) != set(inputs.parent_lookup):
            raise ValueError(f"fold {fold}: representation roster does not cover accepted parents")
    return coverage


class RepresentationInputs:
    def __init__(self, inputs, directory):
        directory = Path(directory)
        self.manifest = json.loads((directory / "manifest.json").read_text())
        required = {"parents.csv", "e3fp_conformers.csv", *(f"{arm}.npy" for arm in ARMS)}
        if (inputs.source != "tenv" or self.manifest.get("schema") != SCHEMA
                or self.manifest.get("inputs_sha256") != sha256(inputs.directory / "manifest.json")
                or self.manifest.get("widths") != WIDTHS
                or self.manifest.get("coverage") != fold_coverage(inputs)
                or set(self.manifest.get("files", ())) != required):
            raise ValueError("representation freeze does not match accepted tenv inputs")
        for name, digest in self.manifest["files"].items():
            if sha256(directory / name) != digest:
                raise ValueError(f"representation freeze digest mismatch: {name}")
        parents = pd.read_csv(directory / "parents.csv")
        if not parents.equals(inputs.parents[PARENT_COLUMNS].reset_index(drop=True)):
            raise ValueError("representation freeze parent order or structures differ")
        conformers = pd.read_csv(directory / "e3fp_conformers.csv")
        if (not conformers.drug_parent_inchikey.equals(parents.drug_parent_inchikey)
                or not pd.api.types.is_integer_dtype(conformers.conformers)
                or not conformers.conformers.between(1, 3).all()):
            raise ValueError("invalid E3FP conformer ledger")
        self._values = {
            arm: checked_vectors(arm, np.load(directory / f"{arm}.npy", mmap_mode="r"), len(parents))
            for arm in ARMS
        }
        if not np.array_equal(self._values["ECFP4C"], inputs.morgan):
            raise ValueError("ECFP4C differs from accepted historical ECFP4")
        self._row_to_parent = inputs.obs.drug_parent_inchikey.map(inputs.parent_lookup).to_numpy()

    def width(self, arm):
        return self._values[arm].shape[1]

    def vectors(self, arm, rows):
        return self._values[arm][self._row_to_parent[rows]]
