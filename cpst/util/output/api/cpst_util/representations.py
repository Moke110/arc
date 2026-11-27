"""Provenance-bound desalted-parent identities and drug representations."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from ._paths import repository

WIDTHS = {
    "ECFP4C": 2048, "ECFP6C": 2048, "FCFP4C": 2048, "FCFP6C": 2048,
    "E3FP": 1024, "KPGT": 2304, "UniMol": 512,
}
ARMS = tuple(WIDTHS)
SCHEMA = "tsg-drug-emb-representations-v1"
PARENT_COLUMNS = ["drug_parent_inchikey", "drug_parent_smiles"]

# Exact manifest bytes of the accepted ADR 0079 representation releases.
ACCEPTED = {
    "tenv": ("cpst/util/output/representations",
             "bc10c6af897b1b0d939a685aa198e03b22a3392e3891b721f70cc49f1f028678"),
    "tahoe_bulk": ("cpst/util/output/representations/tahoe",
                   "4d00b23d6dcf0ddfbba4f60725b7f32701103a475f278e1db94eeebb968fbeac"),
}


def ordered_parents(observations):
    """Extract a unique parent table in InChIKey order from repeated observations."""
    parents = observations[PARENT_COLUMNS].drop_duplicates()
    if (parents.isna().any().any() or parents.eq("").any().any()
            or parents[PARENT_COLUMNS[0]].duplicated().any()):
        raise ValueError("parent identities contain missing or conflicting structures")
    return parents.sort_values(PARENT_COLUMNS[0], ignore_index=True)


def _sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _checked_vectors(arm, values, rows):
    if (values.shape != (rows, WIDTHS[arm]) or values.dtype != np.float32
            or not np.isfinite(values).all()):
        raise ValueError(f"{arm}: expected {rows} finite float32 vectors of width {WIDTHS[arm]}")
    if arm not in ("KPGT", "UniMol") and not ((values == 0) | (values == 1)).all():
        raise ValueError(f"{arm}: fingerprint vectors must be binary")
    return values


class RepresentationTable:
    """Open one accepted source with vectors indexed by desalted parent InChIKey."""

    def __init__(self, directory=None, *, source, manifest_sha256=None):
        release_path, accepted_digest = ACCEPTED[source]
        directory = (repository() / release_path) if directory is None else Path(directory)
        payload = (directory / "manifest.json").read_bytes()
        self.manifest_sha256 = hashlib.sha256(payload).hexdigest()
        if self.manifest_sha256 != (accepted_digest if manifest_sha256 is None else manifest_sha256):
            raise ValueError("representation manifest digest mismatch")
        self.manifest = json.loads(payload)
        files = self.manifest.get("files", {})
        required = {"parents.csv", "e3fp_conformers.csv", *(f"{arm}.npy" for arm in ARMS)}
        if (self.manifest.get("source", "tenv") != source
                or self.manifest.get("schema") != SCHEMA
                or self.manifest.get("widths") != WIDTHS
                or not re.fullmatch(r"[0-9a-f]{64}", self.manifest.get("inputs_sha256", ""))
                or set(files) != required):
            raise ValueError("representation manifest does not match the source contract")
        for name, expected in files.items():
            if _sha256(directory / name) != expected:
                raise ValueError(f"representation freeze digest mismatch: {name}")

        parents = pd.read_csv(directory / "parents.csv", dtype=str, keep_default_na=False)
        if (list(parents.columns) != PARENT_COLUMNS or parents.empty
                or parents.eq("").any().any() or parents[PARENT_COLUMNS[0]].duplicated().any()):
            raise ValueError("invalid drug parent identity table")
        conformers = pd.read_csv(directory / "e3fp_conformers.csv")
        if (not conformers.drug_parent_inchikey.equals(parents.drug_parent_inchikey)
                or not pd.api.types.is_integer_dtype(conformers.conformers)
                or not conformers.conformers.between(1, 3).all()):
            raise ValueError("invalid E3FP conformer ledger")
        self._parents = parents
        self._rows = {parent: i for i, parent in enumerate(parents.drug_parent_inchikey)}
        self._values = {
            arm: _checked_vectors(arm, np.load(directory / f"{arm}.npy", mmap_mode="r"), len(parents))
            for arm in ARMS
        }

    @property
    def parents(self):
        """Return the ordered parent identities and structures."""
        return self._parents.copy()

    def verified_parent_keys(self, parents):
        """Check structures against accepted identities and return keys in input order."""
        keys = parents[PARENT_COLUMNS[0]].tolist()
        accepted = self._parents.set_index(PARENT_COLUMNS[0])[PARENT_COLUMNS[1]]
        if any(key not in self._rows for key in keys):
            raise ValueError("parent identity missing from the shared representation release")
        structures = parents[PARENT_COLUMNS[1]]
        if structures.isna().any() or accepted.reindex(keys).tolist() != structures.tolist():
            raise ValueError("parent structures differ from the shared representation release")
        return keys

    def width(self, arm):
        return self._values[arm].shape[1]

    def vectors(self, arm, parent_inchikeys):
        """Return vectors in caller key order, including repeated parents."""
        rows = np.asarray([self._rows[key] for key in parent_inchikeys], dtype=np.intp)
        return self._values[arm][rows]
