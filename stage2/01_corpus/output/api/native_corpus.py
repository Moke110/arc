"""Verified source-native integer counts and observation-level control matches."""
from __future__ import annotations

import json
import operator
import re
from collections import defaultdict
from hashlib import file_digest
from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem
from scipy.sparse import csr_matrix

SCHEMA = "stage2-native-corpus-v1"
MAX_ROWS = 256
SOURCE_NAME = re.compile(r"[A-Za-z0-9_-]+")


class CorpusError(ValueError):
    """A source or release violates its data contract."""


def sha256(path):
    with Path(path).open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def _payload(root, source):
    if not SOURCE_NAME.fullmatch(source):
        raise CorpusError(f"invalid source name: {source!r}")
    return root / f"{source}.h5ad"


class Source:
    """A native gene index and a bounded reader for an immutable source payload."""

    def __init__(self, path, index, common_size, observation_path=None):
        self._file = h5py.File(path, "r")
        try:
            self.obs = read_elem(self._file["obs"])
            if observation_path is not None:
                with h5py.File(observation_path, "r") as metadata:
                    overlay = read_elem(metadata["obs"])
                if not self.obs.index.equals(overlay.index):
                    raise CorpusError("observation metadata is ordered differently from counts")
                self.obs[overlay.columns] = overlay
            self.var = read_elem(self._file["var"])
            self._matrix = self._file["X"]
            self.shape = (len(self.obs), len(self.var))
            self._sparse = isinstance(self._matrix, h5py.Group)
            self._common_size = common_size
            self._index = {int(gene): column for gene, column in index.items()}
            self._check_layout()
            self._set_controls()
        except BaseException:
            self.close()
            raise

    def _check_layout(self):
        if self._sparse:
            if self._matrix.attrs["encoding-type"] != "csr_matrix":
                raise CorpusError("sparse counts must be CSR")
            shape = tuple(self._matrix.attrs["shape"])
            dtype = self._matrix["data"].dtype
        else:
            shape, dtype = self._matrix.shape, self._matrix.dtype
        if shape != self.shape or dtype.kind not in "iu":
            raise CorpusError("raw count shape or integer dtype violates the corpus contract")
        for frame in (self.obs, self.var):
            if frame.index.has_duplicates or frame.index.hasnans or "" in frame.index:
                raise CorpusError("missing or duplicate source-native identities")
        if len(set(self._index.values())) != len(self._index) or any(
            not 0 <= gene < self._common_size or not 0 <= column < self.shape[1]
            for gene, column in self._index.items()
        ):
            raise CorpusError("common gene coverage map is out of bounds or ambiguous")
        if len(self._file.get("layers", {})):
            raise CorpusError("count corpora cannot carry expression layers")

    def _set_controls(self):
        roles = self.obs["role"]
        groups = self.obs["matched_control_group"]
        if not roles.isin(["treated", "vehicle", "blank", "unknown"]).all() or groups.isna().any():
            raise CorpusError("invalid roles or missing control-group metadata")
        self._controls = defaultdict(list)
        for index, (role, group) in enumerate(zip(roles, groups, strict=True)):
            if role == "vehicle" and group:
                self._controls[group].append(index)
        self.obs["match_status"] = [
            ("matched" if self._controls.get(group) else "unmatched") if role == "treated" else role
            for role, group in zip(roles, groups, strict=True)
        ]

    @property
    def shared_axis_complete(self):
        return len(self._index) == self._common_size

    def native_columns(self, common_gene_indices):
        """Resolve measured common genes in caller order, failing on unmeasured genes."""
        columns = []
        for index in common_gene_indices:
            index = operator.index(index)
            if index not in self._index:
                raise CorpusError(f"common gene {index} is not measured by this source")
            columns.append(self._index[index])
        return tuple(columns)

    def matched_control_rows(self, index):
        """Resolve a treated observation to vehicle row positions, if any."""
        index = operator.index(index)
        if not 0 <= index < self.shape[0]:
            raise CorpusError("invalid observation row")
        role = self.obs["role"].iloc[index]
        if role != "treated":
            return ()
        group = self.obs["matched_control_group"].iloc[index]
        return tuple(self._controls.get(group, ()))

    def counts(self, start, stop):
        start, stop = operator.index(start), operator.index(stop)
        if not 0 <= start <= stop <= self.shape[0] or stop - start > MAX_ROWS:
            raise CorpusError(f"count window must contain at most {MAX_ROWS} valid rows")
        if not self._sparse:
            block = self._matrix[start:stop]
        else:
            matrix = self._matrix
            pointers = matrix["indptr"][start:stop + 1]
            lo, hi = pointers[0], pointers[-1]
            block = csr_matrix(
                (matrix["data"][lo:hi], matrix["indices"][lo:hi], pointers - lo),
                shape=(stop - start, self.shape[1]),
            ).toarray()
        if np.any(block < 0):
            raise CorpusError("source counts must be nonnegative")
        return block

    def batches(self, batch_size=MAX_ROWS):
        batch_size = operator.index(batch_size)
        if not 1 <= batch_size <= MAX_ROWS:
            raise CorpusError(f"batch size must be within 1 to {MAX_ROWS}")
        for start in range(0, self.shape[0], batch_size):
            yield start, self.counts(start, min(start + batch_size, self.shape[0]))

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class Corpus:
    """A release of independent native gene spaces with one common-axis identity."""

    def __init__(self, root):
        self.root = Path(root)
        self.manifest = json.loads((self.root / "manifest.json").read_text())
        if self.manifest["schema"] != SCHEMA:
            raise CorpusError("unsupported corpus release version")
        self.manifest_sha256 = sha256(self.root / "manifest.json")
        self._sources = {}

    def source(self, name):
        if name not in self._sources:
            record = self.manifest["sources"][name]
            path = _payload(self.root, name)
            if sha256(path) != record["sha256"]:
                raise CorpusError(f"source payload digest differs: {name}")
            observation_path = None
            if "observation_sha256" in record:
                observation_path = path.with_suffix(".obs.h5")
                if sha256(observation_path) != record["observation_sha256"]:
                    raise CorpusError(f"source observation digest differs: {name}")
            source = Source(path, record["common_columns"], self.manifest["common_axis"]["size"],
                            observation_path)
            if list(source.shape) != record["shape"] or (
                source.shared_axis_complete != record["shared_axis_complete"]
            ):
                source.close()
                raise CorpusError(f"source coverage or population differs: {name}")
            self._sources[name] = source
        return self._sources[name]

    def close(self):
        for source in self._sources.values():
            source.close()
        self._sources.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def open_corpus(root):
    return Corpus(root)


def seal_corpus(root, mappings, common_size, axis_sha256, inputs, exclusions, *, artifacts=()):
    """Measure and seal payloads already staged by source-specific adapters."""
    root = Path(root)
    manifest = {
        "schema": SCHEMA,
        "common_axis": {"size": common_size, "sha256": axis_sha256},
        "inputs": inputs,
        "exclusions": exclusions,
        "sources": {},
        "outputs": {name: sha256(root / name) for name in artifacts},
    }
    for name, mapping in sorted(mappings.items()):
        path = _payload(root, name)
        observation_path = path.with_suffix(".obs.h5")
        if not observation_path.exists():
            observation_path = None
        with Source(path, mapping, common_size, observation_path) as source:
            for _ in source.batches():
                pass
            manifest["sources"][name] = {
                "sha256": sha256(path),
                "shape": list(source.shape),
                "common_columns": {str(index): column for index, column in mapping.items()},
                "shared_axis_complete": source.shared_axis_complete,
                "match_status_counts": source.obs.match_status.value_counts().to_dict(),
            }
        manifest["outputs"][path.name] = manifest["sources"][name]["sha256"]
        if observation_path is not None:
            digest = sha256(observation_path)
            manifest["sources"][name]["observation_sha256"] = digest
            manifest["outputs"][observation_path.name] = digest
    path = root / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest
