"""Verified source-native integer counts and observation-level control matches."""
from __future__ import annotations

import json
import operator
import re
from collections import defaultdict
from hashlib import file_digest, sha256 as hash_bytes
from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem, write_elem
from scipy.sparse import csr_matrix

LEGACY_SCHEMA = "stage2-native-corpus-v1"
SCHEMA = "stage2-native-corpus-v2"
DATASET_SCHEMA = "stage2-native-dataset-v1"
DATASETS = frozenset({
    "tenv", "tahoe_bulk", "sctox", "dtoxs_ipsc_conv", "dtoxs_ipsc_dge", "dtoxs_primary",
    "gse247541",
})
MAX_ROWS = 256
SOURCE_NAME = re.compile(r"[A-Za-z0-9_-]+")


class CorpusError(ValueError):
    """A source or release violates its data contract."""


def sha256(path):
    with Path(path).open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def resolve_generation(manifest, current, requested, sources, historical_manifest):
    """Resolve the requested digest against the current or equivalent sealed corpus."""
    if manifest.get("schema") not in {SCHEMA, LEGACY_SCHEMA}:
        raise CorpusError("canonical native corpus provenance is missing")
    if requested == current:
        return current
    predecessor = manifest.get("predecessor", {})
    unchanged = predecessor.get("unchanged_sources", {})
    if (manifest["schema"] == SCHEMA and predecessor.get("manifest_sha256") == requested
            and all(unchanged.get(source) == source for source in sources)
            and historical_manifest.is_file()):
        payload = historical_manifest.read_bytes()
        if (hash_bytes(payload).hexdigest() == requested
                and json.loads(payload).get("schema") == LEGACY_SCHEMA):
            return requested
    raise CorpusError("requested corpus or axis does not match the canonical release inputs")

def _payload(root, source):
    if not SOURCE_NAME.fullmatch(source):
        raise CorpusError(f"invalid source name: {source!r}")
    return root / f"{source}.h5ad"


class Source:
    """A native gene index and a bounded reader for an immutable source payload."""

    def __init__(self, path, index, common_size, observation_path=None, *, sealed=False):
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
            self._set_controls(sealed)
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

    def _set_controls(self, sealed):
        roles = self.obs["role"]
        groups = self.obs["matched_control_group"]
        if not roles.isin(["treated", "vehicle", "blank", "unknown"]).all() or groups.isna().any():
            raise CorpusError("invalid roles or missing control-group metadata")
        if "control_contexts" in self.obs:
            self._contexts = [
                tuple(tuple(context) for context in json.loads(contexts))
                for contexts in self.obs.control_contexts
            ]
        else:
            self._contexts = [(group,) if group else () for group in groups]
        self._controls = defaultdict(list)
        for index, (role, contexts) in enumerate(zip(roles, self._contexts, strict=True)):
            if role == "vehicle":
                for context in contexts:
                    self._controls[context].append(index)
        status = [
            ("matched" if any(self._controls.get(context) for context in contexts) else "unmatched")
            if role == "treated" else role
            for role, contexts in zip(roles, self._contexts, strict=True)
        ]
        if sealed and ("match_status" not in self.obs or self.obs.match_status.tolist() != status):
            raise CorpusError("stored match status differs from observation control contexts")
        self.obs["match_status"] = status

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
        return tuple(sorted({
            row for context in self._contexts[index] for row in self._controls.get(context, ())
        }))

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
        if self.manifest["schema"] not in {SCHEMA, LEGACY_SCHEMA}:
            raise CorpusError("unsupported corpus release version")
        self.manifest_sha256 = sha256(self.root / "manifest.json")
        self._sources = {}

    def resolve_generation(self, requested, sources, *, historical_manifest=None):
        if historical_manifest is None:
            historical_manifest = self.root.resolve().parent / "historical/native-v1/manifest.json"
        return resolve_generation(self.manifest, self.manifest_sha256, requested, sources,
                                  Path(historical_manifest))

    def source(self, name):
        if name not in self._sources:
            record = self.manifest["sources"][name]
            path = _payload(self.root, name)
            if self.manifest["schema"] == SCHEMA:
                path = _payload(self.root / name, name)
                manifest_path = path.parent / "manifest.json"
                if sha256(manifest_path) != record["manifest_sha256"]:
                    raise CorpusError(f"dataset manifest digest differs: {name}")
                dataset = json.loads(manifest_path.read_text())
                expected = {key: value for key, value in record.items() if key != "manifest_sha256"}
                output_prefix = f"{name}/"
                expected_outputs = {
                    filename.removeprefix(output_prefix): digest
                    for filename, digest in self.manifest["outputs"].items()
                    if filename.startswith(output_prefix) and filename != f"{name}/manifest.json"
                }
                if (dataset["schema"], dataset["dataset"], dataset["source"]) != (
                    DATASET_SCHEMA, name, expected
                ) or dataset["outputs"] != expected_outputs or (
                    expected_outputs.get(path.name) != record["sha256"]
                ):
                    raise CorpusError(f"dataset manifest disagrees with release inventory: {name}")
            if sha256(path) != record["sha256"]:
                raise CorpusError(f"source payload digest differs: {name}")
            observation_path = None
            if "observation_sha256" in record:
                observation_path = path.with_suffix(".obs.h5")
                if sha256(observation_path) != record["observation_sha256"]:
                    raise CorpusError(f"source observation digest differs: {name}")
            source = Source(path, record["common_columns"], self.manifest["common_axis"]["size"],
                            observation_path, sealed=self.manifest["schema"] == SCHEMA)
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


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def seal_corpus(root, mappings, common_size, axis_sha256, inputs, exclusions, *,
                supplements=None, predecessor=None):
    """Publish and seal staged source payloads in their dataset directories."""
    root = Path(root)
    supplements = supplements or {}
    manifest = {
        "schema": SCHEMA,
        "common_axis": {"size": common_size, "sha256": axis_sha256},
        "inputs": inputs,
        "exclusions": exclusions,
        "sources": {},
        "outputs": {},
    }
    if predecessor is not None:
        manifest["predecessor"] = predecessor
    for name, mapping in sorted(mappings.items()):
        staged_path = _payload(root, name)
        dataset_root = root / name
        dataset_root.mkdir()
        path = staged_path.rename(_payload(dataset_root, name))
        observation_path = staged_path.with_suffix(".obs.h5")
        if not observation_path.exists():
            observation_path = None
        with Source(path, mapping, common_size, observation_path) as source:
            for _ in source.batches():
                pass
            observations = source.obs
            record = {
                "shape": list(source.shape),
                "common_columns": {str(index): column for index, column in mapping.items()},
                "shared_axis_complete": source.shared_axis_complete,
                "match_status_counts": observations.match_status.value_counts().to_dict(),
            }
        with h5py.File(path, "r+") as handle:
            del handle["obs"]
            write_elem(handle, "obs", observations)
        if observation_path is not None:
            observation_path.unlink()
        record["sha256"] = sha256(path)
        outputs = {path.name: record["sha256"]}
        for filename in supplements.get(name, ()):
            supplement = (root / filename).rename(dataset_root / filename)
            outputs[filename] = sha256(supplement)
        dataset_path = dataset_root / "manifest.json"
        _write_json(dataset_path, {
            "schema": DATASET_SCHEMA,
            "dataset": name,
            "source": record,
            "outputs": outputs,
        })
        manifest["sources"][name] = {**record, "manifest_sha256": sha256(dataset_path)}
        manifest["outputs"].update({
            f"{name}/{filename}": digest for filename, digest in outputs.items()
        })
        manifest["outputs"][f"{name}/manifest.json"] = manifest["sources"][name]["manifest_sha256"]
    _write_json(root / "manifest.json", manifest)
    return manifest
