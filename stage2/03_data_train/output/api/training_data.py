"""Bounded reads of the plate-matched Stage 2 training release."""
from __future__ import annotations

from hashlib import file_digest
import json
import operator
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

DEFAULT = Path(__file__).resolve().parents[1]
SCHEMA = "stage2-plate-matched-training-v1"
DATASET_SCHEMA = "stage2-plate-matched-training-dataset-v1"
MAX_ROWS = 256
SCENARIOS = ("cold-drug", "cold-drug-dose")
FOLD_COLUMNS = tuple(f"{scenario.replace('-', '_')}_{fold}" for scenario in SCENARIOS for fold in range(3))
STRIP_COLUMNS = tuple(f"{column}_reason" for column in FOLD_COLUMNS[3:])
STRIPS = {"excluded-held-parent": "held_parent", "excluded-held-dose": "held_dose"}
FILES = {"training.h5", "genes.csv", "controls.csv", "treatments.csv", "members.csv", "exclusions.csv"}


def sha256(path):
    with path.open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def _read(path):
    return json.loads(path.read_bytes())


def _csv(path, **kwargs):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip", **kwargs)


class TrainingData:
    """One dataset with C, M_C and M_T in the frozen prediction-axis order."""

    def __init__(self, root, source):
        path = root / "manifest.json"
        pair = _read(path)
        self.release_sha256 = sha256(path)
        sources = pair.get("sources", {})
        if pair.get("schema") != SCHEMA or source not in sources:
            raise ValueError(f"source is missing from the training release: {source}")
        if not source or Path(source).name != source or source in {".", ".."}:
            raise ValueError(f"invalid training source: {source}")
        self.root = root / source
        self.manifest_sha256 = sha256(self.root / "manifest.json")
        if self.manifest_sha256 != sources[source]:
            raise ValueError(f"training dataset manifest digest differs: {source}")
        self.manifest = _read(self.root / "manifest.json")
        if (self.manifest.get("schema") != DATASET_SCHEMA or self.manifest.get("source") != source
                or self.manifest.get("inputs", {}).get("raw_sha256") != pair.get("raw_sha256")
                or self.manifest.get("inputs", {}).get("axis_sha256") != pair.get("axis_sha256")
                or set(self.manifest.get("outputs", ())) != FILES
                or {path.name for path in self.root.iterdir()} != FILES | {"manifest.json"}):
            raise ValueError(f"training dataset inventory or generation differs: {source}")
        for name, digest in self.manifest["outputs"].items():
            if sha256(self.root / name) != digest:
                raise ValueError(f"training payload digest differs: {source}/{name}")
        self.genes = _csv(self.root / "genes.csv")
        identity_types = {"plate": str, "cell_line": str}
        self.controls = _csv(self.root / "controls.csv", dtype=identity_types)
        self.treatments = _csv(self.root / "treatments.csv", dtype=identity_types)
        self._handle = h5py.File(self.root / "training.h5", "r")
        try:
            self._check_shapes()
        except BaseException:
            self.close()
            raise

    def _check_shapes(self):
        width = len(self.genes)
        n_controls, n_treatments = len(self.controls), len(self.treatments)
        if (self.genes.gene_index.duplicated().any() or width != self.manifest["axis_size"]
                or self.controls.duplicated(["plate", "cell_line"]).any()
                or self.treatments.duplicated(["plate", "cell_line", "drug_parent", "dose_um"]).any()):
            raise ValueError("training matrix identities differ from the manifest")
        if (n_controls, n_treatments) != (
                self.manifest["counts"]["controls"], self.manifest["counts"]["treatments"]):
            raise ValueError("training matrix populations differ from the manifest")
        positions = self.treatments.control_row.to_numpy()
        if not np.issubdtype(positions.dtype, np.integer) or ((positions < 0) | (positions >= n_controls)).any():
            raise ValueError("invalid matched control row")
        matches = self.controls.iloc[positions][["plate", "cell_line"]].to_numpy()
        if not np.array_equal(matches, self.treatments[["plate", "cell_line"]].to_numpy()):
            raise ValueError("treatment points to a control from another plate or cell line")
        for column in FOLD_COLUMNS:
            allowed = {"train", "valid", "test"} if column not in FOLD_COLUMNS[3:] else {
                "train", "valid", "test", "excluded"
            }
            if not self.treatments[column].isin(allowed).all():
                raise ValueError(f"invalid partition in {column}")
        for column, strip in zip(FOLD_COLUMNS[3:], STRIP_COLUMNS, strict=True):
            labels, reasons = self.treatments[column], self.treatments[strip]
            if not ((labels.eq("excluded") & reasons.isin(STRIPS.values()))
                    | (labels.ne("excluded") & reasons.eq(""))).all():
                raise ValueError(f"invalid exclusion reason in {strip}")
        for name, shape in {
            "C": (width,), "M_C": (n_controls, width),
            "M_T/treatment": (n_treatments, width), "M_T/delta": (n_treatments, width),
        }.items():
            dataset = self._handle[name]
            if dataset.shape != shape or dataset.dtype != np.float32:
                raise ValueError(f"unsupported training matrix shape or dtype: {name}")

    def control(self):
        """Return the equal-weight dataset-wide control C."""
        return self._handle["C"][:]

    def read_controls(self, start, stop):
        """Read consecutive rows of M_C."""
        return self._read_window("M_C", start, stop)

    def read_treatments(self, start, stop, layer="treatment"):
        """Read consecutive rows of one M_T expression layer."""
        if layer not in {"treatment", "delta"}:
            raise ValueError(f"unknown treatment layer: {layer}")
        return self._read_window(f"M_T/{layer}", start, stop)

    def _read_window(self, name, start, stop):
        matrix = self._handle[name]
        start, stop = operator.index(start), operator.index(stop)
        if not 0 <= start <= stop <= matrix.shape[0] or stop - start > MAX_ROWS:
            raise ValueError(f"training window must contain at most {MAX_ROWS} valid rows")
        return matrix[start:stop]

    def membership(self):
        """Read source profile membership, distinguishing control and treatment rows."""
        return _csv(self.root / "members.csv", dtype={"profile_id": str})

    def exclusions(self):
        """Read the source profile identities and reasons omitted from this release."""
        return _csv(self.root / "exclusions.csv", dtype={"profile_id": str})

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def open_training(root=None, *, source):
    """Open an accepted training dataset and verify its manifest and payloads."""
    return TrainingData(DEFAULT if root is None else Path(root), source)
