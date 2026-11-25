"""Bounded reads of the plate-pooled Stage 2 tenv training releases."""
from __future__ import annotations

from contextlib import ExitStack
from hashlib import file_digest
import json
import operator
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

DEFAULT = Path(__file__).resolve().parents[1]
MAX_ROWS = 256
STRIPS = {"excluded-held-parent": "held_parent", "excluded-held-dose": "held_dose"}
TENV_BASELINE_AXIS = "hvg-5k"
TENV_AXIS_STEPS = {
    "hvg-3k": "03.3_tenv_hvg-3k",
    "hvg-4k": "03.4_tenv_hvg-4k",
    TENV_BASELINE_AXIS: "03.1_tenv_hvg-5k",
    "denoise-3k": "03.5_tenv_denoise-3k",
    "denoise-4k": "03.6_tenv_denoise-4k",
    "denoise-5k": "03.7_tenv_denoise-5k",
}
TENV_DEFAULT = DEFAULT.parent / TENV_AXIS_STEPS[TENV_BASELINE_AXIS] / "output"
TENV_SCHEMA = "stage2-tenv-plate-matched-training-v2"
TENV_FOLD_COLUMNS = tuple(f"{prefix}_{fold}"
                          for prefix in ("cold_drug", "cold_dose_drug") for fold in range(3))
TENV_STRIP_COLUMNS = tuple(f"{column}_reason" for column in TENV_FOLD_COLUMNS[3:])
TENV_MATRIX_FILES = {
    "C": "tenv-ref-vector.h5", "M_C": "plate-ref-mtx.h5",
    "M_T/treatment": "pert-mtx.h5", "M_T/delta": "pert-mtx.h5",
}
TENV_FILES = {"genes.csv", "controls.csv", "treatments.csv", "members.csv", "exclusions.csv"} | set(TENV_MATRIX_FILES.values())


def sha256(path):
    with path.open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def _read(path):
    return json.loads(path.read_bytes())


def _csv(path, **kwargs):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip", **kwargs)


class TenvTrainingData:
    """An audited axis release containing one control row per tenv plate."""

    control_key = ("plate",)
    fold_columns = TENV_FOLD_COLUMNS
    strip_columns = TENV_STRIP_COLUMNS
    strip_label = ""

    def __init__(self, root):
        self.root = root
        path = root / "manifest.json"
        self.manifest_sha256 = self.release_sha256 = sha256(path)
        self.manifest = _read(path)
        if (self.manifest.get("schema") != TENV_SCHEMA or self.manifest.get("source") != "tenv"
                or set(self.manifest.get("outputs", ())) != TENV_FILES
                or {path.name for path in root.iterdir()} != TENV_FILES | {"manifest.json"}):
            raise ValueError("tenv training dataset inventory differs")
        for name, digest in self.manifest["outputs"].items():
            if sha256(root / name) != digest:
                raise ValueError(f"tenv training payload digest differs: {name}")

        self.genes = _csv(self.root / "genes.csv")
        identity_types = dict.fromkeys((*self.control_key, "drug_parent"), str)
        self.controls = _csv(self.root / "controls.csv", dtype=identity_types)
        self.treatments = _csv(self.root / "treatments.csv", dtype=identity_types)
        self._handles = ExitStack()
        try:
            handles = {file: self._handles.enter_context(h5py.File(self.root / file, "r"))
                       for file in dict.fromkeys(TENV_MATRIX_FILES.values())}
            self._arrays = {key: handles[file][key] for key, file in TENV_MATRIX_FILES.items()}
            self._check_shapes()
        except BaseException:
            self.close()
            raise

    def _check_shapes(self):
        width = len(self.genes)
        n_controls, n_treatments = len(self.controls), len(self.treatments)
        control_key = list(self.control_key)
        treatment_key = [*control_key, "drug_parent", "dose_um"]
        if (self.genes.gene_index.duplicated().any() or width != self.manifest["axis_size"]
                or self.controls.duplicated(control_key).any()
                or self.treatments.duplicated(treatment_key).any()):
            raise ValueError("training matrix identities differ from the manifest")
        if (n_controls, n_treatments) != (
                self.manifest["counts"]["controls"], self.manifest["counts"]["treatments"]):
            raise ValueError("training matrix populations differ from the manifest")
        positions = self.treatments.control_row.to_numpy()
        if not np.issubdtype(positions.dtype, np.integer) or ((positions < 0) | (positions >= n_controls)).any():
            raise ValueError("invalid matched control row")
        matches = self.controls.iloc[positions][control_key].to_numpy()
        if not np.array_equal(matches, self.treatments[control_key].to_numpy()):
            raise ValueError("treatment points to a control from another plate or cell line")
        for column in self.fold_columns:
            allowed = {"train", "valid", "test"}
            if column in self.fold_columns[3:]:
                allowed.add(self.strip_label)
            if not self.treatments[column].isin(allowed).all():
                raise ValueError(f"invalid partition in {column}")
        for column, strip in zip(self.fold_columns[3:], self.strip_columns, strict=True):
            labels, reasons = self.treatments[column], self.treatments[strip]
            if not ((labels.eq(self.strip_label) & reasons.isin(STRIPS.values()))
                    | (labels.ne(self.strip_label) & reasons.eq(""))).all():
                raise ValueError(f"invalid exclusion reason in {strip}")
        for name, shape in {
            "C": (width,), "M_C": (n_controls, width),
            "M_T/treatment": (n_treatments, width), "M_T/delta": (n_treatments, width),
        }.items():
            dataset = self._arrays[name]
            if dataset.shape != shape or dataset.dtype != np.float32:
                raise ValueError(f"unsupported training matrix shape or dtype: {name}")
        for matrix, table in (("M_C", self.controls), ("M_T/treatment", self.treatments)):
            handle = self._arrays[matrix].file
            for column in table:
                values = handle[column]
                stored = values.asstr()[:] if h5py.check_string_dtype(values.dtype) else values[:]
                if not np.array_equal(stored, table[column].to_numpy()):
                    raise ValueError(f"tenv matrix row metadata differs: {matrix}/{column}")

    def control(self):
        """Return the equal-weight dataset-wide control C."""
        return self._arrays["C"][:]

    def read_controls(self, start, stop):
        """Read consecutive rows of M_C."""
        return self._read_window("M_C", start, stop)

    def read_treatments(self, start, stop, layer="treatment"):
        """Read consecutive rows of one M_T expression layer."""
        if layer not in {"treatment", "delta"}:
            raise ValueError(f"unknown treatment layer: {layer}")
        return self._read_window(f"M_T/{layer}", start, stop)

    def _read_window(self, name, start, stop):
        matrix = self._arrays[name]
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
        self._handles.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def open_tenv_training(root=None):
    """Verify and open the split-file tenv contract on a published prediction axis."""
    return TenvTrainingData(TENV_DEFAULT if root is None else Path(root))


def tenv_axis_output(axis_name, root=None):
    """Resolve a named tenv release under a Step 03 root."""
    if axis_name not in TENV_AXIS_STEPS:
        raise ValueError(f"unknown tenv training axis: {axis_name}")
    step = DEFAULT.parent if root is None else Path(root)
    return step / TENV_AXIS_STEPS[axis_name] / "output"


def open_tenv_axis_training(axis_name, root=None):
    """Verify and open one named tenv axis under a Step 03 root."""
    return open_tenv_training(tenv_axis_output(axis_name, root))
