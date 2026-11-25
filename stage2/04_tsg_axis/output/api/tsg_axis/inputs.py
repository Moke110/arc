"""Bounded Stage 03 reads and the identity-checked Step 04 prediction boundary."""
from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
import json
from operator import index
from pathlib import Path
import sys

import numpy as np

STAGE = Path(__file__).resolve().parents[4]
for relative in ("03_data_train/output/api", "util/output/api"):
    sys.path.insert(0, str(STAGE / relative))

from stage2_util.representations import RepresentationTable, WIDTHS
from training_data import open_tenv_training


class _ExpressionRows:
    """Gather arbitrary row order with a fixed-size cache of bounded reads."""

    block_rows = 64

    def __init__(self, release, layer):
        self.release = release
        self.layer = layer
        self._read = lru_cache(maxsize=24)(self._read_block)

    def _read_block(self, block):
        start = block * self.block_rows
        return self.release.read_treatments(start, min(start + self.block_rows, len(self.release.treatments)),
                                            self.layer)

    def __getitem__(self, rows):
        rows = np.asarray(rows, dtype=np.int64)
        if rows.ndim != 1 or np.any((rows < 0) | (rows >= len(self.release.treatments))):
            raise ValueError("invalid treatment row positions")
        result = np.empty((len(rows), len(self.release.genes)), dtype=np.float32)
        blocks, offsets = np.divmod(rows, self.block_rows)
        for block in np.unique(blocks):
            mask = blocks == block
            result[mask] = self._read(int(block))[offsets[mask]]
        return result


@dataclass(frozen=True)
class _MatchedControlRows:
    plate_values: np.ndarray
    match_rows: np.ndarray

    def __getitem__(self, rows):
        return self.plate_values[self.match_rows[rows]]


@dataclass(frozen=True)
class _ModelInput:
    control: _MatchedControlRows
    expression: _ExpressionRows
    _features: np.ndarray | None
    _doses: np.ndarray

    def fingerprints(self, rows):
        if self._features is None:
            raise ValueError("frozen KPGT input is missing")
        return self._features[rows]

    def doses(self, rows):
        return self._doses[rows]


class AxisInputs:
    """Bind pooled conditions and their frozen partitions to the native adapter."""

    def __init__(self, release, representations=None):
        self.release = release
        self.genes = release.genes
        self.treatments = release.treatments
        self.doses = self.treatments.dose_um.to_numpy(dtype=np.float64)
        self.parents = self.treatments.drug_parent.to_numpy(dtype=str)
        match_rows = self.treatments.control_row.to_numpy()
        self._expression = _ExpressionRows(release, "treatment")
        self._measured_delta = _ExpressionRows(release, "delta")
        plates = np.concatenate([
            release.read_controls(start, min(start + 64, len(release.controls)))
            for start in range(0, len(release.controls), 64)
        ])
        self._controls = _MatchedControlRows(plates, match_rows)
        self._features = None
        if representations is not None:
            parents, inverse = np.unique(self.parents, return_inverse=True)
            features = representations.vectors("KPGT", parents)
            if features.shape != (len(parents), WIDTHS["KPGT"]):
                raise ValueError("invalid frozen tenv KPGT width")
            self._features = features[inverse]

    def rows(self, fold, partition):
        if fold not in range(3) or partition not in {"train", "valid", "test"}:
            raise ValueError("unknown cold-drug partition")
        return np.flatnonzero(self.treatments[f"cold_drug_{fold}"].eq(partition).to_numpy())

    def model_input(self, partition):
        if partition not in {"train", "valid", "test"}:
            raise ValueError("unknown model partition")
        return _ModelInput(self._controls, self._expression, self._features, self.doses)

    def condition_ids(self, rows):
        columns = self.treatments.iloc[rows][["plate", "drug_parent", "dose_um"]]
        return tuple(json.dumps((plate, parent, float(dose)), ensure_ascii=False, separators=(",", ":"))
                     for plate, parent, dose in columns.itertuples(index=False, name=None))

    def condition_table(self, rows):
        table = self.treatments.iloc[rows][["plate", "drug_parent", "dose_um"]].copy()
        table.insert(0, "condition_id", self.condition_ids(rows))
        return table


@contextmanager
def open_inputs(training=None, representations=None):
    """Verify the Stage 03 dataset and the accepted tenv KPGT release."""
    with open_tenv_training(training) as release:
        vectors = RepresentationTable(representations, source="tenv")
        yield AxisInputs(release, vectors), vectors


@dataclass(frozen=True)
class PredictionTargets:
    condition_ids: tuple[str, ...]
    predicted_delta: np.ndarray
    measured_delta: np.ndarray


def score_fixed_prediction(data, fold, prediction, condition_ids, *, rows=None):
    """Pair held-out predictions with plate-matched truth in registered row order."""
    test = data.rows(fold, "test")
    rows = test if rows is None else rows
    invalid_rows = "prediction rows must be ordered members of the cold-drug test partition"
    positions = []
    try:
        if np.ndim(rows) != 1:
            raise TypeError("expected a vector of row positions")
        for row in rows:
            if isinstance(row, (bool, np.bool_)):
                raise TypeError("boolean row position")
            positions.append(index(row))
    except (TypeError, ValueError) as exc:
        raise ValueError(invalid_rows) from exc
    if not positions or min(positions) < 0 or max(positions) >= len(data.treatments):
        raise ValueError(invalid_rows)
    rows = np.asarray(positions, dtype=np.int64)
    if np.any(rows[1:] <= rows[:-1]) or not np.isin(rows, test).all():
        raise ValueError(invalid_rows)
    expected = data.condition_ids(rows)
    if tuple(condition_ids) != expected:
        raise ValueError("prediction condition identities differ from the registered held-out rows")
    values = np.asarray(prediction)
    if (values.shape != (len(rows), len(data.genes)) or values.dtype != np.float32
            or not np.isfinite(values).all()):
        raise ValueError("predicted expression must be finite float32 on the released gene axis")
    return PredictionTargets(expected, values - data.model_input("test").control[rows], data._measured_delta[rows])
