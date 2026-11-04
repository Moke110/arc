"""Frozen-input data access for both stages.

Reads ticket #3's frozen input contract:

``cm_transfer.h5ad`` / ``folds.parquet`` / ``consensus_controls.npz`` /
``consensus_membership.parquet`` / ``tahoe_transfer.h5ad`` /
``tahoe_control_groups.parquet`` / ``tahoe_pairs.npz`` / ``locked_config.json`` /
``manifest.json``.

Leakage boundaries implemented here (and tested in
``tests/test_data_leakage.py``):

* a training loader only ever exposes the ``role=train`` treated rows of its
  fold (plus the #11 eligibility-excluded rows, which have no fold and are never
  scored); asking it for a test/valid row raises;
* the CM validation/test inference input is the fold's **fixed consensus** c*;
  matched controls are only used to compute the truth Δ, and only inside the
  evaluation path (:class:`CMControlPool`);
* c* is reconstructed from the fold's recorded membership and compared to the
  stored array before use;
* the dose label is the author's continuous rFCFP4 x log10(dose_uM+1) encoding;
  no binning, no unit conversion.

The whole-module design keeps *noise policy* out of the data layer: training
noise is drawn by the trainer from its seeded generator; validation and
inference use zero noise.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import cache as cache_mod
from . import config as config_mod
from . import h5io
from . import manifest as manifest_mod
from . import paths
from . import tables

__all__ = [
    "DataError",
    "FrozenInput",
    "FoldTable",
    "load_folds",
    "DrugStructures",
    "CMDataset",
    "CMPredictSource",
    "CMControlPool",
    "TahoeDataset",
    "TahoePredictSource",
    "morgan_fingerprints",
    "dose_scale",
    "REQUIRED_CM_COLUMNS",
    "REQUIRED_TAHOE_COLUMNS",
]

log = logging.getLogger("ptlib.data")

_CM_H5AD = "cm_transfer.h5ad"
_TAHOE_H5AD = "tahoe_transfer.h5ad"
_FOLDS = "folds.parquet"
_CONSENSUS = "consensus_controls.npz"
_MEMBERSHIP = "consensus_membership.parquet"
_TAHOE_GROUPS = "tahoe_control_groups.parquet"
_TAHOE_PAIRS = "tahoe_pairs.npz"

REQUIRED_CM_COLUMNS = (
    "well", "plate", "drug_key", "dose_um", "is_control", "vehicle",
    "qc_pass", "main_set",
)
REQUIRED_TAHOE_COLUMNS = (
    "plate", "cell_line_id", "dose_um", "is_control", "qc_pass", "tahoe_split",
)

_ROLES = ("train", "valid", "test")


class DataError(RuntimeError):
    """The frozen input does not satisfy the contracted schema."""


# --------------------------------------------------------------------------- #
# small column helpers
# --------------------------------------------------------------------------- #
def _column(obs: pd.DataFrame, names: Sequence[str], what: str, required: bool = True):
    for name in names:
        if name in obs.columns:
            return name
    if required:
        raise DataError(
            f"obs has no {what} column (looked for {list(names)}; "
            f"present: {sorted(obs.columns)[:40]})"
        )
    return None


def _bool_column(series: pd.Series) -> np.ndarray:
    if series.dtype == bool:
        return series.to_numpy(dtype=bool)
    values = series.to_numpy()
    if values.dtype.kind in "ifb":
        return values.astype(bool)
    return series.astype(str).str.lower().isin(("true", "t", "1", "yes", "y")).to_numpy()


def dose_scale(dose_um: np.ndarray) -> np.ndarray:
    """Author continuous dose encoding: ``log10(dose_uM + 1)``."""
    return np.log10(np.asarray(dose_um, dtype=np.float64) + 1.0).astype(np.float32)


def morgan_fingerprints(unique_smiles: np.ndarray, n_bits: int = 1024) -> np.ndarray:
    """rFCFP4: ``GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True, nBits)``.

    Same-source copy of ``05/train_prnet.py:611`` / the author
    ``Drug_dose_encoder``.  One row per unique SMILES; the dose scaling is
    applied per row by the dataset.  An unparseable SMILES yields a zero row and
    raises nothing (the caller counts it and the run records the count).
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    table = np.zeros((len(unique_smiles), int(n_bits)), dtype=np.float32)
    failures = 0
    for i, smiles in enumerate(unique_smiles):
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            failures += 1
            continue
        bit_string = AllChem.GetMorganFingerprintAsBitVect(
            mol, 2, useFeatures=True, nBits=int(n_bits)
        ).ToBitString()
        if len(bit_string) != int(n_bits):
            raise DataError(f"RDKit returned {len(bit_string)} bits, expected {n_bits}")
        table[i] = np.fromiter((ch == "1" for ch in bit_string), dtype=np.float32, count=len(bit_string))
    if failures:
        log.warning("%d/%d SMILES failed to parse in RDKit (zero fingerprint)", failures, len(unique_smiles))
    return table


# --------------------------------------------------------------------------- #
# drug structures
# --------------------------------------------------------------------------- #
class DrugStructures:
    """drug_key -> SMILES, resolved from the frozen input, then frozen repo evidence.

    The frozen input contract does not list a structures table, so this looks in
    the frozen input first and then in *frozen, tracked* repo evidence (the CM
    corpus ``compounds.csv`` and the 05 pretrain-view ``compounds.csv``/obs).
    The source of every resolved drug is recorded in the run manifest.  A
    structure-bearing drug with no SMILES is counted and dropped; a scored
    (valid/test) row with no SMILES fails closed.
    """

    def __init__(self, table: Mapping[str, str], sources: Sequence[str], unresolved: Sequence[str]):
        self._table = {str(k): str(v) for k, v in table.items() if isinstance(v, str) and v}
        self.sources = list(sources)
        self.unresolved = list(unresolved)

    def __len__(self) -> int:
        return len(self._table)

    def smiles_for(self, drug_key: Any) -> str | None:
        return self._table.get(str(drug_key))

    @classmethod
    def resolve(
        cls,
        data_dir: Path,
        cm_h5ad: Path,
        tahoe_h5ad: Path,
        cm_obs: pd.DataFrame | None = None,
        tahoe_obs: pd.DataFrame | None = None,
    ) -> "DrugStructures":
        table: dict[str, str] = {}
        sources: list[str] = []

        candidates: list[tuple[str, Path | None, tuple[str, ...], tuple[str, ...]]] = [
            ("frozen input drug_smiles.csv", data_dir / "drug_smiles.csv",
             ("drug_key", "key", "compound_name"), ("smiles", "SMILES", "canonical_smiles")),
            ("frozen input compounds.csv", data_dir / "compounds.csv",
             ("drug_key", "parent_inchikey", "compound_name"), ("smiles", "SMILES")),
            ("repo CM corpus compounds.csv",
             paths.ROOT / "stage2/_frozen/cm-data-v1/corpus/v1/compounds.csv",
             ("drug_key",), ("smiles",)),
            ("repo Tahoe pretrain-view compounds.csv",
             paths.ROOT / "stage2/05_prnet_pretrain/output/pretrain_view/compounds.csv",
             ("parent_inchikey", "drug", "compound_name"), ("smiles",)),
        ]
        for label, path, keys, value_cols in candidates:
            if path is None or not Path(path).exists():
                continue
            try:
                frame = pd.read_csv(path)
            except (OSError, ValueError) as exc:  # pragma: no cover - defensive
                log.warning("cannot read structures from %s: %s", path, exc)
                continue
            keys_present = [k for k in keys if k in frame.columns]
            value_col = next((c for c in value_cols if c in frame.columns), None)
            if not keys_present or value_col is None:
                continue
            added = 0
            for key_col in keys_present:
                for key, value in zip(frame[key_col].astype(str), frame[value_col]):
                    if key in ("", "nan", "None") or key in table:
                        continue
                    if isinstance(value, str) and value and value.lower() != "nan":
                        table[key] = value
                        added += 1
            if added:
                sources.append(f"{label} ({added} keys)")

        # Last resort: a SMILES column carried directly on the frozen obs.
        for label, obs in (("cm obs.smiles", cm_obs), ("tahoe obs.smiles", tahoe_obs)):
            if obs is None or "smiles" not in obs.columns or "drug_key" not in obs.columns:
                continue
            added = 0
            for key, value in zip(obs["drug_key"].astype(str), obs["smiles"]):
                if key in table or key in ("", "nan"):
                    continue
                if isinstance(value, str) and value and value.lower() != "nan":
                    table[key] = value
                    added += 1
            if added:
                sources.append(f"{label} ({added} keys)")

        return cls(table, sources, [])


# --------------------------------------------------------------------------- #
# frozen input
# --------------------------------------------------------------------------- #
@dataclass
class FrozenInput:
    """A verified frozen input directory + its locked config."""

    data_dir: Path
    config: dict[str, Any]
    config_notes: dict[str, Any]
    artifact_hashes: dict[str, str]
    manifest_sha256: str
    hash_record: dict[str, Any]
    var_names: np.ndarray
    gene_measured: np.ndarray
    cm_obs: pd.DataFrame
    tahoe_obs: pd.DataFrame
    cm_h5ad: Path
    tahoe_h5ad: Path
    table_cache: Path

    @property
    def n_genes(self) -> int:
        return int(len(self.var_names))

    @property
    def structure(self) -> dict[str, Any]:
        return dict(self.config["structure"])

    def record(self) -> dict[str, Any]:
        return {
            "dir": str(self.data_dir),
            "artifacts": dict(self.artifact_hashes),
            "manifest_sha256": self.manifest_sha256,
            "locked_config": {
                "path": str(self.data_dir / "locked_config.json"),
                "sha256": manifest_mod.sha256_file(self.data_dir / "locked_config.json"),
                "mapping": self.config_notes.get("mapping", []),
                "defaulted": self.config_notes.get("defaulted", []),
            },
            "gene_axis": {
                "n_genes": self.n_genes,
                "measured": int(np.asarray(self.gene_measured, dtype=bool).sum()),
                "gene_measured_source": (
                    "var.gene_measured" if bool(np.asarray(self.gene_measured).any()) else "all-true"
                ),
            },
            "hash_verification": {
                "verified": len(self.hash_record.get("verified", {})),
                "cached": list(self.hash_record.get("cached", [])),
                "bytes": self.hash_record.get("bytes"),
            },
        }

    def structures(self) -> DrugStructures:
        return DrugStructures.resolve(
            self.data_dir, self.cm_h5ad, self.tahoe_h5ad, self.cm_obs, self.tahoe_obs
        )


def open_frozen_input(
    data_dir: Path | str,
    hash_cache: Path | str | None = None,
    verify: bool = True,
    table_cache: Path | str | None = None,
) -> FrozenInput:
    """Verify the frozen artifacts (fail closed) and open their metadata."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise DataError(f"frozen input directory not found: {data_dir}")
    manifest_file = data_dir / "manifest.json"
    if not manifest_file.exists():
        raise DataError(f"frozen input manifest is missing: {manifest_file}")
    raw = manifest_mod.read_json(manifest_file)
    artifact_hashes = manifest_mod.normalize_artifact_manifest(raw)
    if verify:
        hash_record = manifest_mod.verify_artifacts(
            data_dir, artifact_hashes, hash_cache=hash_cache
        )
    else:  # pragma: no cover - explicit escape hatch for tests only
        log.warning("frozen artifact hashing DISABLED (verify=False); this run is not evidence")
        hash_record = {"verified": {}, "cached": [], "bytes": 0}

    config, notes = config_mod.read_locked_config(data_dir / "locked_config.json")
    cm_h5ad = data_dir / _CM_H5AD
    tahoe_h5ad = data_dir / _TAHOE_H5AD
    cm_meta = h5io.open_obs_var(cm_h5ad)
    tahoe_meta = h5io.open_obs_var(tahoe_h5ad)

    n_genes = int(config["structure"]["n_genes"])
    for label, meta in (("cm_transfer", cm_meta), ("tahoe_transfer", tahoe_meta)):
        if meta.n_vars != n_genes:
            raise DataError(
                f"{label}.h5ad has {meta.n_vars} genes but locked_config structure.n_genes={n_genes}"
            )
    var_names = np.asarray(cm_meta.var.index.astype(str))
    tahoe_var = np.asarray(tahoe_meta.var.index.astype(str))
    if not np.array_equal(var_names, tahoe_var):
        raise DataError("cm_transfer and tahoe_transfer are not on the same gene axis/order")
    gene_measured = np.ones(n_genes, dtype=bool)
    if "gene_measured" in cm_meta.var.columns:
        gene_measured = _bool_column(cm_meta.var["gene_measured"])

    return FrozenInput(
        data_dir=data_dir,
        config=config,
        config_notes=notes,
        artifact_hashes=artifact_hashes,
        manifest_sha256=manifest_mod.sha256_file(manifest_file),
        hash_record=hash_record,
        var_names=var_names,
        gene_measured=gene_measured,
        cm_obs=cm_meta.obs,
        tahoe_obs=tahoe_meta.obs,
        cm_h5ad=cm_h5ad,
        tahoe_h5ad=tahoe_h5ad,
        table_cache=(
            Path(table_cache) if table_cache is not None else tables.DEFAULT_TABLE_CACHE
        ),
    )


# --------------------------------------------------------------------------- #
# folds
# --------------------------------------------------------------------------- #
@dataclass
class FoldTable:
    """The explicit ``drug_key -> fold/role`` table (main-set drugs only)."""

    frame: pd.DataFrame
    n_folds: int
    report: dict[str, Any] = field(default_factory=dict)

    @property
    def drugs(self) -> np.ndarray:
        return np.asarray(self.frame["drug_key"].astype(str).unique())

    def role_map(self, fold: int) -> dict[str, str]:
        sub = self.frame[self.frame["fold"] == int(fold)]
        return dict(zip(sub["drug_key"].astype(str), sub["role"].astype(str)))

    def drugs_for(self, fold: int, role: str) -> np.ndarray:
        sub = self.frame[(self.frame["fold"] == int(fold)) & (self.frame["role"] == str(role))]
        return np.asarray(sub["drug_key"].astype(str))

    def validate(self, n_folds: int) -> dict[str, Any]:
        for column in ("drug_key", "fold", "role"):
            if column not in self.frame.columns:
                raise DataError(f"folds.parquet has no {column!r} column")
        frame = self.frame
        roles = set(frame["role"].astype(str))
        unknown = roles - set(_ROLES)
        if unknown:
            raise DataError(f"folds.parquet has unknown role(s) {sorted(unknown)}")
        folds = sorted(int(x) for x in frame["fold"].unique())
        if folds != list(range(int(n_folds))):
            raise DataError(
                f"folds.parquet covers folds {folds}; locked config says n_folds={n_folds}"
            )
        for fold in folds:
            sub = frame[frame["fold"] == fold]
            if sub["drug_key"].duplicated().any():
                raise DataError(f"fold {fold} assigns a drug more than one role")
            for role in _ROLES:
                if int((sub["role"] == role).sum()) == 0 and role == "test":
                    raise DataError(f"fold {fold} has no test drug")
        for role in ("test", "valid"):
            counts = frame[frame["role"] == role].groupby("drug_key").size()
            if (counts != 1).any():
                bad = counts[counts != 1].index.tolist()[:5]
                raise DataError(f"main-set drugs must be {role} exactly once per drug; {bad}")
        counts = frame.groupby("drug_key").size()
        if (counts != len(folds)).any():
            bad = counts[counts != len(folds)].index.tolist()[:5]
            raise DataError(f"every main-set drug must appear in every fold; {bad}")
        return {
            "n_folds": len(folds),
            "n_main_drugs": int(frame["drug_key"].nunique()),
            "per_fold": {
                str(fold): {
                    role: int((frame[(frame["fold"] == fold) & (frame["role"] == role)]).shape[0])
                    for role in _ROLES
                }
                for fold in folds
            },
        }

    @classmethod
    def load(
        cls,
        path: Path | str,
        n_folds: int,
        cache_dir: Path | str | None = None,
        expected_sha256: str | None = None,
    ) -> "FoldTable":
        frame = tables.read_table(path, expected_sha256=expected_sha256, cache_dir=cache_dir)
        table = cls(frame=frame, n_folds=int(n_folds))
        table.report = table.validate(n_folds)
        return table


def load_folds(frozen: FrozenInput) -> FoldTable:
    """Load ``folds.parquet`` through the parquet/CSV-cache reader."""
    return FoldTable.load(
        frozen.data_dir / _FOLDS,
        int(frozen.config["n_folds"]),
        cache_dir=frozen.table_cache,
        expected_sha256=frozen.artifact_hashes.get(_FOLDS),
    )


# --------------------------------------------------------------------------- #
# consensus controls c*
# --------------------------------------------------------------------------- #
def _load_consensus_array(path: Path | str, fold: int, n_folds: int, n_genes: int) -> np.ndarray:
    payload = np.load(path)
    keys = list(payload.keys())

    def _coerce(array: Any, label: str) -> np.ndarray:
        arr = np.asarray(array, dtype=np.float32)
        if arr.shape != (n_genes,):
            raise DataError(f"consensus array {label} has shape {arr.shape}, expected {(n_genes,)}")
        return arr

    for key in (f"fold{fold}", f"fold_{fold}", f"c_star_fold{fold}", f"c_star_{fold}"):
        if key in keys:
            return _coerce(payload[key], key)
    if "folds" in keys and "c_star" in keys:
        folds = np.asarray(payload["folds"]).astype(int).tolist()
        matrix = np.asarray(payload["c_star"], dtype=np.float32)
        if matrix.shape != (len(folds), n_genes):
            raise DataError(f"consensus c_star shape {matrix.shape} != {(len(folds), n_genes)}")
        if fold not in folds:
            raise DataError(f"consensus c_star has no fold {fold} (folds={folds})")
        return matrix[folds.index(fold)]
    if "c_star" in keys:
        matrix = np.asarray(payload["c_star"], dtype=np.float32)
        if matrix.ndim == 2 and matrix.shape == (n_folds, n_genes):
            return matrix[int(fold)]
        if matrix.ndim == 2 and matrix.shape == (n_genes, n_folds):
            return matrix[:, int(fold)]
    raise DataError(
        f"consensus_controls.npz has no per-fold c* for fold {fold}; keys={keys}. "
        "Expected c_star (n_folds, n_genes) + folds, or c_star_fold<k>."
    )


# --------------------------------------------------------------------------- #
# CM dataset
# --------------------------------------------------------------------------- #
@dataclass
class CMDataset:
    """One fold's CM training/validation view with the leakage boundaries enforced.

    * ``train_rows`` — role=train treated wells of this fold (plus #11's expanded
      eligibility-excluded wells when ``train_drugs='expanded'``).  The control
      input is the same-plate matched control mean built from the fold's
      consensus membership.
    * ``valid_rows`` — role=valid treated wells; the model input is c*.
    * ``test_rows`` — role=test treated wells; **never read by this class in a
      training run** (``forbidden_rows`` fails closed on any attempt).
    """

    frozen: FrozenInput
    fold: int
    fold_table: FoldTable
    structures: DrugStructures
    train_drugs: str = "main"
    load_roles: tuple[str, ...] = ("train", "valid")
    forbidden_rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))

    # filled by ``prepare``
    obs: pd.DataFrame = field(default=None)  # type: ignore[assignment]
    c_star: np.ndarray = field(default=None)  # type: ignore[assignment]
    membership: pd.DataFrame = field(default=None)  # type: ignore[assignment]
    plate_control: dict[str, np.ndarray] = field(default_factory=dict)
    train_rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    valid_rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    test_rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    expanded_rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    counts: dict[str, Any] = field(default_factory=dict)
    consensus_check: dict[str, Any] = field(default_factory=dict)

    # derived
    _matrix: np.ndarray | None = None
    _row_pos: dict[int, int] = field(default_factory=dict)
    _fps: np.ndarray | None = None
    _fp_ix: np.ndarray | None = None
    _dose_scale: np.ndarray | None = None
    _control: np.ndarray | None = None
    _rows_read: int = 0

    # -- construction ------------------------------------------------------
    def prepare(self) -> "CMDataset":
        obs = self.frozen.cm_obs
        required = {name: name for name in REQUIRED_CM_COLUMNS}
        for name in required:
            if name not in obs.columns:
                raise DataError(f"cm_transfer.h5ad obs has no {name!r} column")
        self.obs = obs

        treated = ~_bool_column(obs["is_control"])
        qc_pass = _bool_column(obs["qc_pass"])
        main_set = _bool_column(obs["main_set"]) if "main_set" in obs.columns else np.ones(len(obs), bool)
        dose = pd.to_numeric(obs["dose_um"], errors="coerce").to_numpy(dtype=np.float64)
        plate = obs["plate"].astype(str).to_numpy()
        drug_key = obs["drug_key"].astype(str).to_numpy()
        okay = treated & qc_pass & np.isfinite(dose) & (dose >= 0)

        # structures
        smiles = np.array(
            [self.structures.smiles_for(key) or "" for key in drug_key], dtype=object
        )
        has_smiles = np.array([bool(s) for s in smiles], dtype=bool)

        # folds authority
        role_map = self.fold_table.role_map(self.fold)
        role_of_row = np.array([role_map.get(key, "") for key in drug_key], dtype=object)

        rows_by_role: dict[str, np.ndarray] = {}
        missing_structure: dict[str, int] = {}
        for role in ("train", "valid", "test"):
            mask = okay & (role_of_row == role)
            rows = np.flatnonzero(mask).astype(np.int64)
            missing = int((~has_smiles[rows]).sum())
            if missing:
                missing_structure[role] = missing
                if role in ("valid", "test"):
                    bad = sorted({drug_key[r] for r in rows if not has_smiles[r]})[:5]
                    raise DataError(
                        f"fold {self.fold} has {missing} {role} row(s) without a resolvable SMILES "
                        f"(drug_key examples {bad}); the dose label cannot be encoded. "
                        "Fix the frozen structures, do not drop scored rows silently."
                    )
                rows = rows[has_smiles[rows]]
            rows_by_role[role] = rows

        expanded_rows = np.empty(0, dtype=np.int64)
        expanded_counts: dict[str, Any] = {}
        if self.train_drugs not in ("main", "expanded"):
            raise DataError(f"train_drugs must be 'main' or 'expanded', got {self.train_drugs!r}")
        if self.train_drugs == "expanded":
            in_folds = np.isin(drug_key, self.fold_table.drugs)
            eligibility = (
                obs["eligibility_decision"].astype(str).to_numpy()
                if "eligibility_decision" in obs.columns
                else np.full(len(obs), "unknown", dtype=object)
            )
            mask = okay & ~main_set & ~in_folds & has_smiles
            expanded_rows = np.flatnonzero(mask).astype(np.int64)
            expanded_counts = {
                "rows": int(expanded_rows.size),
                "eligibility_decision": dict(
                    pd.Series(eligibility[expanded_rows]).value_counts().astype(int).to_dict()
                ),
                "missing_smiles": int((okay & ~main_set & ~in_folds & ~has_smiles).sum()),
            }
            if expanded_rows.size == 0:
                raise DataError(
                    "train_drugs='expanded' but no eligibility-excluded QC-pass treated well "
                    "with a structure was found (#11 needs them); refusing to pretend otherwise"
                )

        self.train_rows = np.concatenate([rows_by_role["train"], expanded_rows])
        self.train_rows = np.unique(self.train_rows)
        self.valid_rows = rows_by_role["valid"]
        self.test_rows = rows_by_role["test"]
        self.expanded_rows = expanded_rows
        self._role_rows = rows_by_role
        if self.train_rows.size == 0:
            raise DataError(f"fold {self.fold} has no training rows")
        if self.valid_rows.size == 0:
            raise DataError(f"fold {self.fold} has no validation rows")
        self.forbidden_rows = np.unique(self.test_rows)

        # consensus c*: reconstruct from the recorded membership and verify.
        self.membership = tables.read_table(
            self.frozen.data_dir / _MEMBERSHIP,
            expected_sha256=self.frozen.artifact_hashes.get(_MEMBERSHIP),
            cache_dir=self.frozen.table_cache,
        )
        matrix_rows = np.unique(np.concatenate([self.train_rows, self.valid_rows]))
        self.rows = matrix_rows
        self._row_pos = {int(row): i for i, row in enumerate(matrix_rows)}
        self._load_consensus(membership_required_plates=set(plate[self.train_rows]))

        # per-row dose scale and drug labels
        self._dose_scale = dose_scale(dose[self.rows])
        smiles_sel = smiles[self.rows]
        unique_smiles, fp_ix = np.unique(smiles_sel, return_inverse=True)
        n_bits = int(self.frozen.structure["drug_dim"])
        self._fps = morgan_fingerprints(unique_smiles, n_bits)
        self._fp_ix = fp_ix.astype(np.int32)
        self._plate = plate[self.rows]
        self._drug = drug_key[self.rows]

        # training control = per-plate matched control mean from membership;
        # validation control = the fold's fixed consensus.
        control = np.empty((self.rows.size, self.frozen.n_genes), dtype=np.float32)
        per_plate_rows = np.flatnonzero(np.isin(self.rows, self.train_rows))
        for row in per_plate_rows:
            key = plate[self.rows[row]]
            if key not in self.plate_control:
                raise DataError(
                    f"training well on plate {key!r} has no consensus-membership control well; "
                    "the frozen training control pool is incomplete (fail closed)"
                )
            control[row] = self.plate_control[key]
        valid_plate_rows = np.flatnonzero(np.isin(self.rows, self.valid_rows))
        control[valid_plate_rows] = self.c_star
        self._control = control

        self.counts = {
            "fold": int(self.fold),
            "train_rows_main": int(rows_by_role["train"].size),
            "train_rows_expanded": int(expanded_rows.size),
            "train_rows": int(self.train_rows.size),
            "valid_rows": int(self.valid_rows.size),
            "test_rows": int(self.test_rows.size),
            "n_plates": int(pd.unique(self._plate).size),
            "n_conditions_valid": int(self.condition_ids(self.valid_rows).size),
            "missing_smiles_train_dropped": int(missing_structure.get("train", 0)),
            "expanded": expanded_counts,
        }
        return self

    def _load_consensus(self, membership_required_plates: set[str]) -> None:
        n_genes = self.frozen.n_genes
        n_folds = int(self.frozen.config.get("n_folds", 5))
        self.c_star = _load_consensus_array(
            self.frozen.data_dir / _CONSENSUS, self.fold, n_folds, n_genes
        )
        frame = self.membership
        needed = ("fold", "plate", "row")
        for column in needed:
            if column not in frame.columns:
                raise DataError(
                    f"consensus_membership.parquet has no {column!r} column (has {list(frame.columns)})"
                )
        fold_frame = frame[frame["fold"].astype(int) == int(self.fold)]
        if fold_frame.empty:
            raise DataError(f"consensus membership has no rows for fold {self.fold}")
        rows = np.asarray(fold_frame["row"], dtype=np.int64)
        if rows.min() < 0 or rows.max() >= len(self.frozen.cm_obs):
            raise DataError("consensus membership references rows outside cm_transfer.h5ad")
        log1p, sorted_rows = h5io.read_layer_rows(
            self.frozen.cm_h5ad, np.unique(rows), n_genes
        )
        pos = {int(row): i for i, row in enumerate(sorted_rows)}
        matrix = log1p[[pos[int(row)] for row in rows]]
        plates = fold_frame["plate"].astype(str).to_numpy()

        has_weight = "weight" in fold_frame.columns
        weights = (
            np.asarray(fold_frame["weight"], dtype=np.float64) if has_weight else None
        )
        if weights is not None:
            recomputed = np.tensordot(weights, matrix.astype(np.float64), axes=(0, 0))
        else:
            unique_plates = sorted(set(plates))
            plate_means = np.stack(
                [matrix[plates == p].astype(np.float64).mean(axis=0) for p in unique_plates]
            )
            recomputed = plate_means.mean(axis=0)
            weights = None

        diff = float(np.max(np.abs(recomputed.astype(np.float32) - self.c_star)))
        if diff > 1e-4:
            raise DataError(
                f"fold {self.fold}: c* cannot be reconstructed from consensus_membership "
                f"(max abs diff {diff:.3g}); the frozen consensus and its membership disagree"
            )

        unique_plates = sorted(set(plates))
        self.plate_control = {}
        for plate in unique_plates:
            mask = plates == plate
            if weights is not None:
                w = weights[mask]
                self.plate_control[str(plate)] = (
                    (matrix[mask].astype(np.float64) * w[:, None]).sum(axis=0) / w.sum()
                ).astype(np.float32)
            else:
                self.plate_control[str(plate)] = (
                    matrix[mask].astype(np.float64).mean(axis=0).astype(np.float32)
                )
        missing_plates = membership_required_plates - set(self.plate_control)
        if missing_plates:
            raise DataError(
                f"fold {self.fold}: training wells on plate(s) {sorted(missing_plates)} have no "
                "control well in the consensus membership"
            )
        self.consensus_check = {
            "fold": int(self.fold),
            "n_plates": len(unique_plates),
            "n_member_rows": int(len(rows)),
            "max_abs_diff": diff,
            "used_weights": bool(has_weight),
        }

    # -- matrices ----------------------------------------------------------
    def _matrix_array(self) -> np.ndarray:
        if self._matrix is None:
            layer, sorted_rows = h5io.read_layer_rows(
                self.frozen.cm_h5ad, self.rows, self.frozen.n_genes
            )
            if not np.array_equal(sorted_rows, self.rows):
                raise DataError("CM row read returned an unexpected row order")
            self._matrix = layer
        return self._matrix

    def _positions(self, rows: np.ndarray) -> np.ndarray:
        """Map obs row ids to this fold's train/valid positions, failing closed."""
        rows = np.asarray(rows, dtype=np.int64)
        forbidden = np.intersect1d(rows, self.forbidden_rows, assume_unique=False)
        if forbidden.size:
            raise DataError(
                f"training loader refused to read {forbidden.size} held-out (test) treated row(s) "
                f"of fold {self.fold}: {forbidden[:5].tolist()}"
            )
        missing = [int(row) for row in rows if int(row) not in self._row_pos]
        if missing:
            raise DataError(
                f"row(s) are not in fold {self.fold}'s train/valid row set: {missing[:5]}"
            )
        return np.array([self._row_pos[int(row)] for row in rows], dtype=np.int64)

    def target(self, rows: np.ndarray) -> np.ndarray:
        pos = self._positions(rows)
        matrix = self._matrix_array()
        self._rows_read += int(pos.size)
        return matrix[pos]

    # -- batching ----------------------------------------------------------
    def condition_key(self, rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        drug = self.frozen.cm_obs["drug_key"].astype(str).to_numpy()[rows]
        dose = pd.to_numeric(self.frozen.cm_obs["dose_um"], errors="coerce").to_numpy()[rows]
        return np.asarray([f"{d}|{dose:g}" for d, dose in zip(drug, dose)], dtype=object)

    def condition_ids(self, rows: np.ndarray) -> np.ndarray:
        keys = self.condition_key(rows)
        unique, inverse = np.unique(keys, return_inverse=True)
        self._last_condition_names = unique
        return inverse

    def batch(self, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(control, label, target)`` for treated *rows* (model input)."""
        rows = np.asarray(rows, dtype=np.int64)
        pos = self._positions(rows)
        control = self._control[pos]
        label = self._fps[self._fp_ix[pos]].copy()
        label *= self._dose_scale[pos, None]
        target = self.target(rows)
        return control, label, target

    def baseline(self, rows: np.ndarray) -> np.ndarray:
        return self._control[self._positions(rows)]

    @property
    def rows_read(self) -> int:
        """Number of treated-expression row reads this loader performed (test hook)."""
        return int(self._rows_read)


# --------------------------------------------------------------------------- #
# CM prediction source (target-free)
# --------------------------------------------------------------------------- #
class CMPredictSource:
    """The target-free CM inference view: obs rows + the fold's fixed c*.

    Deliberately has no access to treated expression: ``predict.py`` can only
    ask for row identities, drug labels and the c* baseline.  This is what makes
    "the model never receives a test/held-out control pool" and "prediction is
    target-free" structural rather than aspirational.
    """

    def __init__(self, frozen: FrozenInput, fold_table: FoldTable, fold: int, structures: DrugStructures):
        self.frozen = frozen
        self.fold_table = fold_table
        self.fold = int(fold)
        self.structures = structures
        obs = frozen.cm_obs
        for column in ("is_control", "qc_pass", "drug_key", "dose_um", "plate", "main_set"):
            if column not in obs.columns:
                raise DataError(f"cm_transfer.h5ad obs has no {column!r} column")
        n_genes = frozen.n_genes
        n_folds = int(frozen.config.get("n_folds", 5))
        self.c_star = _load_consensus_array(
            frozen.data_dir / _CONSENSUS, self.fold, n_folds, n_genes
        )
        self.fold_table.validate(n_folds)

    def rows_for(self, role: str) -> np.ndarray:
        if role not in _ROLES:
            raise DataError(f"role must be one of {_ROLES}, got {role!r}")
        obs = self.frozen.cm_obs
        treated = ~_bool_column(obs["is_control"])
        qc_pass = _bool_column(obs["qc_pass"])
        dose = pd.to_numeric(obs["dose_um"], errors="coerce").to_numpy(dtype=np.float64)
        drug_key = obs["drug_key"].astype(str).to_numpy()
        role_map = self.fold_table.role_map(self.fold)
        role_of_row = np.array([role_map.get(key, "") for key in drug_key], dtype=object)
        smiles = np.array(
            [self.structures.smiles_for(key) or "" for key in drug_key], dtype=object
        )
        mask = treated & qc_pass & np.isfinite(dose) & (dose >= 0) & (role_of_row == role)
        mask &= np.array([bool(s) for s in smiles], dtype=bool)
        return np.flatnonzero(mask).astype(np.int64)

    def condition_key(self, rows: np.ndarray) -> np.ndarray:
        obs = self.frozen.cm_obs
        drug = obs["drug_key"].astype(str).to_numpy()[rows]
        dose = pd.to_numeric(obs["dose_um"], errors="coerce").to_numpy()[rows]
        return np.asarray([f"{d}|{dose:g}" for d, dose in zip(drug, dose)], dtype=object)

    def meta(self, rows: np.ndarray) -> pd.DataFrame:
        obs = self.frozen.cm_obs.reset_index(drop=True)
        frame = obs.iloc[np.asarray(rows, dtype=np.int64)][
            ["well", "plate", "drug_key", "dose_um"]
        ].copy()
        frame.insert(0, "row", np.asarray(rows, dtype=np.int64))
        frame["condition"] = self.condition_key(np.asarray(rows))
        return frame.reset_index(drop=True)

    def label(self, rows: np.ndarray) -> np.ndarray:
        """rFCFP4 x log10(dose+1) for the requested rows."""
        obs = self.frozen.cm_obs.reset_index(drop=True)
        keys = obs["drug_key"].astype(str).to_numpy()[np.asarray(rows, dtype=np.int64)]
        smiles = np.array([self.structures.smiles_for(k) or "" for k in keys], dtype=object)
        if any(not s for s in smiles):
            raise DataError("cannot encode a dose label for a drug without a SMILES")
        dose = pd.to_numeric(
            obs["dose_um"], errors="coerce"
        ).to_numpy()[np.asarray(rows, dtype=np.int64)]
        unique, inverse = np.unique(smiles, return_inverse=True)
        fps = morgan_fingerprints(unique, int(self.frozen.structure["drug_dim"]))
        label = fps[inverse].copy()
        label *= dose_scale(dose)[:, None]
        return label

    def baseline(self, rows: np.ndarray) -> np.ndarray:
        return np.broadcast_to(self.c_star, (len(rows), self.frozen.n_genes)).copy()


# --------------------------------------------------------------------------- #
# CM matched-control pool (evaluation only)
# --------------------------------------------------------------------------- #
class CMControlPool:
    """Matched controls for the experimental truth Δ — **evaluation only**.

    The truth Δ of a validation/test well is ``treated log1p - mean(log1p of the
    matched DMSO wells)``, where "matched" means the same plate (and the same
    vehicle); it is computed here and never exposed as a model input.  The model
    input at validation/test is always the fold's fixed c*.

    This reads the control wells' expression once (a few hundred rows).  It is
    instantiated only by ``ptlib.evaluate`` / the trainer's validation scorer.
    """

    def __init__(self, frozen: FrozenInput):
        self.frozen = frozen
        obs = frozen.cm_obs
        is_control = _bool_column(obs["is_control"])
        qc_pass = _bool_column(obs["qc_pass"])
        vehicle_col = _column(obs, ("vehicle", "solvent"), "vehicle/solvent", required=False)
        vehicle = (
            obs[vehicle_col].astype(str).to_numpy()
            if vehicle_col is not None
            else np.full(len(obs), "DMSO", dtype=object)
        )
        plate = obs["plate"].astype(str).to_numpy()
        mask = is_control & qc_pass & np.isin(vehicle, ("DMSO", "dmso", "Dmso"))
        rows = np.flatnonzero(mask).astype(np.int64)
        if rows.size == 0:
            raise DataError("no QC-pass DMSO control wells for the matched-control truth Δ")
        layer, sorted_rows = h5io.read_layer_rows(frozen.cm_h5ad, rows, frozen.n_genes)
        pos = {int(row): i for i, row in enumerate(sorted_rows)}
        self.control_rows = rows
        self.plate_means: dict[str, np.ndarray] = {}
        for plate_key in pd.unique(plate[rows]):
            selected = [pos[int(row)] for row in rows if plate[row] == plate_key]
            self.plate_means[str(plate_key)] = (
                layer[selected].astype(np.float64).mean(axis=0).astype(np.float32)
            )

    def mean_for_rows(self, rows: np.ndarray) -> np.ndarray:
        obs = self.frozen.cm_obs
        plate = obs["plate"].astype(str).to_numpy()
        out = np.empty((len(rows), self.frozen.n_genes), dtype=np.float32)
        for i, row in enumerate(np.asarray(rows, dtype=np.int64)):
            key = str(plate[row])
            if key not in self.plate_means:
                raise DataError(
                    f"no QC-pass matched DMSO control on plate {key!r} for the truth Δ "
                    "(evaluation path fails closed; it never borrows another plate's control)"
                )
            out[i] = self.plate_means[key]
        return out


# --------------------------------------------------------------------------- #
# Tahoe dataset
# --------------------------------------------------------------------------- #
class _TahoeMatrix:
    """Row access for the pretrain matrix, in one of three modes.

    ``cache`` — derive-once dense float32 ``.npy`` read through a shared memmap
    (the #5 memory rule; page cache is shared across concurrent seeds);
    ``csr`` — the whole frozen CSR layer in process memory (the 05 pattern);
    ``dense`` — an in-memory dense subset (smoke/subset runs).
    """

    def __init__(self, mode: str, source=None, rows: np.ndarray | None = None, shape=None):
        self.mode = mode
        self.source = source
        self.rows_index = None if rows is None else {int(r): i for i, r in enumerate(rows)}
        self.shape = shape

    def read(self, rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        if self.mode == "dense":
            index = np.array([self.rows_index[int(r)] for r in rows], dtype=np.int64)
            return np.asarray(self.source[index], dtype=np.float32)
        if self.mode == "cache":
            # the cache holds all 67k rows in frozen obs order: index directly
            return np.asarray(self.source.rows(rows), dtype=np.float32)
        if self.mode == "csr":
            # the full sparse layer holds all obs rows in frozen obs order
            return self.source[rows].toarray().astype(np.float32, copy=False)
        raise DataError(f"unknown Tahoe matrix mode {self.mode!r}")


@dataclass
class TahoeDataset:
    """Pretrain view (Tahoe) with matched ``(plate, cell_line)`` control inputs.

    Selection: treated rows with a structure, QC pass, finite non-negative dose
    and a resolvable SMILES.  Split authority is the frozen ``tahoe_split``
    column (plate13/14 = valid per the locked config; not a cold-drug claim).
    ``max_train_rows`` / ``max_valid_rows`` sample a deterministic subset for the
    smoke run.  Control means never cross cell lines: the group key is
    ``(plate, cell_line_id)``.
    """

    frozen: FrozenInput
    structures: DrugStructures
    max_train_rows: int | None = None
    max_valid_rows: int | None = None
    seed: int = 20261021
    cache: cache_mod.TahoeCache | None = None
    full: bool = False
    load_targets: bool = True

    obs: pd.DataFrame = field(default=None)  # type: ignore[assignment]
    rows_train: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    rows_valid: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    rows: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    control_values: np.ndarray = field(default=None)  # type: ignore[assignment]
    group_index: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    counts: dict[str, Any] = field(default_factory=dict)
    matrix: _TahoeMatrix | None = None
    _row_pos: dict[int, int] = field(default_factory=dict)
    _fps: np.ndarray | None = None
    _fp_ix: np.ndarray | None = None
    _dose_scale: np.ndarray | None = None
    _condition: np.ndarray | None = None
    cache_record: dict[str, Any] = field(default_factory=dict)

    def prepare(self) -> "TahoeDataset":
        obs = self.frozen.tahoe_obs
        for column in ("plate", "cell_line_id", "dose_um", "is_control", "tahoe_split"):
            if column not in obs.columns:
                raise DataError(f"tahoe_transfer.h5ad obs has no {column!r} column")
        self.obs = obs.reset_index(drop=True)
        n = len(obs)
        is_control = _bool_column(obs["is_control"])
        qc_pass = (
            _bool_column(obs["qc_pass"]) if "qc_pass" in obs.columns else np.ones(n, bool)
        )
        has_structure = (
            _bool_column(obs["has_structure"])
            if "has_structure" in obs.columns
            else np.ones(n, bool)
        )
        dose = pd.to_numeric(obs["dose_um"], errors="coerce").to_numpy(dtype=np.float64)
        split = obs["tahoe_split"].astype(str).to_numpy()
        drug_col = _column(obs, ("drug_key", "compound_name"), "drug_key/compound_name")
        drug_key = obs[drug_col].astype(str).to_numpy()

        smiles = np.array([self.structures.smiles_for(k) or "" for k in drug_key], dtype=object)
        has_smiles = np.array([bool(s) for s in smiles], dtype=bool)
        if "smiles" in obs.columns:
            for i, (s, k) in enumerate(zip(obs["smiles"], drug_key)):
                if not has_smiles[i] and isinstance(s, str) and s and s.lower() != "nan":
                    smiles[i] = s
                    has_smiles[i] = True

        # Frozen usability mask when the artifact carries it: ticket #3's
        # ``tahoe_usable`` = ~is_control AND has_structure AND qc_pass AND
        # matched_control_ok.  Re-deriving it here silently re-includes the 4,757
        # treated rows that have no QC-pass matched control (they have no baseline,
        # and the contract says those are reported, never imputed).  The synthetic
        # fixture has no such column, so the derived mask is its fallback; every
        # fixture row carries a matched control, so both paths agree there.
        frozen_usable = "tahoe_usable" in obs.columns
        usable = (
            _bool_column(obs["tahoe_usable"])
            if frozen_usable
            else (~is_control) & qc_pass & has_structure
        )
        encodable = (
            np.isfinite(dose)
            & (dose >= 0)
            & has_smiles
            & np.isin(split, ("train", "valid"))
        )
        if frozen_usable and (usable & ~encodable).any():
            raise DataError(
                f"tahoe_transfer.h5ad marks {int((usable & ~encodable).sum())} tahoe_usable "
                "row(s) without a finite dose, a resolvable SMILES or a train/valid split; "
                "the frozen usability contract and the loader disagree"
            )
        okay = usable & encodable

        # matched controls: group rows per (plate, cell_line_id); a selected row
        # without a control group has no baseline.  The frozen usable mask already
        # excludes those, so a drop below means the side files disagree.
        groups, index_map = self._control_groups(obs)
        keys = np.asarray(
            [
                f"{plate}|{line}"
                for plate, line in zip(
                    obs["plate"].astype(str).to_numpy(),
                    obs["cell_line_id"].astype(str).to_numpy(),
                )
            ],
            dtype=object,
        )
        has_group = np.isin(keys, list(index_map))
        if frozen_usable and (okay & ~has_group).any():
            # The frozen mask says these rows are usable, yet no QC-pass matched
            # control group exists for their (plate, cell_line).  Never impute a
            # baseline — the two side files disagree, so fail closed.
            raise DataError(
                f"tahoe_transfer.h5ad marks {int((okay & ~has_group).sum())} row(s) "
                "tahoe_usable whose (plate, cell_line_id) has no matched-control group in "
                "tahoe_control_groups.parquet; the frozen artifacts disagree"
            )
        okay = okay & has_group
        counts = {
            "n_obs": int(n),
            "n_controls": int(is_control.sum()),
            "n_perturbed": int((~is_control).sum()),
            "dropped_qc_fail": int((~is_control & ~qc_pass).sum()),
            "dropped_no_structure": int((~is_control & ~has_structure).sum()),
            "dropped_no_smiles": int((~is_control & has_structure & ~has_smiles).sum()),
            "dropped_bad_dose": int((~is_control & ~np.isfinite(dose)).sum()),
            "dropped_bad_split": int((~is_control & ~np.isin(split, ("train", "valid"))).sum()),
            "usable_mask_source": (
                "obs.tahoe_usable" if frozen_usable else "derived(qc_pass&has_structure)"
            ),
            "dropped_no_matched_control": int((usable & encodable & ~has_group).sum()),
        }
        rng = np.random.default_rng(int(self.seed))
        rows_train = np.flatnonzero(okay & (split == "train")).astype(np.int64)
        rows_valid = np.flatnonzero(okay & (split == "valid")).astype(np.int64)
        if self.max_train_rows is not None and rows_train.size > int(self.max_train_rows):
            rows_train = np.sort(
                rng.choice(rows_train, size=int(self.max_train_rows), replace=False)
            )
        if self.max_valid_rows is not None and rows_valid.size > int(self.max_valid_rows):
            rows_valid = np.sort(
                rng.choice(rows_valid, size=int(self.max_valid_rows), replace=False)
            )
        if rows_train.size < 2:
            raise DataError(f"Tahoe pretrain has {rows_train.size} training row(s)")
        if rows_valid.size == 0:
            raise DataError("Tahoe pretrain has no validation rows")
        counts["n_train"] = int(rows_train.size)
        counts["n_valid"] = int(rows_valid.size)
        counts["smoke_subset"] = bool(self.max_train_rows or self.max_valid_rows)

        self.rows_train, self.rows_valid = rows_train, rows_valid
        self.rows = np.unique(np.concatenate([rows_train, rows_valid]))
        self._row_pos = {int(row): i for i, row in enumerate(self.rows)}

        self.group_index = np.array([index_map[keys[row]] for row in self.rows], dtype=np.int64)
        self.control_values = self._control_means(groups, index_map)

        self._dose_scale = dose_scale(dose[self.rows])
        # Hoist the per-row column conversion out of the comprehension: calling
        # ``obs[...].astype(str).to_numpy()`` inside it is O(n_rows x n_obs) and
        # took 134 s on the real 67k-row view (the synthetic fixture hid it).
        cell_line = obs["cell_line_id"].astype(str).to_numpy()
        self._condition = np.asarray(
            [
                f"{drug_key[row]}|{dose[row]:g}|{cell_line[row]}"
                for row in self.rows
            ],
            dtype=object,
        )
        unique_smiles, fp_ix = np.unique(smiles[self.rows], return_inverse=True)
        self._fps = morgan_fingerprints(unique_smiles, int(self.frozen.structure["drug_dim"]))
        self._fp_ix = fp_ix.astype(np.int32)

        if self.load_targets:
            self._prepare_matrix(obs)
        counts["control_groups"] = int(len(index_map))
        counts["pairs_check"] = self._pairs_check
        self.counts = counts
        return self

    # -- controls ----------------------------------------------------------
    def _read_control_groups(self, obs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(member_rows, group_ids)`` for the matched-control groups."""
        path = self.frozen.data_dir / _TAHOE_GROUPS
        if not path.exists():
            raise DataError(f"tahoe_control_groups.parquet is missing: {path}")
        frame = tables.read_table(
            path,
            expected_sha256=self.frozen.artifact_hashes.get(_TAHOE_GROUPS),
            cache_dir=self.frozen.table_cache,
        )
        # Long form: one row per (group, member observation).
        if {"group_index", "row"} <= set(frame.columns) or {"group_key", "row"} <= set(frame.columns):
            rows = np.asarray(frame["row"], dtype=np.int64)
            if "group_index" in frame.columns:
                group_ids = np.asarray(frame["group_index"], dtype=np.int64)
            else:
                group_ids = np.asarray(
                    pd.factorize(frame["group_key"].astype(str))[0], dtype=np.int64
                )
            return rows, group_ids
        raise DataError(
            "tahoe_control_groups.parquet is not in the contracted long form "
            f"(group_index|group_key, row); columns={list(frame.columns)}"
        )

    def _control_groups(self, obs: pd.DataFrame) -> tuple[dict[int, str], dict[str, int]]:
        member_rows, group_ids = self._read_control_groups(obs)
        if member_rows.size == 0:
            raise DataError("tahoe_control_groups.parquet lists no control rows")
        # Map each group id to a canonical (plate, cell_line_id) key, then drop
        # any group that is not a QC-pass DMSO control of the same line/plate.
        is_control = _bool_column(obs["is_control"])
        plate = obs["plate"].astype(str).to_numpy()
        line = obs["cell_line_id"].astype(str).to_numpy()
        qc_pass = _bool_column(obs["qc_pass"]) if "qc_pass" in obs.columns else np.ones(len(obs), bool)
        bad = (~is_control[member_rows]) | (~qc_pass[member_rows])
        if "vehicle" in obs.columns:
            vehicle = obs["vehicle"].astype(str).to_numpy()
            bad = bad | ~np.isin(vehicle[member_rows], ("DMSO_TF", "DMSO", "dmso", "DMSO-tf"))
        if bad.any():
            raise DataError(
                f"tahoe_control_groups references {int(bad.sum())} non-QC-pass/non-DMSO control "
                "row(s); matched controls must be the frozen DMSO_TF wells"
            )
        key_of_group: dict[int, str] = {}
        members_by_key: dict[str, list[int]] = {}
        for row, gid in zip(member_rows, group_ids):
            key = f"{plate[row]}|{line[row]}"
            existing = key_of_group.get(int(gid))
            if existing is not None and existing != key:
                raise DataError(
                    f"Tahoe control group {gid} spans multiple (plate, cell_line) keys "
                    f"({existing!r}, {key!r}); controls must never cross cell lines"
                )
            key_of_group[int(gid)] = key
            members_by_key.setdefault(key, []).append(int(row))
        index_map = {key: i for i, key in enumerate(sorted(members_by_key))}
        self._control_members = {key: np.asarray(sorted(v)) for key, v in members_by_key.items()}
        self._pairs_check = self._verify_pairs(obs, key_of_group)
        return key_of_group, index_map

    def _verify_pairs(self, obs: pd.DataFrame, key_of_group: Mapping[int, str]) -> dict[str, Any]:
        """Cross-check ``tahoe_pairs.npz`` against the group keys (when present)."""
        path = self.frozen.data_dir / _TAHOE_PAIRS
        if not path.exists():
            return {"pairs_present": False}
        payload = np.load(path)
        if "row" not in payload or "group_index" not in payload:
            raise DataError(
                f"tahoe_pairs.npz must hold 'row' and 'group_index'; keys={list(payload.keys())}"
            )
        pair_rows = np.asarray(payload["row"], dtype=np.int64)
        pair_groups = np.asarray(payload["group_index"], dtype=np.int64)
        plate = obs["plate"].astype(str).to_numpy()
        line = obs["cell_line_id"].astype(str).to_numpy()
        mismatches = 0
        for row, gid in zip(pair_rows, pair_groups):
            expected = f"{plate[row]}|{line[row]}"
            actual = key_of_group.get(int(gid))
            if actual != expected:
                mismatches += 1
        if mismatches:
            raise DataError(
                f"tahoe_pairs.npz disagrees with the control-group keys for {mismatches} treated "
                "row(s); the matched-control mapping is inconsistent"
            )
        return {"pairs_present": True, "pairs_checked": int(pair_rows.size)}

    def _control_means(self, member_rows: np.ndarray, index_map: Mapping[str, int]) -> np.ndarray:
        keys = sorted(self._control_members)
        rows = np.unique(np.concatenate([self._control_members[k] for k in keys]))
        layer, sorted_rows = h5io.read_layer_rows(self.frozen.tahoe_h5ad, rows, self.frozen.n_genes)
        pos = {int(row): i for i, row in enumerate(sorted_rows)}
        values = np.empty((len(keys), self.frozen.n_genes), dtype=np.float32)
        for key in keys:
            member_pos = [pos[int(row)] for row in self._control_members[key]]
            values[index_map[key]] = layer[member_pos].astype(np.float64).mean(axis=0).astype(np.float32)
        self._control_keys = keys
        return values

    # -- matrix ------------------------------------------------------------
    def _prepare_matrix(self, obs: pd.DataFrame) -> None:
        if self.cache is not None:
            self.matrix = _TahoeMatrix("cache", source=self.cache)
            self.cache_record = {
                "mode": "cache",
                "path": str(self.cache.path),
                "verify": self.cache.verify_mode,
                "sample_max_abs_diff": self.cache.sample_max_abs_diff,
                "layer_sha256": self.cache.record.get("layer_sha256"),
                "file_sha256": self.cache.record.get("file_sha256"),
                "source_sha256": self.cache.record.get("source_sha256"),
            }
            return
        if self.full:
            sparse = h5io.read_layer_all(
                self.frozen.tahoe_h5ad, len(obs), self.frozen.n_genes
            )
            if not isinstance(sparse, np.ndarray):
                sparse = sparse.tocsr()
            self.matrix = _TahoeMatrix("csr", source=sparse)
            self.cache_record = {"mode": "csr-full"}
            return
        layer, sorted_rows = h5io.read_layer_rows(
            self.frozen.tahoe_h5ad, self.rows, self.frozen.n_genes
        )
        if not np.array_equal(sorted_rows, self.rows):
            raise DataError("Tahoe row read returned an unexpected row order")
        self.matrix = _TahoeMatrix("dense", source=layer, rows=self.rows)
        self.cache_record = {"mode": "dense-subset", "rows": int(self.rows.size)}

    # -- batching ----------------------------------------------------------
    def target(self, rows: np.ndarray) -> np.ndarray:
        if self.matrix is None:
            raise DataError("TahoeDataset was built without targets (predict mode)")
        return self.matrix.read(rows)

    def baseline(self, rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        pos = np.array([self._row_pos[int(row)] for row in rows], dtype=np.int64)
        return self.control_values[self.group_index[pos]]

    def batch(self, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rows = np.asarray(rows, dtype=np.int64)
        pos = np.array([self._row_pos[int(row)] for row in rows], dtype=np.int64)
        control = self.control_values[self.group_index[pos]]
        label = self._fps[self._fp_ix[pos]].copy()
        label *= self._dose_scale[pos, None]
        target = self.target(rows)
        return control, label, target

    def condition_ids(self, rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        keys = self._condition[[self._row_pos[int(r)] for r in rows]]
        unique, inverse = np.unique(keys, return_inverse=True)
        self._last_condition_names = unique
        return inverse


class TahoePredictSource:
    """Target-free Tahoe inference view (obs + matched control inputs)."""

    def __init__(self, frozen: FrozenInput, structures: DrugStructures, cache=None, full: bool = False):
        self.frozen = frozen
        dataset = TahoeDataset(
            frozen=frozen,
            structures=structures,
            cache=cache,
            full=full,
            load_targets=False,
        ).prepare()
        self.dataset = dataset

    def rows_for(self, role: str) -> np.ndarray:
        if role == "train":
            return self.dataset.rows_train
        if role == "valid":
            return self.dataset.rows_valid
        raise DataError("Tahoe has no test split; role must be 'train' or 'valid'")

    def baseline(self, rows: np.ndarray) -> np.ndarray:
        return self.dataset.baseline(rows)

    def label(self, rows: np.ndarray) -> np.ndarray:
        dataset = self.dataset
        pos = np.array([dataset._row_pos[int(r)] for r in np.asarray(rows)], dtype=np.int64)
        label = dataset._fps[dataset._fp_ix[pos]].copy()
        label *= dataset._dose_scale[pos, None]
        return label

    def meta(self, rows: np.ndarray) -> pd.DataFrame:
        obs = self.dataset.obs
        frame = obs.iloc[np.asarray(rows, dtype=np.int64)][
            ["plate", "cell_line_id", "dose_um"]
        ].copy()
        frame.insert(0, "row", np.asarray(rows, dtype=np.int64))
        if "observation_id" in obs.columns:
            frame.insert(1, "observation_id", obs["observation_id"].to_numpy()[np.asarray(rows)])
        frame["condition"] = self.dataset._condition[
            [self.dataset._row_pos[int(r)] for r in np.asarray(rows)]
        ]
        return frame.reset_index(drop=True)
