"""Integer raw corpus on the frozen gene axis, with matched controls as rows (ADR 0066).

The corpus is the model-facing raw truth: integer counts on the frozen 18,164
``gene_index`` axis, with **controls as ordinary rows** rather than a separate
layer, and no normalization, no delta, no HVG subsetting and no filter. Row sets
are ``tenv`` 4,206 treated + 272 vehicle wells and ``tahoe_bulk`` 65,318 treated
+ 1,400 ``DMSO_TF`` rows.

Every row carries its ``role`` (``treated``/``vehicle``), source profile/well
identity, plate, cell line, desalted drug-parent structure, micromolar and
original dose, the source item and its item-manifest digest, and its
``matched_control_group``. A row without a compatible control is excluded and
its reason is recorded in the release exclusions.

Matched-control correspondence (ADR 0061 semantics retained):

* ``tenv`` — same plate; the control is the same-plate DMSO vehicle. The 80
  Media wells are vehicle rows but are not matched to any treated row.
* ``tahoe_bulk`` — same plate **and** cell line; the control is ``DMSO_TF``.

The public reader (:class:`Corpus`) is the only supported access path; it exposes
bounded dense ``int32`` row windows and the group->control-row resolution that
the view builder consumes.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from . import raw_counts
from .release import (ReleaseError, check_manifest, output_digests, publish,
                      read_json, replay, sha256, write_json)

SCHEMA = "stage2-bases-corpus-v1"
AXIS_SIZE = 18164
BATCH = 256
FAM153B_INDEX = 5091
FAM153B_HGNC = "HGNC:27323"
SOURCES = ("tenv", "tahoe_bulk")
ROLE_VALUES = ("treated", "vehicle")

EXPECTED = {
    "tenv": {"profiles": 4478, "treated": 4206, "vehicle": 272,
             "excluded": {"technical_blank": 32},
             "matched_control_rows": 192, "matched_groups": 16,
             "vehicle_groups": 32},
    "tahoe_bulk": {"profiles": 66718, "treated": 65318, "vehicle": 1400,
                   "excluded": {"missing_or_invalid_structure": 300},
                   "matched_control_rows": 1400, "matched_groups": 700,
                   "vehicle_groups": 700},
}

RECIPE = {
    "tenv": ("Project the native tenv raw counts through tenv_axis_mapping.csv "
             "onto the frozen 18,164 axis; keep the 4,206 treated wells and the "
             "272 vehicle wells as rows; desalted parent structure and "
             "micromolar dose from the native release; same-plate DMSO is the "
             "matched control. int32; no normalization, no delta."),
    "tahoe_bulk": ("Take the projected Tahoe raw counts (67,018 profiles) and "
                   "keep the 65,318 treated rows with a resolvable desalted "
                   "parent structure plus the 1,400 DMSO_TF rows as vehicle "
                   "rows; matched control is the same-plate, same-cell-line "
                   "DMSO_TF row. int32; no normalization, no delta."),
    "controls": ("Controls are rows, not layers. matched_control_group is a "
                 "stable string; the reader resolves the group to its vehicle "
                 "rows. No compatible control means the row is excluded with a "
                 "recorded reason."),
}

PRODUCT_FILES = ("tenv.h5ad", "tahoe_bulk.h5ad", "exclusions.csv")


class CorpusError(ReleaseError):
    """A corpus release violates its contract."""


# --------------------------------------------------------------------------- #
# shared structure/identity helpers (ADR 0061 desalting recipe)
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=4096)
def parent_structure(smiles):
    """Desalted parent InChIKey + isomeric SMILES from one source SMILES."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError("invalid_structure")
    exact = Chem.MolToInchiKey(mol)
    parent = rdMolStandardize.Uncharger().uncharge(
        rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol))
    )
    Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
    key = Chem.MolToInchiKey(parent)
    if len(key) != 27:
        raise ValueError("unresolved_parent_structure")
    return key, Chem.MolToSmiles(parent, isomericSmiles=True), exact


def dose_micromolar(value, unit="uM"):
    factors = {"uM": 1.0, "µM": 1.0, "μM": 1.0, "nM": 0.001, "mM": 1000.0, "M": 1e6}
    if unit not in factors:
        raise ValueError("unresolvable_dose_unit")
    dose = float(value) * factors[unit]
    if not np.isfinite(dose) or dose <= 0:
        raise ValueError("invalid_treated_dose")
    return dose


def tahoe_structures(drug_metadata):
    """Resolve every source drug name to its desalted parent (deterministic order)."""
    frame = pd.read_parquet(drug_metadata)
    if "drug" not in frame or "canonical_smiles" not in frame:
        raise CorpusError("Tahoe drug metadata lacks drug/canonical_smiles")
    structures = {}
    for name, part in frame.dropna(subset=["drug"]).groupby("drug", sort=True):
        values = sorted({str(value).strip() for value in part.canonical_smiles.dropna()
                         if str(value).strip()})
        resolved = []
        for value in values:
            try:
                resolved.append((value, *parent_structure(value)))
            except ValueError:
                continue
        if len({item[1] for item in resolved}) > 1:
            raise CorpusError(f"ambiguous Tahoe parent structure: {name}")
        structures[str(name).strip()] = resolved
    return structures


def _structure_row(value):
    try:
        key, smiles, exact = parent_structure(value)
    except ValueError:
        return None
    return {"drug_parent_inchikey": key, "drug_parent_smiles": smiles,
            "source_inchikey": exact}


# --------------------------------------------------------------------------- #
# row metadata and matched-control correspondence
# --------------------------------------------------------------------------- #
def _tenv_rows(obs, item):
    records, exclusions = [], []
    for source_row, row in enumerate(obs.itertuples(index=False)):
        profile = str(obs.index[source_row])
        if bool(row.is_blank):
            exclusions.append({"source_row": source_row, "profile_id": profile,
                               "reason": "technical_blank"})
            continue
        vehicle = bool(row.is_control)
        record = {
            "profile_id": profile,
            "source_item": item["name"],
            "source_item_manifest_sha256": item["sha256"],
            "source_row": source_row,
            "role": "vehicle" if vehicle else "treated",
            "plate": str(row.plate),
            "cell_line": "hiPSC-CM",
            "cell_line_id": "",
            "compound_name": str(row.compound_name),
            "vehicle": str(row.vehicle),
            "dose_um": float(row.dose_um),
            "dose_original_value": float(row.dose_um),
            "dose_original_unit": "uM",
            "drug_parent_inchikey": "",
            "drug_parent_smiles": "",
            "source_inchikey": "",
            "n_cells": 0,
            "matched_control_group": f"tenv::{row.plate}::{row.vehicle}",
        }
        if vehicle:
            record["role"] = "vehicle"
        else:
            record["matched_control_group"] = f"tenv::{row.plate}::DMSO"
            if str(row.vehicle) != "DMSO":
                exclusions.append({"source_row": source_row, "profile_id": profile,
                                   "reason": "incompatible_matched_control_conditions"})
                continue
            try:
                dose_micromolar(row.dose_um)
            except (TypeError, ValueError) as error:
                exclusions.append({"source_row": source_row, "profile_id": profile,
                                   "reason": str(error)})
                continue
            structure = _structure_row(row.smiles)
            if structure is None:
                exclusions.append({"source_row": source_row, "profile_id": profile,
                                   "reason": "missing_or_invalid_structure"})
                continue
            record.update(structure)
        records.append(record)
    return records, exclusions


def _tahoe_rows(obs, item, structures):
    records, exclusions = [], []
    for source_row, row in enumerate(obs.itertuples(index=False)):
        profile = str(obs.index[source_row])
        vehicle = bool(row.is_control)
        record = {
            "profile_id": profile,
            "source_item": item["name"],
            "source_item_manifest_sha256": item["sha256"],
            "source_row": source_row,
            "role": "vehicle" if vehicle else "treated",
            "plate": str(row.plate),
            "cell_line": str(row.cell_name),
            "cell_line_id": str(row.cell_line_id),
            "compound_name": str(row.drug),
            "vehicle": "",
            "dose_um": float(row.dose_um),
            "dose_original_value": float(row.dose_um),
            "dose_original_unit": "uM",
            "drug_parent_inchikey": "",
            "drug_parent_smiles": "",
            "source_inchikey": "",
            "n_cells": int(row.n_cells),
            "matched_control_group": f"tahoe::{row.plate}::{row.cell_line_id}",
        }
        if vehicle:
            record["vehicle"] = "DMSO_TF"
        else:
            candidates = structures.get(str(row.drug).strip(), [])
            if not candidates:
                exclusions.append({"source_row": source_row, "profile_id": profile,
                                   "reason": "missing_or_invalid_structure"})
                continue
            _smiles, key, parent_smiles, exact = candidates[0]
            record.update({"drug_parent_inchikey": key,
                           "drug_parent_smiles": parent_smiles,
                           "source_inchikey": exact})
            try:
                dose_micromolar(row.dose_um)
            except (TypeError, ValueError) as error:
                exclusions.append({"source_row": source_row, "profile_id": profile,
                                   "reason": str(error)})
                continue
        records.append(record)
    return records, exclusions


def _control_report(frame, source):
    """Resolve groups to vehicle rows; drop treated rows with no compatible control."""
    vehicle = frame.role.eq("vehicle")
    groups = {}
    for index in np.flatnonzero(vehicle.to_numpy()):
        groups.setdefault(str(frame.matched_control_group.iloc[index]), []).append(int(index))
    keep = np.ones(len(frame), dtype=bool)
    exclusions = []
    if (~vehicle).any():
        treated = np.flatnonzero((~vehicle).to_numpy())
        for index in treated:
            group = str(frame.matched_control_group.iloc[index])
            if not groups.get(group):
                keep[index] = False
                exclusions.append({"source_row": int(frame.source_row.iloc[index]),
                                   "profile_id": str(frame.profile_id.iloc[index]),
                                   "reason": "missing_matched_control"})
    matched_groups = {str(frame.matched_control_group.iloc[index])
                      for index in np.flatnonzero((~vehicle.to_numpy()) & keep)}
    matched_rows = sum(len(groups.get(group, ())) for group in matched_groups)
    report = {
        "rule": ("same plate; same-plate DMSO vehicle" if source == "tenv"
                 else "same plate and cell line; DMSO_TF vehicle"),
        "vehicle_rows": int(vehicle.sum()),
        "vehicle_groups": int(len(groups)),
        "matched_control_rows": int(matched_rows),
        "matched_groups": int(len(matched_groups)),
        "treated_without_control": int((~keep & ~vehicle.to_numpy()).sum()),
    }
    return keep, report, exclusions


def _write_h5ad(path, obs, var, source, raw, keep, projection):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        write_elem(handle, "obs", obs)
        write_elem(handle, "var", var)
        write_elem(handle, "uns", {
            "schema": SCHEMA, "source": source,
            "counts": "raw integer counts on the frozen 18,164 axis; no normalization, no delta",
            "controls": "vehicle rows are ordinary rows; matched_control_group resolves them",
        })
        write_elem(handle, "layers", {})
        matrix = handle.create_dataset(
            "X", shape=(len(obs), len(var)), dtype="int32",
            chunks=(min(64, len(obs)), len(var)), compression="gzip",
            compression_opts=1, shuffle=True)
        matrix.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
        position = 0
        for start, block in raw.batches(BATCH):
            local = np.flatnonzero(keep[start:start + len(block)])
            if not len(local):
                continue
            selected = block[local]
            if projection is not None:
                selected = selected[:, projection]
            matrix[position:position + len(selected)] = selected
            position += len(selected)
        if position != len(obs):
            raise CorpusError(f"{source}: corpus row accounting differs")
    temporary.replace(path)


def _var_frame(axis):
    var = axis.copy()
    var.index = pd.Index([str(value) for value in var.gene_index],
                         name="gene_index_label")
    return var


# --------------------------------------------------------------------------- #
# publication
# --------------------------------------------------------------------------- #
def _relative(path, repo_root):
    path = Path(path)
    if repo_root is None:
        return path.name
    root = Path(repo_root).absolute()
    try:
        return str(path.absolute().relative_to(root))
    except ValueError as error:
        raise CorpusError(f"corpus input outside the repository: {path}") from error


def _item(name, raw_manifest, repo_root):
    declared = raw_manifest["sources"][name if name == "tenv" else "tahoe"].get(
        "item_manifest")
    if not declared:
        if repo_root is not None:
            raise CorpusError(f"{name}: raw release has no item-manifest evidence")
        return {"name": "fixture", "sha256": ""}
    digest = declared["sha256"]
    if repo_root is not None:
        actual = sha256(Path(repo_root) / declared["path"])
        if actual != digest:
            raise CorpusError(f"{name}: source item manifest digest differs")
    return {"name": "GSE262419" if name == "tenv" else "Tahoe-100M",
            "sha256": digest, "path": declared["path"]}


def build(candidate, *, raw_counts_dir, axis_path, drug_metadata, repo_root=None,
          fixture=False):
    """Write a complete corpus candidate outside ``output/``."""
    candidate = Path(candidate)
    candidate.mkdir(parents=True, exist_ok=True)
    raw_root = Path(raw_counts_dir)
    raw_manifest = read_json(raw_root / "manifest.json")
    if raw_manifest.get("schema") != raw_counts.SCHEMA:
        raise CorpusError("unexpected raw-count release schema")
    if sha256(axis_path) != raw_manifest["axis"]["sha256"]:
        raise CorpusError("axis digest differs from the raw-count release")
    axis = pd.read_csv(axis_path, dtype={"symbol": str, "hgnc_id": str,
                                         "ensembl_id": str}, keep_default_na=False)
    if not fixture and len(axis) != AXIS_SIZE:
        raise CorpusError(f"expected {AXIS_SIZE} frozen genes, got {len(axis)}")
    axis = axis.reset_index(drop=True)
    axis.insert(0, "gene_index", np.arange(len(axis), dtype=np.int64))
    var = _var_frame(axis)
    if len(axis) > FAM153B_INDEX and str(axis.iloc[FAM153B_INDEX].hgnc_id) != FAM153B_HGNC:
        raise CorpusError("frozen FAM153B gene_index moved")
    inputs = {}
    for path in (raw_root / "manifest.json", axis_path):
        _add_input(path, inputs, repo_root)
    names = ["tenv.h5ad", "tenv_probe_mapping.csv", "tenv_axis_mapping.csv",
             "tahoe.h5ad", "token_mapping.csv"]
    for name in names:
        _add_input(raw_root / name, inputs, repo_root)
    _add_input(drug_metadata, inputs, repo_root)
    exclusions = []
    summaries = {}
    controls = {}
    for source in SOURCES:
        raw_name = "tahoe" if source == "tahoe_bulk" else "tenv"
        item = _item(raw_name, raw_manifest, repo_root)
        opener = (raw_counts.open_tenv if raw_name == "tenv"
                  else raw_counts.open_tahoe)
        with opener(raw_root) as raw:
            obs = raw.obs
            if raw_name == "tenv":
                records, dropped = _tenv_rows(obs, item)
                bridge = raw.bridge()
                projection = bridge.native_column.to_numpy()
            else:
                records, dropped = _tahoe_rows(
                    obs, item, tahoe_structures(drug_metadata))
                projection = None
            frame = pd.DataFrame(records)
            keep, report, unmatched = _control_report(frame, source)
            exclusions.extend({"source": source, **entry} for entry in dropped + unmatched)
            frame = frame.loc[keep].reset_index(drop=True)
            # keep is over raw payload rows (blanks/structure rows already absent
            # from `frame`, so recompute a per-raw-row mask from source_row).
            mask = np.zeros(len(obs), dtype=bool)
            mask[frame.source_row.to_numpy()] = True
            obs_frame = frame.set_index("profile_id", drop=False)
            obs_frame.index.name = "profile_id"
            _write_h5ad(candidate / f"{source}.h5ad", obs_frame, var, source,
                        raw, mask, projection)
            summaries[source] = {
                "counts": {"profiles": int(len(frame)),
                           "treated": int((frame.role == "treated").sum()),
                           "vehicle": int((frame.role == "vehicle").sum()),
                           "raw_profiles": int(len(obs))},
                "excluded_by_reason": _reason_counts(dropped + unmatched),
                "controls": report,
                "output": f"{source}.h5ad",
            }
            controls[source] = report
    exclusion_frame = pd.DataFrame(exclusions,
                                   columns=["source", "source_row", "profile_id",
                                            "reason"])
    exclusion_frame.to_csv(candidate / "exclusions.csv", index=False, lineterminator="\n")
    if not fixture:
        for source in SOURCES:
            expected = EXPECTED[source]
            summary = summaries[source]
            if (summary["counts"]["profiles"], summary["counts"]["treated"],
                    summary["counts"]["vehicle"]) != (
                    expected["profiles"], expected["treated"], expected["vehicle"]):
                raise CorpusError(f"{source}: unexpected corpus population: {summary}")
            if summary["excluded_by_reason"] != expected["excluded"]:
                raise CorpusError(f"{source}: unexpected exclusion accounting: {summary}")
            if controls[source]["matched_control_rows"] != expected["matched_control_rows"]:
                raise CorpusError(f"{source}: unexpected matched control count")
            if controls[source]["matched_groups"] != expected["matched_groups"]:
                raise CorpusError(f"{source}: unexpected matched control group count")
            if controls[source]["treated_without_control"]:
                raise CorpusError(f"{source}: retained treated rows lack a control")
    outputs = {}
    for name in PRODUCT_FILES:
        path = candidate / name
        record = {"sha256": sha256(path), "bytes": path.stat().st_size}
        if name.endswith(".h5ad"):
            with h5py.File(path) as handle:
                record.update({"rows": int(handle["X"].shape[0]),
                               "columns": int(handle["X"].shape[1]),
                               "dtype": str(handle["X"].dtype)})
        else:
            record.update({"rows": int(len(exclusion_frame))})
        outputs[name] = record
    manifest = {
        "schema": SCHEMA,
        "recipe": RECIPE,
        "fixture": bool(fixture),
        "axis": {"size": int(len(axis)),
                 "path": (raw_manifest["axis"]["path"] if not fixture else _relative(axis_path, repo_root)),
                 "sha256": raw_manifest["axis"]["sha256"]},
        "inputs": inputs,
        "outputs": outputs,
        "sources": summaries,
        "controls": {source: controls[source] for source in SOURCES},
        "counts": {"axis_size": int(len(axis)),
                   "tenv_profiles": summaries["tenv"]["counts"]["profiles"],
                   "tahoe_bulk_profiles": summaries["tahoe_bulk"]["counts"]["profiles"]},
        "rdkit_version": rdBase.rdkitVersion,
    }
    write_json(candidate / "manifest.json", manifest)
    return manifest


def _add_input(path, inputs, repo_root, key=None):
    path = Path(path)
    key = _relative(path, repo_root) if key is None else key
    if key in inputs:
        raise CorpusError(f"duplicate corpus input path: {key}")
    inputs[key] = {"sha256": sha256(path), "bytes": path.stat().st_size}


def _reason_counts(entries):
    counts = {}
    for entry in entries:
        counts[entry["reason"]] = counts.get(entry["reason"], 0) + 1
    return {str(key): int(value) for key, value in sorted(counts.items())}


def verify(candidate, fixture=False):
    """Re-check a corpus candidate's schema, roles, controls and index space."""
    candidate = Path(candidate)
    manifest_path = candidate / "manifest.json"
    if not manifest_path.is_file():
        raise CorpusError("candidate has no manifest.json")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != SCHEMA:
        raise CorpusError("unexpected corpus schema")
    if bool(manifest.get("fixture")) is not bool(fixture):
        raise CorpusError("corpus fixture flag differs from the verification mode")
    check_manifest(candidate, manifest)
    axis_size = int(manifest["counts"]["axis_size"])
    expected = None if fixture else EXPECTED
    exclusions = pd.read_csv(candidate / "exclusions.csv", keep_default_na=False)
    if not {"source", "source_row", "profile_id", "reason"} <= set(exclusions.columns):
        raise CorpusError("exclusions.csv lacks its reason columns")
    for source in SOURCES:
        summary = manifest["sources"][source]
        path = candidate / summary["output"]
        with h5py.File(path, "r") as handle:
            obs = read_elem(handle["obs"])
            var = read_elem(handle["var"])
            if set(handle["layers"]):
                raise CorpusError(f"{source}: corpus must not carry normalized layers")
            if handle["X"].dtype.kind not in "iu":
                raise CorpusError(f"{source}: corpus counts must be integer")
            if tuple(handle["X"].shape) != (len(obs), len(var)):
                raise CorpusError(f"{source}: count shape differs from metadata")
            if handle["X"][:min(8, len(obs))].min() < 0:
                raise CorpusError(f"{source}: negative corpus counts")
        if len(var) != axis_size or var.gene_index.tolist() != list(range(axis_size)):
            raise CorpusError(f"{source}: var is not the ordered frozen axis")
        obs_columns = set(obs.columns)
        for forbidden in ("control", "delta", "delta_mean", "normalized"):
            if forbidden in obs_columns:
                raise CorpusError(f"{source}: corpus obs carries a {forbidden} column")
        if not set(obs.role) <= set(ROLE_VALUES):
            raise CorpusError(f"{source}: invalid role values")
        if obs.role.eq("vehicle").sum() != summary["counts"]["vehicle"]:
            raise CorpusError(f"{source}: vehicle row count differs from the manifest")
        if obs.role.eq("treated").sum() != summary["counts"]["treated"]:
            raise CorpusError(f"{source}: treated row count differs from the manifest")
        if obs.index.has_duplicates or obs.index.hasnans:
            raise CorpusError(f"{source}: duplicate or missing profile identities")
        if obs.matched_control_group.eq("").any():
            raise CorpusError(f"{source}: row without a matched_control_group")
        groups = {}
        for position in np.flatnonzero(obs.role.eq("vehicle").to_numpy()):
            key = str(obs.matched_control_group.iloc[position])
            groups.setdefault(key, []).append(int(position))
        treated = np.flatnonzero(obs.role.eq("treated").to_numpy())
        for position in treated:
            key = str(obs.matched_control_group.iloc[position])
            if not groups.get(key):
                raise CorpusError(
                    f"{source}: treated row {position} has no vehicle control row")
        for field in ("source_item", "source_item_manifest_sha256",
                      "drug_parent_inchikey", "dose_um", "dose_original_value"):
            if field not in obs_columns:
                raise CorpusError(f"{source}: missing corpus field {field}")
        if obs.loc[obs.role.eq("treated"), "drug_parent_inchikey"].eq("").any():
            raise CorpusError(f"{source}: treated row lacks a desalted parent structure")
        if expected is not None:
            if obs.role.eq("treated").sum() != expected[source]["treated"]:
                raise CorpusError(f"{source}: unexpected treated total")
            if obs.role.eq("vehicle").sum() != expected[source]["vehicle"]:
                raise CorpusError(f"{source}: unexpected vehicle total")
        counts = exclusions[exclusions.source.eq(source)].reason.value_counts().to_dict()
        if {str(k): int(v) for k, v in counts.items()} != summary["excluded_by_reason"]:
            raise CorpusError(f"{source}: exclusion accounting differs")
    if len(exclusions) != sum(
            sum(summary["excluded_by_reason"].values())
            for summary in manifest["sources"].values()):
        raise CorpusError("exclusions.csv row count differs from the manifest")
    return {"result": "PASS", "manifest_sha256": sha256(manifest_path),
            "outputs": sorted(output_digests(manifest))}


def publish_release(output, *, root=None, **kwargs):
    return publish(output, lambda stage: build(stage, **kwargs),
                   lambda stage: verify(stage, fixture=kwargs.get("fixture", False)),
                   root=root)


def replay_release(output, *, root=None, **kwargs):
    return replay(output, lambda stage: build(stage, **kwargs),
                  lambda stage: verify(stage, fixture=kwargs.get("fixture", False)),
                  root=root)


# --------------------------------------------------------------------------- #
# public reader
# --------------------------------------------------------------------------- #
class CorpusSource:
    """One row set of the corpus with bounded count access and control resolution."""

    def __init__(self, root, source):
        self.root = Path(root)
        self.source = source
        self.path = self.root / f"{source}.h5ad"
        self._handle = h5py.File(self.path, "r")
        self.obs = read_elem(self._handle["obs"])
        self.var = read_elem(self._handle["var"])
        self._matrix = self._handle["X"]
        self.shape = (len(self.obs), len(self.var))
        if list(self.var.gene_index) != list(range(self.shape[1])):
            raise CorpusError(f"{source}: var is not the ordered frozen axis")
        self._groups = {}
        for position in np.flatnonzero(self.obs.role.eq("vehicle").to_numpy()):
            key = str(self.obs.matched_control_group.iloc[position])
            self._groups.setdefault(key, []).append(int(position))

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def gene_index(self):
        """Frozen gene_index vector; the corpus is never renumbered."""
        return np.asarray(self.var.gene_index, dtype=np.int64)

    def counts(self, start=0, stop=None):
        stop = self.shape[0] if stop is None else int(stop)
        return np.asarray(self._matrix[start:stop], dtype=np.int32)

    def batches(self, batch_size=BATCH):
        if batch_size <= 0:
            raise CorpusError("batch_size must be positive")
        for start in range(0, self.shape[0], batch_size):
            stop = min(start + batch_size, self.shape[0])
            yield start, self.counts(start, stop)

    def control_groups(self):
        return {key: np.asarray(rows, dtype=np.int64)
                for key, rows in self._groups.items()}

    def control_rows(self, group):
        if group not in self._groups:
            raise CorpusError(f"unknown matched control group: {group}")
        return np.asarray(self._groups[group], dtype=np.int64)

    def matched_control_rows(self, index):
        """Vehicle row positions matched to corpus row ``index``."""
        if not 0 <= int(index) < self.shape[0]:
            raise CorpusError(f"corpus row out of range: {index}")
        group = str(self.obs.matched_control_group.iloc[int(index)])
        return self.control_rows(group)


class Corpus:
    """Read-only corpus release: two row sets on the frozen axis."""

    def __init__(self, root):
        self.root = Path(root)
        self.manifest = read_json(self.root / "manifest.json")
        if self.manifest.get("schema") != SCHEMA:
            raise CorpusError("unsupported corpus schema")
        self.sources = tuple(sorted(self.manifest["sources"]))
        if set(self.sources) != set(SOURCES):
            raise CorpusError(f"unexpected corpus sources: {self.sources}")
        self.axis_size = int(self.manifest["counts"]["axis_size"])

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def source(self, name):
        if name not in SOURCES:
            raise CorpusError(f"unknown corpus source: {name}")
        return CorpusSource(self.root, name)


def open_corpus(root):
    return Corpus(root)
