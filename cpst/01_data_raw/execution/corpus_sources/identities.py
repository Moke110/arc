"""Assign corpus observation roles, control groups and parent identities."""
from functools import lru_cache

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from cpst_util.derivation.release import ReleaseError


class CorpusError(ReleaseError):
    """Observation metadata violates the source contract."""


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
