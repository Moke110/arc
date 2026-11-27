"""Pool native tenv observations and validate the resulting training release."""
from __future__ import annotations

from collections import defaultdict
from contextlib import ExitStack
import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
sys.path[:0] = [str(path) for path in (
    STEP / "output/api", STAGE / "util/output/api", STAGE / "01_data_raw/output/api",
    STAGE / "02_gene_axis/02.1_hvg-deg/output/api",
)]

from native_corpus import open_corpus
from prediction_axis import open_axes
from cpst_util.derivation import release
from cpst_util.publication import publication_lock, publish_candidate
from training_data import (TENV_FILES, TENV_FOLD_COLUMNS, TENV_MATRIX_FILES, TENV_SCHEMA,
                           TENV_STRIP_COLUMNS, sha256)

from fold_ledger import FOLD_RECEIPT, FOLD_ROOT, load_folds

from frozen_folds import group_labels

TREATMENT_KEY = ("plate", "cell_line", "drug_parent", "dose_um")
CONTROL_KEY = TREATMENT_KEY[:2]
TENV_RECIPE = {
    "normalization": "ordered-axis integer sum, CP10K on selected-axis total, natural log1p in float32",
    "C": "equal-weight mean of normalized M_C rows",
    "delta": "normalized M_T treatment minus matched normalized M_C",
}


@dataclass(frozen=True)
class PooledGroup:
    key: tuple
    pooled_row: int
    members: list[int]
    library_total: int

    @property
    def control_key(self):
        return self.key[:2]


def _treatment_reason(row):
    if row.view_exclusion_reason:
        return row.view_exclusion_reason
    if not isinstance(row.drug_parent_inchikey, str) or not row.drug_parent_inchikey:
        return "unresolved_drug"
    if pd.isna(row.dose_um) or not np.isfinite(row.dose_um) or row.dose_um <= 0:
        return "unresolved_dose"
    if not row.plate or not row.cell_line:
        return "unresolved_control_context"
    return ""


def _population(source, solvent):
    text_columns = ["view_exclusion_reason", "plate", "cell_line", "drug_parent_inchikey", "vehicle"]
    obs = source.obs[[*text_columns, "role", "dose_um"]].copy()
    obs[text_columns] = obs[text_columns].astype("string").fillna("")
    controls = defaultdict(list)
    treatments = defaultdict(list)
    exclusions = {}
    keys = {}
    for i, row in enumerate(obs.itertuples(index=False)):
        if row.role == "vehicle":
            if row.vehicle != solvent:
                exclusions[i] = "incompatible_vehicle"
            elif not row.plate or not row.cell_line:
                exclusions[i] = "unresolved_control_context"
            else:
                key = (row.plate, row.cell_line)
                controls[key].append(i)
                keys[i] = key
        elif row.role != "treated":
            exclusions[i] = row.view_exclusion_reason or f"role_{row.role}"
    for i, row in enumerate(obs.itertuples(index=False)):
        if row.role != "treated":
            continue
        reason = _treatment_reason(row)
        key = (row.plate, row.cell_line)
        if not reason and not any(keys.get(member) == key for member in source.matched_control_rows(i)):
            reason = "no_matched_vehicle"
        if reason:
            exclusions[i] = reason
        else:
            treatments[(*key, row.drug_parent_inchikey, row.dose_um)].append(i)
    return obs, controls, treatments, exclusions


def _pool_counts(path, projection, n_rows, groups):
    group_index = np.full(n_rows, -1, dtype=np.int64)
    for index, members in enumerate(groups):
        group_index[members] = index
    pooled = np.lib.format.open_memmap(
        path, mode="w+", dtype=np.int64, shape=(len(groups), len(projection.columns))
    )
    pooled[:] = 0
    for start, block in projection.batches():
        indices = group_index[start:start + len(block)]
        valid = indices >= 0
        np.add.at(pooled, indices[valid], block[valid])
    return pooled


def _select_groups(pooled, controls, treatments, exclusions):
    active_controls = {}
    active_treatments = []
    for group_index, (key, members) in enumerate(controls.items()):
        total = pooled[group_index].sum(dtype=np.int64)
        if total > 0:
            active_controls[key] = PooledGroup(key, group_index, members, total)
        else:
            exclusions.update(dict.fromkeys(members, "zero_control_counts"))
    for group_index, (key, members) in enumerate(treatments.items(), len(controls)):
        total = pooled[group_index].sum(dtype=np.int64)
        if key[:2] not in active_controls:
            exclusions.update(dict.fromkeys(members, "zero_control_counts"))
        elif total <= 0:
            exclusions.update(dict.fromkeys(members, "zero_treatment_counts"))
        else:
            active_treatments.append(PooledGroup(key, group_index, members, total))
    if not active_controls or not active_treatments:
        raise ValueError("training release requires usable controls and treatments")
    return active_controls, active_treatments


def _normalize(counts, totals):
    return np.log1p(counts.astype(np.float64) * (10000.0 / np.asarray(totals)[:, None])).astype(np.float32)


def _write_tables(root, obs, controls, treatments, exclusions, folds, control_rows):
    control_table = [(*group.key, len(group.members), group.library_total)
                     for group in controls.values()]
    controls_frame = pd.DataFrame(control_table, columns=[*CONTROL_KEY, "member_count", "library_total"])
    sizes = [len(group.members) for group in treatments]
    members = [member for group in treatments for member in group.members]
    groups = np.repeat(np.arange(len(treatments)), sizes)
    labels = group_labels(folds, obs.iloc[members], groups, strip_label="")
    treatment_table = [
        (*group.key, control_rows[group.control_key], size, group.library_total, *label)
        for group, size, label in zip(treatments, sizes, labels, strict=True)
    ]
    treatments_frame = pd.DataFrame(treatment_table, columns=[
        *TREATMENT_KEY, "control_row", "member_count", "library_total", *TENV_FOLD_COLUMNS, *TENV_STRIP_COLUMNS,
    ])
    controls_frame = controls_frame.drop(columns="cell_line")
    treatments_frame = treatments_frame.drop(columns="cell_line")
    if controls_frame.plate.duplicated().any():
        raise ValueError("tenv training requires one cell line per plate")
    controls_frame.to_csv(root / "controls.csv", index=False)
    treatments_frame.to_csv(root / "treatments.csv", index=False)
    with (root / "members.csv").open("w", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(("profile_id", "role", "row"))
        for role, groups in (("control", controls.values()), ("treatment", treatments)):
            for row, group in enumerate(groups):
                writer.writerows((obs.index[member], role, row) for member in group.members)
    with (root / "exclusions.csv").open("w", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(("profile_id", "role", "reason"))
        writer.writerows((obs.index[index], obs.role.iloc[index], reason)
                         for index, reason in sorted(exclusions.items()))
    return controls_frame, treatments_frame


def _write_row_metadata(handle, table):
    for column, series in table.items():
        if pd.api.types.is_string_dtype(series.dtype):
            handle.create_dataset(column, data=series.to_numpy(dtype=object), dtype=h5py.string_dtype())
        else:
            handle.create_dataset(column, data=series.to_numpy())


def _write_matrices(root, pooled, controls, treatments, control_rows, tables):
    width = pooled.shape[1]
    control_groups = list(controls.values())
    matrix = np.empty((len(controls), width), dtype=np.float32)
    for start in range(0, len(controls), 64):
        groups = control_groups[start:start + 64]
        rows = [group.pooled_row for group in groups]
        totals = [group.library_total for group in groups]
        matrix[start:start + len(groups)] = _normalize(pooled[rows], totals)
    with ExitStack() as stack:
        files = TENV_MATRIX_FILES
        handles = {file: stack.enter_context(h5py.File(root / file, "w"))
                   for file in dict.fromkeys(files.values())}
        reference = handles[files["C"]]
        plate = handles[files["M_C"]]
        perturbation = handles[files["M_T/treatment"]]
        reference.create_dataset("C", data=matrix.mean(axis=0, dtype=np.float64).astype(np.float32))
        plate.create_dataset("M_C", data=matrix)
        _write_row_metadata(plate, tables[0])
        _write_row_metadata(perturbation, tables[1])
        shape = (len(treatments), width)
        chunk = (min(64, shape[0]), width)
        treatment = perturbation.create_dataset("M_T/treatment", shape, np.float32,
                                                chunks=chunk, compression="gzip", compression_opts=1)
        delta = perturbation.create_dataset("M_T/delta", shape, np.float32,
                                            chunks=chunk, compression="gzip", compression_opts=1)
        for start in range(0, len(treatments), 64):
            block = treatments[start:start + 64]
            rows = [group.pooled_row for group in block]
            totals = [group.library_total for group in block]
            values = _normalize(pooled[rows], totals)
            matched = matrix[[control_rows[group.control_key] for group in block]]
            treatment[start:start + len(block)] = values
            delta[start:start + len(block)] = values - matched


def _build_dataset(root, corpus, axis, fold_axis, fold_root, fold_receipt):
    folds = load_folds(fold_root, fold_receipt, "tenv", fold_axis.sha256)
    projection = axis.project(corpus, "tenv")
    obs, controls, treatments, exclusions = _population(projection.source, "DMSO")
    if not controls or not treatments:
        raise ValueError("tenv has no candidate controls or treatments")
    root.mkdir(exist_ok=True)
    cache = root / "pooled.npy"
    pooled = _pool_counts(cache, projection, len(obs), [*controls.values(), *treatments.values()])
    try:
        active_controls, active_treatments = _select_groups(pooled, controls, treatments, exclusions)
        control_rows = {key: row for row, key in enumerate(active_controls)}
        tables = _write_tables(root, obs, active_controls, active_treatments, exclusions, folds,
                               control_rows)
        _write_matrices(root, pooled, active_controls, active_treatments, control_rows,
                        tables)
    finally:
        del pooled
        cache.unlink()
    pd.DataFrame({"gene_index": axis.gene_indices, "symbol": axis.symbols}).to_csv(
        root / "genes.csv", index=False
    )
    source_record = corpus.manifest["sources"]["tenv"]
    release.write_json(root / "manifest.json", {
        "schema": TENV_SCHEMA,
        "source": "tenv",
        "axis_size": len(axis.gene_indices),
        "recipe": {**TENV_RECIPE, "vehicle": "DMSO"},
        "inputs": {
            "raw_source": "tenv",
            "raw_sha256": corpus.manifest_sha256,
            "raw_dataset_sha256": source_record["manifest_sha256"],
            "raw_counts_sha256": source_record["sha256"],
            "axis_sha256": axis.sha256,
            "folds": {**folds.evidence, "axis_sha256": fold_axis.sha256},
        },
        "counts": {
            "raw_rows": len(obs), "candidate_controls": len(controls),
            "candidate_treatments": len(treatments), "controls": len(active_controls),
            "treatments": len(active_treatments), "excluded_raw_rows": len(exclusions),
        },
        "outputs": {name: sha256(root / name) for name in TENV_FILES},
    })


def _verify_dataset(data, obs, axis):
    if (tuple(data.genes.gene_index) != axis.gene_indices
            or tuple(data.genes.symbol) != axis.symbols):
        raise ValueError("training gene order differs")
    members = data.membership()
    excluded = data.exclusions()
    identities = pd.concat([members.profile_id, excluded.profile_id], ignore_index=True)
    if (identities.duplicated().any() or set(identities) != set(obs.index)
            or not members.role.isin({"control", "treatment"}).all()
            or excluded.reason.eq("").any()
            or len(excluded) != data.manifest["counts"]["excluded_raw_rows"]):
        raise ValueError("training membership and exclusions do not partition the raw observations")
    for role, table in (("control", data.controls), ("treatment", data.treatments)):
        counts = members.loc[members.role.eq(role), "row"].value_counts().sort_index()
        if (list(counts.index) != list(table.index)
                or not np.array_equal(counts.to_numpy(), table.member_count.to_numpy())):
            raise ValueError(f"{role} member counts differ from the ledger")
    control_sum = np.zeros(len(data.genes), dtype=np.float64)
    controls = []
    for start in range(0, len(data.controls), 64):
        block = data.read_controls(start, min(start + 64, len(data.controls)))
        if not np.isfinite(block).all() or (block < 0).any():
            raise ValueError("normalized control values must be finite and nonnegative")
        control_sum += block.sum(axis=0, dtype=np.float64)
        controls.append(block)
    if not np.array_equal(data.control(), (control_sum / len(data.controls)).astype(np.float32)):
        raise ValueError("C is not the equal-weight mean of M_C")
    control_matrix = np.concatenate(controls)
    for start in range(0, len(data.treatments), 64):
        stop = min(start + 64, len(data.treatments))
        treatment = data.read_treatments(start, stop)
        control = control_matrix[data.treatments.control_row.iloc[start:stop].to_numpy()]
        delta = data.read_treatments(start, stop, "delta")
        if (not np.isfinite(treatment).all() or (treatment < 0).any()
                or not np.array_equal(delta, treatment - control)):
            raise ValueError("M_T contains invalid expression or unmatched Delta")
