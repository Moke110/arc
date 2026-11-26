"""Pool raw tenv counts on the common CM axis for shared disease direction truth."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from prediction_axis import COMMON_SHA256, COMMON_SIZE, sha256


def disease_gene_sets(axis_release):
    """Load all eligible published DEG1000 lists with their release digests."""
    root = Path(axis_release)
    manifest = json.loads((root / "manifest.json").read_bytes())
    prefix = "evidence/disease-fusion/"

    def admitted(name):
        path = root / (prefix + name)
        if Path(name).name != name or sha256(path) != manifest["outputs"][prefix + name]["sha256"]:
            raise ValueError(f"disease DEG1000 evidence digest differs: {name}")
        return path

    if manifest["common_axis"] != {"size": COMMON_SIZE, "sha256": COMMON_SHA256}:
        raise ValueError("disease reference must use the frozen CM axis")
    catalog = admitted("catalog.json")
    disease = {}
    pins = {"catalog.json": sha256(catalog)}
    for row in json.loads(catalog.read_bytes())["diseases"]:
        if row["eligible"]:
            name = row["category"]
            path = admitted(row["deg1000_file"])
            frame = pd.read_csv(path)
            genes = frame.gene_index.to_numpy()
            if (name in disease or len(genes) != row["deg1000"] or len(genes) not in range(1, 1001)
                    or frame.disease_rank.tolist() != list(range(1, len(genes) + 1))
                    or genes.dtype.kind not in "iu" or frame.gene_index.duplicated().any()
                    or np.any((genes < 0) | (genes >= COMMON_SIZE))):
                raise ValueError(f"invalid released disease DEG1000 list: {name}")
            disease[name] = genes
            pins[path.name] = sha256(path)
    if not disease:
        raise ValueError("no eligible disease DEG1000 lists")
    return disease, pins


def _normalize_counts(counts):
    totals = counts.sum(axis=1, dtype=np.int64)
    if (totals <= 0).any():
        raise ValueError("pooled common-axis counts must have positive library totals")
    return np.log1p(counts.astype(np.float64) * (10000 / totals[:, None])).astype(np.float32)


def write_cm_directions(source, release, path, *, common_size=COMMON_SIZE):
    """Write one signed measured Delta row per baseline training condition."""
    if not source.shared_axis_complete:
        raise ValueError("tenv must cover the complete common CM gene axis")
    columns = source.native_columns(range(common_size))
    controls, treatments = len(release.controls), len(release.treatments)
    members = release.membership()
    member_rows = source.obs.index.get_indexer(members.profile_id)
    roles = members.role
    rows = members.row.to_numpy()
    control_member = roles.eq("control").to_numpy()
    treatment_member = roles.eq("treatment").to_numpy()
    if (members.profile_id.duplicated().any() or (member_rows < 0).any() or (rows < 0).any()
            or (rows[control_member] >= controls).any() or (rows[treatment_member] >= treatments).any()
            or not np.all(control_member | treatment_member)):
        raise ValueError("common-axis direction rows require valid raw memberships")
    group_indices = np.full(len(source.obs), -1, np.int64)
    group_indices[member_rows] = rows + controls * treatment_member
    path = Path(path)
    pooled_path = path.with_suffix(".pool.npy")
    pooled = np.lib.format.open_memmap(pooled_path, mode="w+", dtype=np.int64,
                                      shape=(controls + treatments, common_size))
    try:
        pooled[:] = 0
        for start, block in source.batches(64):
            groups = group_indices[start:start + len(block)]
            selected = groups >= 0
            np.add.at(pooled, groups[selected], block[selected][:, columns])
        reference = _normalize_counts(pooled[:controls])
        matches = release.treatments.control_row.to_numpy()
        directions = np.lib.format.open_memmap(path, mode="w+", dtype=np.int8,
                                               shape=(treatments, common_size))
        for start in range(0, treatments, 64):
            stop = min(start + 64, treatments)
            measured = _normalize_counts(pooled[controls + start:controls + stop])
            directions[start:stop] = np.sign(measured - reference[matches[start:stop]])
        directions.flush()
    finally:
        del pooled
        pooled_path.unlink()
