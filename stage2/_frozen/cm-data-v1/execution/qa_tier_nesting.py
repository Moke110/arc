# FROZEN COPY — archived evidence, NOT executable code and NOT a step entry point.
# Origin: stage2/02_cm_data/execution/cmdata/qa.py (cross-view tier-nesting gates) @ commit 7e7d689 (2026-11-04).
# See ./README.md. Do not wire this into any build; do not edit it.
"""Frozen v1 cross-view QA gates (archival excerpt).

Extracted from the retired ``stage2/02_cm_data/execution/cmdata/qa.py``: the
gates that assert one shared var axis, the temp_env tier nesting
(``atc-moa ⊆ atc-moa-mesh ⊆ all-smiles``) and the shared temp_env controls.
The original module's source spot checks, Stage 1 regression and report writer
are not frozen here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from paths_v1 import TEMP_ENV_VIEWS

if TYPE_CHECKING:  # archived writer module; not importable without the retired package
    from views import ViewData


CONTRACT_OBS = [
    "well", "dataset_id", "compound_name", "drug_key", "inchikey", "smiles",
    "dose_um", "dose_source", "time_h", "is_control", "vehicle", "plate",
    "control_group", "cell_type", "donor", "replicate", "atc", "moa",
    "mesh_classes", "pharm_atc_moa", "pharm_mesh", "annotation_source",
    "lib_size", "n_genes_detected",
]


@dataclass
class GateResult:
    name: str
    passed: bool
    detail: str = ""
    skipped: bool = False


def cross_view_gates(views: dict[str, ViewData],
                     universe: GeneUniverse) -> list[GateResult]:
    results = []

    # 1. identical var
    base = views[list(views)[0]].adata.var
    ok = all(view.adata.var.equals(base) for view in views.values())
    results.append(GateResult("var identical across all views", ok,
                              f"{len(base):,} genes"))

    # 2. temp_env tier nesting + value identity + shared controls
    narrow_to_wide = [("temp_env_atc-moa.h5ad", "temp_env_atc-moa-mesh.h5ad"),
                      ("temp_env_atc-moa-mesh.h5ad", "temp_env_all-smiles.h5ad")]
    for narrow, wide in narrow_to_wide:
        if narrow not in views or wide not in views:
            continue
        a, b = views[narrow].adata, views[wide].adata
        a_pos = {well: i for i, well in enumerate(a.obs["well"])}
        b_pos = {well: i for i, well in enumerate(b.obs["well"])}
        missing = sorted(set(a_pos) - set(b_pos))
        ok = not missing
        detail = f"{narrow} rows absent from {wide}: {missing[:5]}"
        if ok:
            for well, row_a in a_pos.items():
                row_b = b_pos[well]
                if not np.array_equal(a.X[row_a], b.X[row_b]):
                    ok, detail = False, f"X differs for {well}"
                    break
                if not np.allclose(a.layers["log1p_cp10k"][row_a],
                                   b.layers["log1p_cp10k"][row_b], atol=1e-6):
                    ok, detail = False, f"layer differs for {well}"
                    break
                for column in CONTRACT_OBS:
                    if a.obs.iloc[row_a][column] != b.obs.iloc[row_b][column]:
                        ok, detail = False, f"{column} differs for {well}"
                        break
                if not ok:
                    break
        results.append(GateResult(f"tier nesting {narrow} ⊆ {wide}", ok, detail))

    # 3. the 272 temp_env controls are identical across the three views
    temp_views = [views[name] for name in TEMP_ENV_VIEWS if name in views]
    if len(temp_views) >= 2:
        reference = temp_views[0].adata
        ref_ctrl = reference.obs["is_control"].to_numpy()
        ref_wells = reference.obs.loc[ref_ctrl, "well"].tolist()
        ok, detail = True, (f"{len(ref_wells)} control rows")
        for view in temp_views[1:]:
            adata = view.adata
            ctrl = adata.obs["is_control"].to_numpy()
            wells = adata.obs.loc[ctrl, "well"].tolist()
            if set(wells) != set(ref_wells):
                ok, detail = False, f"{view.name}: control wells differ"
                break
            position = {well: i for i, well in enumerate(adata.obs["well"])}
            for well in ref_wells:
                ref_row = reference.obs.index[reference.obs["well"] == well][0]
                ref_pos = reference.obs.index.get_loc(ref_row)
                row = position[well]
                if not np.array_equal(reference.X[ref_pos], adata.X[row]):
                    ok, detail = False, f"{view.name}: control counts differ at {well}"
                    break
                for column in ("control_group", "vehicle", "compound_name",
                               "dose_um"):
                    if (reference.obs.iloc[ref_pos][column]
                            != adata.obs.iloc[row][column]):
                        ok = False
                        detail = f"{view.name}: control {column} differs at {well}"
                        break
                if not ok:
                    break
            if not ok:
                break
        results.append(GateResult("temp_env controls identical across views", ok, detail))
    return results


