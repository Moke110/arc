"""View construction: corpus gene axis + one tier filter -> load-ready h5ad.

The five views share one writer, one schema gate and one control layer.  A view
is a filtered projection of the resolved source rows; the build never rewrites
a count, so a row retained by two gse262419 views is value-identical by
construction (and asserted in ``qa.py``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from .drugs import TIER_LABELS, CompoundTable
from .genes import GeneUniverse


class ViewError(RuntimeError):
    """A view violated the h5ad contract."""


@dataclass
class ViewData:
    name: str
    adata: ad.AnnData
    n_perturbed: int
    n_controls: int
    n_drugs: int


def map_counts(counts: np.ndarray, gene_map, universe: GeneUniverse) -> np.ndarray:
    """Sum source rows into the fixed universe order (duplicate probes summed)."""
    target = gene_map.targets(universe)
    keep = target >= 0
    if not keep.any():
        return np.zeros((counts.shape[0], len(universe)), dtype=np.int32)
    selected = counts[:, keep].astype(np.int64)
    target = target[keep]
    aggregation = sp.csr_matrix(
        (np.ones(len(target), dtype=np.int64), (np.arange(len(target)), target)),
        shape=(len(target), len(universe)))
    mapped = selected @ aggregation
    return np.asarray(mapped, dtype=np.int32)


def _replicate_index(meta: pd.DataFrame) -> np.ndarray:
    """Within-group replicate index, deterministic by well id."""
    group = np.where(
        meta["is_control"].to_numpy(),
        meta["control_group"].astype(str) + "|control|" + meta["vehicle"].astype(str),
        meta["control_group"].astype(str) + "|" + meta["drug_key"].astype(str)
        + "|" + meta["dose_um"].astype(str))
    order = np.lexsort((meta["well"].astype(str).to_numpy(), group))
    counts = np.zeros(len(meta), dtype=np.int32)
    last_group, running = None, 0
    for position in order:
        key = group[position]
        if key != last_group:
            running = 0
            last_group = key
        counts[position] = running
        running += 1
    return counts


def _select_rows(meta: pd.DataFrame, tier: str) -> np.ndarray:
    tier_column = TIER_LABELS[tier]
    controls = meta["is_control"].to_numpy()
    perturbed = ((~meta["is_control"]) & (~meta["is_blank"])
                 & (meta["entity_type"] == "small_molecule")
                 & meta[tier_column].astype(bool)).to_numpy()
    return np.flatnonzero(controls | perturbed)


def build_view(name: str, datasets: list[str], tier: str,
               sources: dict[str, object], universe: GeneUniverse,
               table: CompoundTable) -> ViewData:
    pieces, obs_rows = [], []
    for dataset_id in datasets:
        source = sources[dataset_id]
        meta = source.meta.copy()
        meta["control_group"] = (meta["dataset_id"].astype(str) + "::"
                                 + meta["plate"].astype(str))
        meta["replicate"] = _replicate_index(meta)
        selected = _select_rows(meta, tier)
        if len(selected) == 0:
            continue
        if source.mapped_counts is not None:
            counts = source.mapped_counts[selected]
        else:
            counts = map_counts(source.counts[selected], source.gene_map, universe)
        pieces.append(counts)
        rows = meta.iloc[selected].copy()
        rows["source_row"] = selected
        obs_rows.append(rows)

    if not pieces:
        raise ViewError(f"{name}: no rows selected")

    counts = np.vstack(pieces).astype(np.int32, copy=False)
    meta = pd.concat(obs_rows, ignore_index=True)

    lib_size = counts.sum(axis=1, dtype=np.int64)
    if (lib_size <= 0).any():
        bad = meta.loc[lib_size <= 0, "well"].tolist()[:5]
        raise ViewError(f"{name}: zero library size for wells {bad}")
    n_detected = (counts > 0).sum(axis=1).astype(np.int32)

    obs = pd.DataFrame({
        "well": meta["well"].astype(str).values,
        "dataset_id": meta["dataset_id"].astype(str).values,
        "compound_name": meta["compound_name"].astype(str).values,
        "drug_key": meta["drug_key"].fillna("").astype(str).values,
        "inchikey": meta["inchikey"].fillna("").astype(str).values,
        "smiles": meta["smiles"].fillna("").astype(str).values,
        "dose_um": meta["dose_um"].astype(float).values,
        "dose_source": meta["dose_source"].fillna("").astype(str).values,
        "time_h": meta["time_h"].astype(float).values,
        "is_control": meta["is_control"].astype(bool).values,
        "vehicle": meta["vehicle"].astype(str).values,
        "plate": meta["plate"].astype(str).values,
        "control_group": meta["control_group"].astype(str).values,
        "cell_type": meta["cell_type"].astype(str).values,
        "donor": meta["donor"].astype(str).values,
        "replicate": meta["replicate"].astype(np.int32).values,
        "atc": meta["atc"].fillna("").astype(str).values,
        "moa": meta["moa"].fillna("").astype(str).values,
        "mesh_classes": meta["mesh_classes"].fillna("").astype(str).values,
        "pharm_atc_moa": meta["pharm_atc_moa"].astype(bool).values,
        "pharm_mesh": meta["pharm_mesh"].astype(bool).values,
        "annotation_source": meta["annotation_source"].fillna("").astype(str).values,
        "lib_size": lib_size.astype(np.int64),
        "n_genes_detected": n_detected,
    }, index=pd.Index(meta["well"].astype(str).values, name="well"))

    var = universe.var.copy()
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["log1p_cp10k"] = np.log1p(
        counts / lib_size[:, None] * 1e4).astype(np.float32)

    control_mask = obs["is_control"].to_numpy()
    groups = obs.loc[control_mask, "control_group"]
    if groups.empty:
        raise ViewError(f"{name}: no control rows")
    means = []
    for group in sorted(groups.unique()):
        rows = control_mask & (obs["control_group"].to_numpy() == group)
        means.append(counts[rows].mean(axis=0))
    control_mean = pd.DataFrame(np.vstack(means), index=sorted(groups.unique()),
                                columns=var.index)
    control_mean.index.name = "control_group"
    adata.uns["control_mean"] = control_mean

    view = ViewData(
        name=name,
        adata=adata,
        n_perturbed=int((~control_mask).sum()),
        n_controls=int(control_mask.sum()),
        n_drugs=int(obs.loc[~control_mask, "drug_key"].nunique()),
    )
    validate_view(view, universe, tier, table)
    return view


def validate_view(view: ViewData, universe: GeneUniverse, tier: str,
                  table: CompoundTable) -> None:
    adata = view.adata
    obs, var, X = adata.obs, adata.var, adata.X
    name = view.name

    # 1. gene axis
    if list(var.index) != universe.symbols:
        raise ViewError(f"{name}: var is not the frozen corpus gene axis")
    if var.index.duplicated().any() or (var.index.astype(str) == "").any():
        raise ViewError(f"{name}: duplicate or empty var symbols")
    for column in ("entrez_id", "ensembl_id", "match_type"):
        if column not in var:
            raise ViewError(f"{name}: var missing {column}")

    # 2. row identity and dose
    if obs["well"].duplicated().any():
        dup = obs.loc[obs["well"].duplicated(), "well"].tolist()[:5]
        raise ViewError(f"{name}: duplicate well ids {dup}")
    control = obs["is_control"].to_numpy()
    if (obs.loc[~control, "dose_um"] <= 0).any():
        bad = obs.loc[~control & (obs["dose_um"] <= 0), "well"].tolist()[:5]
        raise ViewError(f"{name}: perturbed rows without a positive dose {bad}")
    if (obs.loc[control, "dose_um"] != 0).any():
        raise ViewError(f"{name}: control rows with a non-zero dose")
    for column in ("drug_key", "smiles", "inchikey"):
        if (obs.loc[~control, column].astype(str) == "").any():
            bad = obs.loc[~control & (obs[column].astype(str) == ""), "well"].tolist()[:5]
            raise ViewError(f"{name}: perturbed rows missing {column} {bad}")
    if (obs.loc[control, "drug_key"].astype(str) != "").any():
        raise ViewError(f"{name}: control rows carry a drug_key")

    # 3. every control group has a control
    groups = obs.loc[control, "control_group"].astype(str)
    if groups.empty:
        raise ViewError(f"{name}: no controls")
    if not (obs.loc[~control, "control_group"].astype(str)
            .isin(set(groups)).all()):
        missing = sorted(set(obs.loc[~control, "control_group"].astype(str))
                         - set(groups))
        raise ViewError(f"{name}: control groups without a control {missing}")

    # 4. counts
    if not np.issubdtype(X.dtype, np.integer):
        raise ViewError(f"{name}: X is {X.dtype}, expected integer raw counts")
    if (X < 0).any():
        raise ViewError(f"{name}: negative counts")
    if not np.isfinite(X).all():
        raise ViewError(f"{name}: non-finite counts")
    if (obs["lib_size"] <= 0).any() or (obs["n_genes_detected"] <= 0).any():
        raise ViewError(f"{name}: empty library")
    if (obs["lib_size"].to_numpy() != X.sum(axis=1)).any():
        raise ViewError(f"{name}: lib_size does not equal the row sum")
    if (obs["n_genes_detected"].to_numpy() != (X > 0).sum(axis=1)).any():
        raise ViewError(f"{name}: n_genes_detected does not match X")

    # 5. convenience layer
    expected = np.log1p(X / obs["lib_size"].to_numpy()[:, None] * 1e4)
    if not np.allclose(adata.layers["log1p_cp10k"], expected, atol=1e-5):
        raise ViewError(f"{name}: log1p_cp10k does not match the raw counts")

    # 6. no integer model encodings
    forbidden = {"pert_idx", "cell_idx", "dose_idx", "pert_dose_idx"}
    if forbidden & set(obs.columns):
        raise ViewError(f"{name}: model encodings leaked into obs "
                        f"{sorted(forbidden & set(obs.columns))}")

    # 7. tier closure
    allowed = table.members(tier)
    used = set(obs.loc[~control, "drug_key"])
    outside = sorted(used - allowed)
    if outside:
        raise ViewError(f"{name}: drug keys outside {tier}: {outside[:5]}")

    # 8. control means
    means = adata.uns.get("control_mean")
    if means is None or len(means) != groups.nunique():
        raise ViewError(f"{name}: uns['control_mean'] does not cover every group")
    if list(means.columns) != universe.symbols:
        raise ViewError(f"{name}: uns['control_mean'] has the wrong gene axis")
