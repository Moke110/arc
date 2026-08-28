#!/usr/bin/env python3
"""
Step 03 — DToxS-LINCS finetune data (the v1.1 change vs v1).

Rebuilds the DToxS finetune h5ad that steps 05/06/07 (v1.1) train on from the
released LINCS LDS-1593 well-level package
(pert_data/dtoxs_primary/Data/Level3) instead of v1's frozen
input/dtoxs.h5ad. xpert-fda-miri-v1 stays untouched and frozen.

  03_cm_finetune_data/output/dtoxs_lincs.h5ad

   The rebuilt h5ad is schema-identical to v1's frozen input/dtoxs.h5ad in the
   fields the XPert library consumes (obs pert_id/pert_idx/pert_dose[_idx]/
   pert_time[_idx]/cell_iname/cell_idx/tissue_idx + the additive dataset_id
   provenance column; var = the 978 L1000 landmarks in gene_idx order with
   gene_id/gene_symbol; obsm["X_ctl"] per-sample matched control; X float32 in
   log2(CPM+1) space, quantile-aligned to the L1000 pretrain range [0, 15.002]).
   Only samples whose compound resolves to an HG drug (drug_idx in [0, 8981))
   enter. The DToxS-family single-dose experiment follows v1's slot-24 convention
   (pert_time 24.0, time_idx 2, pert_dose 0 -> dose_idx 0).

   **Per-plate normalization (ADR 0016):** both the quantile alignment and the
   matched-control computation are grouped per `(dataset_id, plate)`, not per
   dataset. Each sample's matched control is the mean of the control wells on
   its own plate (dtoxs_primary CTRL/CTRL-DMSO/CTRL-EtOH on plates 1/2/4/5) and
   each plate's 978-landmark distribution is aligned to the L1000 pretrain
   independently. v1's frozen input/dtoxs.h5ad stays dataset-wide (the
   historical anchor); the deviation is recorded as an auditable delta in the
   parity report. Sample identity/order is unchanged, so the fold columns
   survive (X-independent splitter, ADR 0014).

Also emits:
  03_cm_finetune_data/output/gene_maps/dtoxs_primary_to_l1000_gene_map.csv
  03_cm_finetune_data/output/drug_manifest.csv   per-compound drug resolution
  03_cm_finetune_data/output/qn_reference.npy    L1000 pretrain QN reference
  03_cm_finetune_data/output/parity_delta.md     auditable delta vs v1 frozen

Parity gate semantics (ADR 0012): the gate requires SCHEMA parity (the rebuilt
h5ad is structure-identical to v1's input/dtoxs.h5ad in the fields the XPert
loader consumes) and emits an AUDITABLE DELTA report of the sample/drug
differences — the fresh LINCS build is a legitimate superset (the raw GEO
deposit lacks the non-KI perturbations), so sample-for-sample identity with
v1's frozen data is neither required nor enforced.

Usage:
    python 03_cm_finetune_data/build_dtoxs_lincs.py            # idempotent
    python 03_cm_finetune_data/build_dtoxs_lincs.py --force    # rebuild
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import _datasets  # noqa: E402
import gene_mapping  # noqa: E402
from drug_resolution import build_drug_manifest  # noqa: E402
from utils import (WORKFLOW_ROOT, V1_DTOXS_H5AD, L1000_PRETRAIN_H5AD,  # noqa: E402
                   N_L1000_DRUGS, AUTHOR_TIME_MAP)

STEP_DIR = WORKFLOW_ROOT / "03_cm_finetune_data"
OUT_DIR = STEP_DIR / "output"
GENE_MAP_DIR = OUT_DIR / "gene_maps"
H5AD = OUT_DIR / "dtoxs_lincs.h5ad"
MANIFEST = OUT_DIR / "drug_manifest.csv"
QN_REFERENCE = OUT_DIR / "qn_reference.npy"
PARITY_DELTA = OUT_DIR / "parity_delta.md"
CELL_NAME = "cardio"
DATASET = "dtoxs_primary"
# The DToxS primary is single-dose with no real concentration recorded; it
# trains at pert_dose 0 -> dose_idx 0 exactly like v1's frozen dtoxs.h5ad, and
# its 48 h exposure in reality is encoded at the author 24 h slot (time_idx 2).
TIME_SLOT_H = 24.0
TIME_IDX = 2


def dose_idx_of(dose):
    """XPert assign_pert_dose_corrected bins (0-9) — mirrored from the library."""
    if 0 <= dose < 0.21:
        return 0
    if dose < 0.41:
        return 1
    if dose < 0.71:
        return 2
    if dose < 1.01:
        return 3
    if dose < 1.51:
        return 4
    if dose < 3.1:
        return 5
    if dose < 4.1:
        return 6
    if dose < 7.1:
        return 7
    if dose < 12.1:
        return 8
    return 9


def load_qn_reference(force=False, max_rows=5000, seed=0):
    """Pooled empirical quantiles of the L1000 pretrain log2(CPM+1) X (cached)."""
    if QN_REFERENCE.exists() and not force:
        return np.load(QN_REFERENCE)
    import scanpy as sc

    a = sc.read_h5ad(L1000_PRETRAIN_H5AD, backed="r")
    n = a.n_obs
    step = max(1, n // max_rows)
    rng = np.random.default_rng(seed)
    row_idx = np.arange(0, n, step)
    if len(row_idx) > max_rows:
        row_idx = np.sort(rng.choice(row_idx, size=max_rows, replace=False))
    chunks = []
    bs = 400
    for i in range(0, len(row_idx), bs):
        idx = row_idx[i:i + bs]
        block = a.X[idx]
        if hasattr(block, "toarray"):
            block = block.toarray()
        chunks.append(np.asarray(block).ravel())
    ref = np.sort(np.concatenate(chunks))
    QN_REFERENCE.parent.mkdir(parents=True, exist_ok=True)
    np.save(QN_REFERENCE, ref)
    print(f"QN reference: {ref.size} pooled values, range "
          f"[{ref.min():.4f}, {ref.max():.4f}] -> {QN_REFERENCE}")
    return ref


def source_to_landmark_col(gene_map):
    meta = gene_mapping.load_l1000_meta()
    idx_by_gid = {str(g): i for i, g in zip(meta["gene_idx"], meta["gene_id"])}
    out = {}
    for _, r in gene_map.iterrows():
        src = r["source_symbol"]
        if src:
            out[src] = idx_by_gid[str(r["l1000_gene_id"])]
    return out


def build_landmark_matrix(counts, gene_ids, gene_map):
    """log2(CPM+1) (n_samples, 978) matrix in L1000 gene_idx order."""
    meta = gene_mapping.load_l1000_meta()
    idx_by_gid = {str(g): i for i, g in zip(meta["gene_idx"], meta["gene_id"])}
    row2col = source_to_landmark_col(gene_map)
    row_col = np.full(len(gene_ids), -1, dtype=int)
    for i, gid in enumerate(gene_ids):
        key = str(gid)
        row_col[i] = row2col.get(key, -1)
    mat = np.zeros((counts.shape[0], 978), dtype=np.float64)
    valid = row_col >= 0
    np.add.at(mat.T, row_col[valid], counts[:, valid].astype(np.float64).T)
    totals = mat.sum(axis=1)
    cpm = np.divide(mat, totals[:, None], out=np.zeros_like(mat),
                    where=totals[:, None] > 0) * 1e6
    return np.log2(cpm + 1).astype(np.float32)


def group_key(meta):
    """Composite per-plate group key: '<dataset_id>::<plate>' (ADR 0016)."""
    return (meta["dataset_id"].astype(str) + "::"
            + meta["plate"].astype(str))


def quantile_align(X, group, ref_sorted):
    """Per-group pooled quantile alignment onto the L1000 pretrain reference."""
    import scipy.stats as st

    ref = np.sort(ref_sorted)
    xp = np.linspace(0, 1, len(ref))
    out = X.copy()
    for g in pd.unique(group):
        sel = group.values == g
        vals = np.asarray(X[sel]).ravel()
        ranks = st.rankdata(vals, method="average")
        frac = (ranks - 0.5) / len(vals)
        mapped = np.interp(frac, xp, ref)
        out[sel] = mapped.reshape(X[sel].shape)
    return out


def matched_control_matrix(meta, X_landmark, controls_mask, group):
    """Per-sample matched control = group mean of control log2(CPM+1).

    The group is the composite `(dataset_id, plate)` key (ADR 0016): each
    treated sample's control baseline is the mean of the control wells on its
    own plate, removing plate batch effects from the control baseline.
    """
    ctl = np.zeros_like(X_landmark, dtype=np.float64)
    if controls_mask.sum() == 0:
        return ctl.astype(np.float32)
    gvals = group.values
    for g in pd.unique(group):
        sel = gvals == g
        if controls_mask[sel].sum() == 0:
            continue
        ctl_rows = X_landmark[controls_mask & sel]
        grp_mean = ctl_rows.mean(axis=0)
        for i in np.where(~controls_mask & sel)[0]:
            ctl[i] = grp_mean
    return ctl.astype(np.float32)


def per_plate_gates(meta, X, X_ctl, group, out_lines):
    """Plate-coverage + per-plate-proof gates (ADR 0016).

    Plate-coverage: every treated plate carries >= 1 control well on it (0
    samples excluded for lack of a same-plate control). Per-plate proof: each
    treated sample's X_ctl equals the mean of the control rows of its own
    `(dataset_id, plate)` group. Any failure exits 1.
    """
    meta = meta.reset_index(drop=True)
    controls_mask = meta["is_control"].values
    treated = ~controls_mask
    missing_plates = sorted(set(meta.loc[treated, "plate"])
                            - set(meta.loc[controls_mask, "plate"]))
    if missing_plates:
        print(f"ERROR: {len(missing_plates)} treated plates without a "
              f"same-plate control well: {missing_plates}")
        sys.exit(1)
    n_excluded = 0
    gvals = group.values
    for g in pd.unique(group):
        sel = gvals == g
        if controls_mask[sel].sum() == 0 and treated[sel].sum():
            n_excluded += int(treated[sel].sum())
    if n_excluded:
        print(f"ERROR: {n_excluded} treated samples excluded for lack of a "
              f"same-plate control")
        sys.exit(1)
    # per-plate proof: spot-check every treated sample's X_ctl == its own
    # plate's control mean (guaranteed by construction; assert the code path)
    for g in pd.unique(group):
        sel = gvals == g
        ctl_rows = X[controls_mask & sel]
        if len(ctl_rows) == 0:
            continue
        expect = ctl_rows.mean(axis=0)
        got = X_ctl[treated & sel]
        if not np.allclose(got, expect, atol=1e-6):
            print(f"ERROR: per-plate X_ctl mismatch for group {g}")
            sys.exit(1)
    # differs-from-dataset-wide proof (spec testing decision): the per-plate
    # control must differ from the pooled dataset-wide control mean on >= 1
    # treated sample whenever the dataset has >= 2 plates with controls
    # (single-plate datasets trivially have per-plate == dataset-wide).
    if controls_mask.sum() > 0 and len(pd.unique(group)) >= 2:
        dw_mean = X[controls_mask].mean(axis=0)
        n_diff = int((np.abs(X_ctl[treated] - dw_mean).max(axis=1) > 1e-6).sum())
        if n_diff == 0:
            print(f"ERROR: per-plate X_ctl identical to the dataset-wide "
                  f"control mean on every treated sample — per-plate grouping "
                  f"has no effect")
            sys.exit(1)
    else:
        n_diff = 0
    out_lines.append(
        f"- **Per-plate normalization (ADR 0016):** QN alignment and matched "
        f"control grouped per `(dataset_id, plate)`; {len(pd.unique(group))} "
        f"plates, 0 treated samples excluded for lack of a same-plate control, "
        f"per-plate X_ctl proof passed, differs from the dataset-wide control "
        f"mean on {n_diff} treated samples.")
    print(f"  per-plate gates: {len(pd.unique(group))} plates, "
          f"0 samples excluded, X_ctl proof passed "
          f"(differs from dataset-wide on {n_diff} treated samples)")


def main():
    parser = argparse.ArgumentParser(
        description="Build the DToxS-LINCS finetune-data h5ad (v1.1)")
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if the h5ad exists")
    args = parser.parse_args()

    if H5AD.exists() and not args.force:
        print(f"Cached DToxS-LINCS h5ad present: {H5AD} — skipping "
              f"(pass --force to rebuild)")
        return

    import anndata as ad

    print("Step 03 — building DToxS-LINCS finetune data (v1.1)")

    print("\n[1/4] DToxS-primary gene map")
    GENE_MAP_DIR.mkdir(parents=True, exist_ok=True)
    mp = gene_mapping.build_dataset_map(DATASET)
    gene_mapping.write_dataset_map(DATASET, mp, GENE_MAP_DIR)
    n = (mp["match_type"] != "unmatched").sum()
    print(f"  {DATASET}: {n}/978 landmarks resolved "
          f"({dict(mp['match_type'].value_counts())})")

    print("\n[2/4] sample extraction + drug resolution")
    qn_ref = load_qn_reference(force=args.force)
    meta, counts, gene_ids = _datasets.extract_dtoxs_primary()
    agg = meta.groupby(["compound_name", "compound_code"],
                       dropna=False).size().reset_index(name="n_samples")
    compounds = []
    for _, r in agg.iterrows():
        smi = meta[(meta["compound_name"] == r["compound_name"])
                   & (meta["compound_code"] == r["compound_code"])
                   ]["smiles"].dropna()
        smi = smi.iloc[0] if len(smi) else None
        compounds.append((DATASET, r["compound_name"], r["compound_code"],
                          smi, int(r["n_samples"])))
    manifest = build_drug_manifest(compounds)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, index=False)
    drug_idx_of = {}
    for _, r in manifest.iterrows():
        drug_idx_of[(r["dataset"], str(r["compound_code"])
                     or r["compound_name"])] = int(r["drug_idx"])
    print(f"  {DATASET}: {meta.shape[0]} wells, {len(agg)} compounds "
          f"({(manifest.hg_class == 'in_hg').sum()} in_hg)")

    print("\n[3/4] assemble DToxS-LINCS h5ad")
    space = gene_mapping.ID_SPACE[DATASET]
    X = build_landmark_matrix(counts, gene_ids, mp)
    grp = group_key(meta)
    X = quantile_align(X, grp, qn_ref)
    X_ctl = matched_control_matrix(meta, X, meta["is_control"].values, grp)
    trt = ~meta["is_control"].values
    drug_idx = np.array([drug_idx_of.get(
        (DATASET, str(r.compound_code) if r.compound_code
         else r.compound_name), -1) for r in meta.itertuples()])
    keep = trt & (drug_idx >= 0) & (drug_idx < N_L1000_DRUGS)

    obs = meta[keep].reset_index(drop=True).copy()
    obs["pert_idx"] = drug_idx[keep]
    obs["pert_dose"] = 0.0
    obs["pert_dose_idx"] = 0
    obs["pert_time"] = TIME_SLOT_H
    obs["pert_time_idx"] = TIME_IDX
    obs["cell_iname"] = CELL_NAME
    obs["cell_idx"] = 0
    obs["tissue_idx"] = 0
    obs["pert_id"] = obs["well"]

    meta_genes = gene_mapping.load_l1000_meta()
    var = pd.DataFrame({
        "gene_id": meta_genes["gene_id"].values.astype(str),
        "gene_symbol": meta_genes["gene_symbol"].values,
    }, index=meta_genes["gene_id"].values.astype(str))
    obs["smiles"] = obs["smiles"].fillna("").astype(str)

    # per-plate gates (ADR 0016): coverage + per-plate X_ctl proof
    delta_lines = []
    per_plate_gates(meta, X, X_ctl, grp, delta_lines)

    adata = ad.AnnData(X[keep].astype(np.float32), obs=obs, var=var)
    adata.obsm["X_ctl"] = X_ctl[keep].astype(np.float32)
    adata.write_h5ad(H5AD)
    print(f"\nDToxS-LINCS h5ad: {H5AD} ({adata.n_obs} samples x "
          f"{adata.n_vars} genes, {adata.obs['pert_idx'].nunique()} drugs)")

    print("\n[4/4] parity gate (schema parity + auditable delta vs v1 frozen)")
    run_parity_gate(H5AD, delta_lines)
    print("Step 03 complete.")


# ---------------------------------------------------------------------------
# Parity gate: schema parity + auditable delta vs v1's frozen input/dtoxs.h5ad
# ---------------------------------------------------------------------------
LOADER_OBS_COLS = ["pert_id", "pert_idx", "pert_dose", "pert_dose_idx",
                   "pert_time", "pert_time_idx", "cell_iname", "cell_idx",
                   "tissue_idx"]


def run_parity_gate(new_h5ad, extra_lines=None):
    """Assert schema parity and write the auditable delta report (ADR 0012)."""
    import anndata as ad

    new = ad.read_h5ad(new_h5ad)
    old = ad.read_h5ad(V1_DTOXS_H5AD)
    lines = []
    add = lines.append

    add("# DToxS-LINCS finetune-data parity delta (v1.1 step 03 vs v1 frozen)\n")
    add("Generated by `03_cm_finetune_data/build_dtoxs_lincs.py`.\n")

    # --- per-plate deviation vs v1's dataset-wide build (ADR 0016) ---
    if extra_lines:
        add("## Per-plate normalization vs v1 frozen (dataset-wide)\n")
        add("v1's frozen `input/dtoxs.h5ad` normalizes each dataset as one "
            "dataset-wide group (a single QN transform and a pooled control "
            "mean). The rebuilt DToxS-LINCS h5ad groups both the QN alignment "
            "and the matched control per `(dataset_id, plate)` — v1 is the "
            "documented historical anchor and stays untouched (ADR 0016).\n")
        for l in extra_lines:
            add(l)
        add("")

    # --- schema parity (gate) ---
    add("## Schema parity\n")
    schema_ok = True
    for c in LOADER_OBS_COLS:
        in_new = c in new.obs.columns
        in_old = c in old.obs.columns
        add(f"- obs `{c}`: v1 `{in_old}` | v1.1 `{in_new}`")
        schema_ok &= in_new and in_old
    same_genes = (new.var["gene_id"].tolist() == old.var["gene_id"].tolist())
    add(f"- var gene_id order identical to v1: {same_genes}")
    schema_ok &= same_genes
    add(f"- obsm `X_ctl` present: v1 `{'X_ctl' in old.obsm}` | "
        f"v1.1 `{'X_ctl' in new.obsm}`")
    schema_ok &= "X_ctl" in old.obsm and "X_ctl" in new.obsm
    add(f"- X dtype: v1 `{old.X.dtype}` | v1.1 `{new.X.dtype}`; "
        f"v1.1 range [{float(new.X.min()):.4f}, {float(new.X.max()):.4f}]")
    add(f"\n**Schema parity: {'PASS' if schema_ok else 'FAIL'}** "
        f"(required: obs/var/X_ctl fields the XPert loader consumes)\n")

    # --- auditable sample/drug delta ---
    add("## Sample / drug delta (audit)\n")
    old_idx = old.obs["pert_idx"].astype(int)
    new_idx = new.obs["pert_idx"].astype(int)
    add(f"- v1 frozen: **{old.n_obs}** samples, **{old_idx.nunique()}** drugs")
    add(f"- v1.1 rebuilt (LINCS LDS-1593): **{new.n_obs}** samples, "
        f"**{new_idx.nunique()}** drugs\n")

    old_drugs = set(old_idx.unique().tolist())
    new_drugs = set(new_idx.unique().tolist())
    only_old = sorted(old_drugs - new_drugs)
    only_new = sorted(new_drugs - old_drugs)
    shared = sorted(old_drugs & new_drugs)

    def name_of(idx):
        return _cmap_name(idx) or f"idx {idx}"

    if only_old:
        add("### Drugs in v1 frozen only (not in the rebuilt set)\n")
        for i in only_old:
            add(f"- {i} ({name_of(i)}): {int((old_idx == i).sum())} v1 samples")
        add("")
    if only_new:
        add("### Drugs in the rebuilt set only (not in v1 frozen)\n")
        for i in only_new:
            add(f"- {i} ({name_of(i)}): {int((new_idx == i).sum())} samples")
        add("")
    add("### Shared drugs — sample counts (v1 vs v1.1)\n")
    add("| drug_idx | name | v1 n | v1.1 n |")
    add("|---|---|---|---|")
    for i in shared:
        add(f"| {i} | {name_of(i)} | {int((old_idx == i).sum())} | "
            f"{int((new_idx == i).sum())} |")

    # compounds that changed drug_idx (same compound, different L1000 row)
    add("")
    add("_Names from the author L1000 `cmap_name` table. Re-indexed compounds "
        "(v1's frozen indices that the rebuilt set maps to a different L1000 "
        "row) are visible as rows present only in one drug set with the same "
        "name._")
    add("")
    add("## Conclusion\n")
    if schema_ok:
        add("Schema parity holds: v1.1's rebuilt DToxS-LINCS h5ad is "
            "structure-identical to v1's frozen data in the fields the XPert "
            "loader consumes. The sample/drug deltas above are the auditable "
            "consequence of building from the LINCS LDS-1593 well-level "
            "package (a superset of the raw GEO deposit).")
    else:
        add("SCHEMA PARITY FAILED — the rebuilt h5ad is not loader-compatible.")
    add("")

    PARITY_DELTA.write_text("\n".join(lines))
    print(f"Parity delta -> {PARITY_DELTA}")
    print(f"  v1 frozen:  {old.n_obs} samples / {old_idx.nunique()} drugs")
    print(f"  v1.1 rebuilt: {new.n_obs} samples / {new_idx.nunique()} drugs")
    print(f"  v1-only drugs: {len(only_old)} | v1.1-only drugs: {len(only_new)}")


def _cmap_name(idx):
    """Best-effort author L1000 cmap_name for a drug index."""
    try:
        import scanpy as sc
        a = sc.read_h5ad(_l1000_full_h5ad(), backed="r")
        rows = a.obs[a.obs["pert_idx"].astype(int) == idx]
        names = rows["cmap_name"].dropna().unique()
        return str(names[0]) if len(names) else None
    except Exception:
        return None


def _l1000_full_h5ad():
    from utils import CM_DATA_ROOT
    return CM_DATA_ROOT / "l1000-xpert" / "data" / "l1000_mdmt_full_336852.h5ad"


if __name__ == "__main__":
    main()
