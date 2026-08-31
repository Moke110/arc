"""
Authoritative DToxS -> L1000 landmark gene mapping (workflow-owned).

Why this exists
---------------
DToxS (3'-DGE, 2015-era symbols) contains all 978 L1000 landmark genes, but for
18 landmarks the DToxS row uses a **retired / alternative symbol** that does
not equal the L1000 `gene_symbol` (e.g. WASHC5 -> DToxS KIAA0196, ELP1 -> IKBKAP,
JPT2 -> HN1L). Direct symbol matching therefore silently zeroes those genes.

Data sources (XPert-author + DToxS ONLY)
----------------------------------------
  1. L1000 landmark metadata  pert_model/xpert/metadata/l1000_gene_info_978.csv
     (author-shipped). Each landmark row carries gene_id (Entrez), the current
     gene_symbol/gene_name_match, AND the 2015-era L1000 feature symbol in
     `gene_name`. The author's own `gene_name` column lists every legacy symbol
     the DToxS table uses (e.g. JPT2's gene_name = HN1L), so alias resolution
     needs no external database.
  2. DToxS raw Level3 read-count tables (symbol-indexed) and the frozen
     finetune h5ad (03_cm_finetune_data/output/dtoxs.h5ad, var gene_symbol -> gene_id).

`corrected_gene_name` is deliberately NOT used as a match key: for JPT2 the
author metadata records it as CRAMP1, a *different* gene (Entrez 57585), while
the DToxS table carries HN1L (JPT2's real alias, Entrez 90861). The mapping
therefore prefers gene_name / gene_name_match, which is what the DToxS-era
symbols match.

Output
------
  02_cm_baseline/output/dtoxs_to_l1000_gene_map.csv
    columns: l1000_gene_symbol, l1000_gene_id, dtoxs_source_symbol,
             match_type (direct | legacy_alias), source_entrez_verified,
             dtoxs_raw_total_all_wells, dtoxs_raw_total_ctrl_wells
    `source_entrez_verified` is the Entrez the source symbol is verified to:
    direct matches resolve through the DToxS h5ad (gene_symbol -> gene_id);
    legacy aliases resolve through the author's own gene_name pairing.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, XPERT_ROOT, DTOXS_ROOT  # noqa: E402

L1000_META = XPERT_ROOT / "metadata" / "l1000_gene_info_978.csv"
OUTPUT = (WORKFLOW_ROOT / "02_cm_baseline" / "output"
          / "dtoxs_to_l1000_gene_map.csv")
H5AD = WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "dtoxs.h5ad"

CONTROL_STATES = {"CTRL", "CTRL-EtOH", "CTRL-DMSO"}


def load_dtoxs_raw_counts():
    """Union all Level3 read-count tables -> (symbol x well) raw counts.

    Returns (raw DataFrame indexed by gene symbol, ctrl_wells list).
    """
    import glob

    import numpy as np

    level3 = DTOXS_ROOT / "Data" / "Level3"
    ctrl_wells = []
    parts = []
    for counts_rel in sorted(
            glob.glob(str(level3 / "**" / "PromoCell-DGE-Read-Counts*.tsv"),
                      recursive=True)):
        counts = pd.read_csv(counts_rel, sep="\t", index_col=0)
        design_rel = counts_rel.replace("PromoCell-DGE-Read-Counts",
                                        "PromoCell-DGE-Experiment-Design")
        design = pd.read_csv(design_rel, sep="\t")
        ctrl = design[design["State"].isin(CONTROL_STATES)]["Well"].tolist()
        ctrl = [w for w in ctrl if w in counts.columns]
        ctrl_wells += ctrl
        parts.append(counts)
    raw = pd.concat(parts, axis=1)
    raw = raw.groupby(raw.index).sum()          # dedupe repeated gene rows
    raw.columns = raw.columns.astype(str)
    return raw, sorted(set(ctrl_wells))


def h5ad_symbol_to_entrez():
    """DToxS h5ad var gene_symbol -> gene_id (frozen, DToxS-provided)."""
    import anndata as ad

    a = ad.read_h5ad(H5AD)
    gs = a.var["gene_symbol"].astype(str)
    gid = a.var["gene_id"].astype(str)
    return dict(zip(gs, gid))


def author_name_index():
    """Author metadata name -> (gene_id, gene_symbol) for all landmarks.

    Uses the author's gene_symbol / gene_name / gene_name_match fields only.
    `corrected_gene_name` is excluded: for JPT2 it is CRAMP1 (a different gene,
    Entrez 57585) while the DToxS-era name is HN1L (Entrez 90861).
    """
    meta = pd.read_csv(L1000_META)
    meta["gene_id"] = meta["gene_id"].astype(str)
    idx = {}
    for _, r in meta.iterrows():
        for c in ("gene_symbol", "gene_name", "gene_name_match"):
            v = r[c]
            if pd.notna(v) and str(v) not in ("nan", ""):
                idx[str(v)] = (str(r["gene_id"]), str(r["gene_symbol"]))
    return idx


def build_mapping(raw, ctrl_wells):
    """Per-landmark authoritative DToxS source symbol (author+DToxS data only).

    For each landmark: the source is the landmark's own gene_symbol when DToxS
    carries it (direct), else the author-recorded legacy name (gene_name /
    gene_name_match) present in DToxS (legacy_alias). A DToxS symbol claimed by
    more than one landmark is treated as unmatched and reported — none exist.
    """
    meta = pd.read_csv(L1000_META)
    meta["gene_id"] = meta["gene_id"].astype(str)
    dtoxs_symbols = set(raw.index)
    name_to_lm = author_name_index()
    h5ad_sym2entrez = h5ad_symbol_to_entrez()
    # index raw totals once
    raw_all = raw.sum(axis=1)
    raw_ctrl = raw[ctrl_wells].sum(axis=1) if ctrl_wells else pd.Series(0, index=raw.index)

    rows = []
    for _, r in meta.sort_values("gene_idx").iterrows():
        sym, gid = r["gene_symbol"], str(r["gene_id"])
        if sym in dtoxs_symbols:
            source, mtype = sym, "direct"
        else:
            # legacy: author-recorded names of this landmark present in DToxS,
            # excluding names claimed by any other landmark.
            cands = []
            for c in ("gene_name", "gene_name_match"):
                v = r[c]
                if pd.notna(v) and str(v) in dtoxs_symbols \
                        and name_to_lm.get(str(v)) == (gid, sym):
                    cands.append(str(v))
            if cands:
                source, mtype = cands[0], "legacy_alias"
            else:
                source, mtype = None, "unmatched"

        if source:
            # direct -> resolve through the DToxS h5ad; legacy -> author pairing
            verified = h5ad_sym2entrez.get(source, gid) if mtype == "direct" else gid
        else:
            verified = ""

        rows.append({
            "l1000_gene_symbol": sym,
            "l1000_gene_id": gid,
            "dtoxs_source_symbol": source if source else "",
            "match_type": mtype,
            "source_entrez_verified": verified,
            "dtoxs_raw_total_all_wells": int(raw_all.get(source, 0)) if source else 0,
            "dtoxs_raw_total_ctrl_wells": int(raw_ctrl.get(source, 0)) if source else 0,
        })
    return pd.DataFrame(rows)


def main():
    raw, ctrl_wells = load_dtoxs_raw_counts()
    mapping = build_mapping(raw, ctrl_wells)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(OUTPUT, index=False)

    n_direct = (mapping["match_type"] == "direct").sum()
    n_alias = (mapping["match_type"] == "legacy_alias").sum()
    n_unmatched = (mapping["match_type"] == "unmatched").sum()
    print(f"DToxS symbols: {len(raw.index)} | control wells: {len(ctrl_wells)}")
    print(f"L1000 landmarks mapped: direct={n_direct}, legacy_alias={n_alias}, "
          f"unmatched={n_unmatched} / {len(mapping)}")
    if n_alias:
        print("\nLegacy-alias matches (DToxS old symbol -> L1000 landmark):")
        for _, r in mapping[mapping["match_type"] == "legacy_alias"].iterrows():
            print(f"  {r['dtoxs_source_symbol']:10s} -> {r['l1000_gene_symbol']:10s} "
                  f"(entrez {r['l1000_gene_id']}, verified {r['source_entrez_verified']})")
    print(f"\nSaved: {OUTPUT}")
    return mapping


if __name__ == "__main__":
    main()
