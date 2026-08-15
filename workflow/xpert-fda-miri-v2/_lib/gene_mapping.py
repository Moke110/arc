"""
Authoritative CM-dataset -> L1000 landmark gene mapping (workflow-owned).

Why this exists
---------------
Each CM finetune-data union dataset (dtoxs_primary, dtoxs_ipsc, gse244740,
gse247541, gse262419) identifies its genes in a different space, and several
use **retired / alternative symbols** for L1000 landmarks (e.g. WASHC5 ->
KIAA0196, ELP1 -> IKBKAP, JPT2 -> HN1L). Direct symbol matching would silently
zero those genes, so every union dataset resolves the 978 landmarks through its
own identifier space with legacy-symbol awareness.

Data sources (XPert-author + dataset data ONLY)
-----------------------------------------------
  1. L1000 landmark metadata  pert_model/xpert/metadata/l1000_gene_info_978.csv
     (author-shipped). Each landmark row carries gene_id (Entrez), gene_symbol,
     gene_name / gene_name_match (the 2015-era L1000 feature symbol, i.e. the
     legacy name a dataset may use), and ensembl_id.
  2. Each union dataset's own gene index, read directly from its raw files
     under pert_data/cm (symbols for the DToxS arms; Ensembl IDs for
     gse244740 / gse247541; SYMBOL_entrezid for gse262419).

`corrected_gene_name` is deliberately NOT used as a match key: for JPT2 the
author metadata records it as CRAMP1, a *different* gene (Entrez 57585), while
the DToxS-era symbol is HN1L (JPT2's real alias, Entrez 90861); and in
gse262419 the corrected names GET1 / TLCD3A are themselves *different* Entrez
genes (7746 / 2328) than the landmarks WRB (7485) / FAM57A (79850). Identity
always follows the author-recorded gene_name / gene_name_match / gene_id /
ensembl_id, never corrected_gene_name.

Output
------
One map CSV per union dataset under
03_cm_finetune_data/output/gene_maps/<dataset>_to_l1000_gene_map.csv
  columns: l1000_gene_symbol, l1000_gene_id, source_symbol, match_type
           (direct | legacy_alias | entrez | ensembl), source_entrez_verified
The DToxS-primary map is also consumed by step 02 (consensus baseline).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, XPERT_ROOT, CM_DATA_ROOT  # noqa: E402

L1000_META = XPERT_ROOT / "metadata" / "l1000_gene_info_978.csv"
GENE_MAP_DIR = (WORKFLOW_ROOT / "03_cm_finetune_data" / "output"
                / "gene_maps")

# Identifier spaces of the union datasets.
ID_SPACE = {
    "dtoxs_primary": "symbol",
    "dtoxs_ipsc": "symbol",
    "gse244740": "ensembl",
    "gse247541": "ensembl",
    "gse262419": "symbol_entrez",
}

CONTROL_STATES = {"CTRL", "CTRL-EtOH", "CTRL-DMSO"}

# gse262419 TempO-Seq panel provably lacks these two landmarks (their probes
# point at *different* Entrez genes, see module docstring); they are zeroed
# for that dataset's samples and reported, not silently mismatched.
GSE262419_PLATFORM_MISSING = {"WRB", "FAM57A"}


def load_l1000_meta():
    """Author L1000 landmark metadata, sorted by gene_idx."""
    meta = pd.read_csv(L1000_META)
    meta["gene_id"] = meta["gene_id"].astype(str)
    return meta.sort_values("gene_idx").reset_index(drop=True)


def author_name_index(meta=None):
    """Author metadata name -> (gene_id, gene_symbol) for all landmarks.

    Uses the author's gene_symbol / gene_name / gene_name_match fields only.
    `corrected_gene_name` is excluded (JPT2 <- CRAMP1 trap).
    """
    if meta is None:
        meta = load_l1000_meta()
    idx = {}
    for _, r in meta.iterrows():
        for c in ("gene_symbol", "gene_name", "gene_name_match"):
            v = r[c]
            if pd.notna(v) and str(v) not in ("nan", ""):
                idx[str(v)] = (str(r["gene_id"]), str(r["gene_symbol"]))
    return idx


def dataset_gene_index(dataset):
    """Source gene identifiers of a union dataset (set of strings).

    Identifier spaces per dataset:
      dtoxs_primary / dtoxs_ipsc  -> gene symbols (raw GEO read-count rows)
      gse244740 / gse247541       -> Ensembl IDs (version suffix stripped)
      gse262419                   -> SYMBOL_entrezid rows, split into the
                                     symbol part and the numeric part (the
                                     numeric part is NOT always the landmark
                                     Entrez — see module docstring).
    Returns (ids, id_space).
    """
    import glob

    if dataset in ("dtoxs_primary", "dtoxs_ipsc"):
        files = sorted(glob.glob(str(CM_DATA_ROOT / dataset
                                     / "*DGE-RNAseq-Read-Counts*.tsv.gz")))
        ids = set()
        for f in files:
            df = pd.read_csv(f, sep="\t", usecols=[0], compression="gzip")
            ids |= set(df.index.astype(str))
        return ids, "symbol"
    if dataset == "gse244740":
        df = pd.read_csv(CM_DATA_ROOT / dataset
                         / "GSE244740_processed_data_counts.txt.gz",
                         sep="\t", usecols=[0], compression="gzip")
        ids = {str(e) for e in df["Ensembl"]}
        return ids, "ensembl"
    if dataset == "gse247541":
        df = pd.read_csv(CM_DATA_ROOT / dataset / "GSE247541_raw_counts.csv.gz",
                         index_col=0, compression="gzip")
        ids = {str(e).split(".")[0] for e in df.index}
        return ids, "ensembl"
    if dataset == "gse262419":
        import glob
        files = sorted(glob.glob(str(CM_DATA_ROOT / dataset
                                     / "GSE262419_Plate*.csv.gz")))
        symbols, numbers = set(), set()
        for f in files:
            genes = pd.read_csv(f, usecols=[0], compression="gzip").iloc[:, 0]
            for g in genes.astype(str):
                parts = g.rsplit("_", 1)
                if len(parts) == 2:
                    symbols.add(parts[0])
                    numbers.add(parts[1])
        return {"symbols": symbols, "entrez": numbers}, "symbol_entrez"
    raise ValueError(f"Unknown union dataset: {dataset}")


def _resolve_symbol_space(meta, name_index, symbols):
    rows = []
    for _, r in meta.iterrows():
        sym, gid = r["gene_symbol"], str(r["gene_id"])
        if sym in symbols:
            source, mtype = sym, "direct"
        else:
            cands = []
            for c in ("gene_name", "gene_name_match"):
                v = r[c]
                if pd.notna(v) and str(v) in symbols \
                        and name_index.get(str(v)) == (gid, sym):
                    cands.append(str(v))
            if cands:
                source, mtype = cands[0], "legacy_alias"
            else:
                source, mtype = None, "unmatched"
        rows.append({"l1000_gene_symbol": sym, "l1000_gene_id": gid,
                     "source_symbol": source if source else "",
                     "match_type": mtype,
                     "source_entrez_verified": gid if source else ""})
    return rows


def resolve_landmarks(dataset, ids=None, meta=None, name_index=None):
    """Per-landmark source identifier for one union dataset.

    Returns a list of dicts with the map CSV columns. `match_type` is
    direct | legacy_alias | entrez | ensembl per the dataset's identifier
    space. gse262419 prefers the numeric (Entrez) part, then the symbol part.
    """
    if ids is None:
        ids, _ = dataset_gene_index(dataset)
    if meta is None:
        meta = load_l1000_meta()
    if name_index is None:
        name_index = author_name_index(meta)
    space = ID_SPACE[dataset]

    if space == "symbol":
        return _resolve_symbol_space(meta, name_index, ids)
    if space == "ensembl":
        rows = []
        for _, r in meta.iterrows():
            sym, gid = r["gene_symbol"], str(r["gene_id"])
            ens = str(r["ensembl_id"]).split(".")[0] if pd.notna(r["ensembl_id"]) else ""
            if ens and ens in ids:
                source, mtype = ens, "ensembl"
            else:
                source, mtype = None, "unmatched"
            rows.append({"l1000_gene_symbol": sym, "l1000_gene_id": gid,
                         "source_symbol": source if source else "",
                         "match_type": mtype,
                         "source_entrez_verified": gid if source else ""})
        return rows
    if space == "symbol_entrez":
        symbols, entrez = ids["symbols"], ids["entrez"]
        rows = []
        for _, r in meta.iterrows():
            sym, gid = r["gene_symbol"], str(r["gene_id"])
            if gid in entrez:
                source, mtype = gid, "entrez"
            elif sym in symbols:
                source, mtype = sym, "direct"
            else:
                cands = []
                for c in ("gene_name", "gene_name_match"):
                    v = r[c]
                    if pd.notna(v) and str(v) in symbols \
                            and name_index.get(str(v)) == (gid, sym):
                        cands.append(str(v))
                if cands:
                    source, mtype = cands[0], "legacy_alias"
                else:
                    source, mtype = None, "unmatched"
            rows.append({"l1000_gene_symbol": sym, "l1000_gene_id": gid,
                         "source_symbol": source if source else "",
                         "match_type": mtype,
                         "source_entrez_verified": gid if source else ""})
        return rows
    raise ValueError(f"Unknown identifier space: {space}")


def build_dataset_map(dataset, ids=None):
    """Build + return the map DataFrame for one union dataset."""
    rows = resolve_landmarks(dataset, ids=ids)
    return pd.DataFrame(rows)


def write_dataset_map(dataset, mapping=None, out_dir=None):
    """Write <dataset>_to_l1000_gene_map.csv under out_dir (or the default)."""
    if mapping is None:
        mapping = build_dataset_map(dataset)
    out_dir = Path(out_dir or GENE_MAP_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{dataset}_to_l1000_gene_map.csv"
    mapping.to_csv(out, index=False)
    return out


def load_dtoxs_primary_raw():
    """Union all raw GEO DGE read-count tables -> (symbol x well) raw counts.

    The raw GEO deposit (pert_data/cm/dtoxs_primary) uses compound well columns
    of the form <State>.<Cell>.<Experiment>.<Dish>.<Plate>.<Well>; design rows
    are joined by re-building that compound key from the design columns. The
    columns are left as the compound keys (well ids repeat across experiments,
    the compound key is unique).

    Returns (raw DataFrame indexed by gene symbol, ctrl_well_keys list).
    """
    import glob

    import numpy as np

    dtoxs_primary = CM_DATA_ROOT / "dtoxs_primary"
    ctrl_wells = []
    parts = []
    for counts_rel in sorted(glob.glob(
            str(dtoxs_primary / "GSE146096_DGE-RNAseq-Read-Counts*.tsv.gz"))):
        counts = pd.read_csv(counts_rel, sep="\t", index_col=0,
                             compression="gzip")
        design_rel = counts_rel.replace("GSE146096_DGE-RNAseq-Read-Counts",
                                        "GSE146096_DGE-RNAseq-Experiment-Design")
        design = pd.read_csv(design_rel, sep="\t", compression="gzip")
        key = (design["State"].astype(str) + "." + design["Cell"].astype(str)
               + "." + design["Experiment"].astype(str) + "."
               + design["Dish"].astype(str) + "." + design["Plate"].astype(str)
               + "." + design["Well"].astype(str))
        valid = key.isin(counts.columns)
        design = design[valid].copy()
        key = key[valid]
        ctrl = key[design["State"].isin(CONTROL_STATES)].tolist()
        ctrl_wells += ctrl
        counts = counts.loc[:, key.values]
        parts.append(counts)
    if not parts:
        raise FileNotFoundError("No GSE146096 read-count tables found under "
                                f"{dtoxs_primary}")
    raw = pd.concat(parts, axis=1)
    raw = raw.groupby(raw.index).sum()          # dedupe repeated gene rows
    raw.columns = raw.columns.astype(str)
    return raw, sorted(set(ctrl_wells))
