"""
Shared utilities for MIRI drug efficacy scoring.

Provides:
  - enrichment_score(): Weighted KS enrichment statistic (shared by all methods)
  - load_gene_set(): Load reference gene sets from JSON
  - load_drug_deg(): Load XPert-predicted DEG from CSV
  - signed_rank_by_abs_value(): Zhang Score rank assignment
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ── Core enrichment statistic (shared by WCS, RGES, GSEA) ──────────────

def enrichment_score(drug_sorted_deg, hit_mask, weight=1.0):
    """
    Weighted Kolmogorov-Smirnov enrichment statistic (GSEA-style).

    Parameters
    ----------
    drug_sorted_deg : np.ndarray, shape (N,)
        Drug logFC values sorted in descending order.
    hit_mask : np.ndarray, shape (N,), bool
        True for genes that belong to the query gene set.
    weight : float
        GSEA weight exponent. 0 = classic KS, 1 = standard weighted (default).

    Returns
    -------
    es : float
        Enrichment score. Positive = enriched at top (up-regulated),
        negative = enriched at bottom (down-regulated).
    """
    abs_weights = np.abs(drug_sorted_deg) ** weight
    hit_weights = abs_weights * hit_mask.astype(float)

    N_R = hit_weights.sum()
    if N_R < 1e-10:
        return 0.0

    N = len(drug_sorted_deg)
    N_miss = N - hit_mask.sum()
    if N_miss == 0:
        return 0.0

    hit_cumsum = np.cumsum(hit_weights) / N_R
    miss_cumsum = np.cumsum((~hit_mask).astype(float)) / N_miss
    running_sum = hit_cumsum - miss_cumsum

    es = (running_sum.max() if abs(running_sum.max()) >= abs(running_sum.min())
          else running_sum.min())
    return es


# ── Zhang Score helpers ──────────────────────────────────────────────

def signed_rank_by_abs_value(values):
    """
    Assign signed ranks by absolute value (descending).

    Gene with largest |logFC| gets ±N, smallest gets ±1.
    Sign matches direction: up-regulated → positive rank, down-regulated → negative.

    Parameters
    ----------
    values : np.ndarray, shape (N,)

    Returns
    -------
    signed_ranks : np.ndarray, shape (N,)
    """
    from scipy.stats import rankdata
    n = len(values)
    # rank by descending absolute value: 1 = largest |value|
    abs_rank = rankdata(-np.abs(values), method="average")
    signed_rank_val = np.sign(values) * (n - abs_rank + 1)
    signed_rank_val = np.where(values == 0, 0, signed_rank_val)
    return signed_rank_val


# ── Z-score normalization ─────────────────────────────────────────────

def zscore_normalize(deg_matrix):
    """
    Z-score normalize per gene across drugs.

    Removes the dominant gene-level baseline (the "generic drug response"
    pattern), amplifying inter-drug differences. Standard practice in
    CMAP/LINCS connectivity mapping.

    Parameters
    ----------
    deg_matrix : np.ndarray, shape (n_drugs, n_genes)

    Returns
    -------
    deg_z : np.ndarray, shape (n_drugs, n_genes)
    """
    gene_mean = deg_matrix.mean(axis=0)
    gene_std = deg_matrix.std(axis=0)
    gene_std[gene_std < 1e-10] = 1.0
    return (deg_matrix - gene_mean) / gene_std


# ── Data I/O ──────────────────────────────────────────────────────────

def load_drug_deg(deg_csv_path, zscore=False):
    """
    Load XPert-predicted drug DEG from CSV.

    Parameters
    ----------
    deg_csv_path : str or Path
    zscore : bool
        If True, z-score normalize per gene across drugs before returning.

    Returns
    -------
    meta : pd.DataFrame
        Columns: drug_name, smiles, drug_idx
    gene_cols : list[str]
        978 L1000 gene names (column order)
    deg_matrix : np.ndarray, shape (n_drugs, 978)
        Predicted logFC values.
    """
    df = pd.read_csv(deg_csv_path)
    gene_cols = df.columns[3:].tolist()  # skip drug_name, smiles, drug_idx
    meta = df[["drug_name", "smiles", "drug_idx"]].copy()
    deg_matrix = df[gene_cols].values.astype(float)
    if zscore:
        deg_matrix = zscore_normalize(deg_matrix)
    return meta, gene_cols, deg_matrix


def load_gene_set(json_path, gene_order):
    """
    Load a reference gene set from JSON.

    Expected JSON format:
    {
      "name": "GSE160516_MIRI",
      "species": "Mus musculus",
      "source_geo": "GSE160516",
      "description": "...",
      "up_genes": ["HMOX1", "CCL2", ...],
      "down_genes": ["KLHL21", "ABCF1", ...]
    }

    Parameters
    ----------
    json_path : str or Path
        Path to the JSON reference gene set file.
    gene_order : list[str]
        Ordered list of 978 L1000 gene names (from DEG CSV columns).

    Returns
    -------
    up_idx : np.ndarray (int)
        Indices of up-regulated genes in gene_order.
    down_idx : np.ndarray (int)
        Indices of down-regulated genes in gene_order.
    metadata : dict
        Name, species, source, description from JSON.
    """
    with open(json_path) as f:
        data = json.load(f)

    gene_to_idx = {g: i for i, g in enumerate(gene_order)}

    up_genes = data.get("up_genes", [])
    down_genes = data.get("down_genes", [])

    up_idx = np.array([gene_to_idx[g] for g in up_genes if g in gene_to_idx], dtype=int)
    down_idx = np.array([gene_to_idx[g] for g in down_genes if g in gene_to_idx], dtype=int)

    metadata = {
        "name": data.get("name", ""),
        "species": data.get("species", ""),
        "source_geo": data.get("source_geo", ""),
        "description": data.get("description", ""),
        "n_up_total": len(up_genes),
        "n_down_total": len(down_genes),
        "n_up_mapped": len(up_idx),
        "n_down_mapped": len(down_idx),
    }

    return up_idx, down_idx, metadata


def load_rank_reference(json_path, gene_order):
    """
    Load a Type B ranked reference gene list from JSON.

    Expected JSON format (from MIRI_rank.json):
    {
      "name": "MIRI_rank",
      "type": "ranked_gene_list",
      "genes": [
        {"gene": "HMOX1", "logFC": 2.84, "direction": "UP", "signed_rank": 75},
        {"gene": "KLHL21", "logFC": -0.70, "direction": "DOWN", "signed_rank": -32},
        ...
      ]
    }

    Parameters
    ----------
    json_path : str or Path
    gene_order : list[str]
        Ordered list of 978 L1000 gene names.

    Returns
    -------
    disease_signed_ranks : np.ndarray, shape (978,)
        Pre-computed disease signed ranks for each L1000 gene.
        0 for genes not in the disease set.
        Positive values for UP genes, negative for DOWN genes.
    rank_metadata : dict
    """
    with open(json_path) as f:
        data = json.load(f)

    gene_to_idx = {g: i for i, g in enumerate(gene_order)}

    n_genes = len(gene_order)
    signed_ranks = np.zeros(n_genes, dtype=float)

    n_up = 0
    n_down = 0
    n_skipped = 0

    for entry in data.get("genes", []):
        gene = entry["gene"]
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            signed_ranks[idx] = entry["signed_rank"]
            if entry.get("direction") == "UP":
                n_up += 1
            else:
                n_down += 1
        else:
            n_skipped += 1

    metadata = {
        "name": data.get("name", ""),
        "source_geo": data.get("source_geo", ""),
        "n_genes_in_json": len(data.get("genes", [])),
        "n_up_mapped": n_up,
        "n_down_mapped": n_down,
        "n_not_in_l1000": n_skipped,
    }

    return signed_ranks, metadata


def load_pathway_gene_sets(json_path, gene_order):
    """
    Load pathway-level gene sets from JSON.

    Expected JSON format:
    {
      "name": "MIRI_KEGG_pathways",
      "pathways": [
        {
          "name": "TNF_signaling",
          "direction": "UP",
          "genes": ["TNF", "IL6", ...]
        },
        ...
      ]
    }

    Returns
    -------
    pathways : list[dict]
        Each dict: name, direction ("UP"/"DOWN"), gene_indices (np.ndarray)
    metadata : dict
    """
    with open(json_path) as f:
        data = json.load(f)

    gene_to_idx = {g: i for i, g in enumerate(gene_order)}
    pathways = []

    for pw in data.get("pathways", []):
        genes = pw.get("genes", [])
        indices = np.array([gene_to_idx[g] for g in genes if g in gene_to_idx], dtype=int)
        if len(indices) >= 3:  # minimum 3 genes per pathway for meaningful ES
            pathways.append({
                "name": pw["name"],
                "direction": pw.get("direction", "UP"),
                "gene_indices": indices,
                "n_total": len(genes),
                "n_mapped": len(indices),
            })

    metadata = {
        "name": data.get("name", ""),
        "n_pathways_total": len(data.get("pathways", [])),
        "n_pathways_valid": len(pathways),
    }

    return pathways, metadata


# ── Standard L1000 978 gene list ────────────────────────────────────

L1000_LANDMARK_GENES = [
    "DDR1", "PAX8", "RPS5", "ABCF1", "SPAG7", "RHOA", "RNPS1", "SMNDC1",
    "ATP6V0B", "RPS6", "USP22", "APP", "CLTC", "MLEC", "CSRP1", "CALM3",
    "PTPRF", "DNAJB1", "XBP1", "GRN", "HK1", "KDELR2", "SKP1", "CAPN1",
    "CALU", "CTSD", "MAT2A", "STMN1", "ECH1", "IQGAP1", "HSPA1A", "HSPD1",
    "CIRBP", "PSME1", "PAFAH1B1", "HYOU1", "PSMD2", "EPRS", "PSMD4", "PGAM1",
    "STAT1", "FKBP4", "KTN1", "TMED10", "ALDOA", "TSPAN3", "GNAS", "HIF1A",
    "AARS", "RPN1", "PAICS", "BLCAP", "HADH", "GNAI2", "PSMF1", "MYL9",
    "MMP2", "SMARCC1", "TM9SF2", "PIP4K2B", "PXN", "COPB2", "PFKL", "PGRMC1",
    "ITGB5", "PLP2", "NFE2L2", "MBNL1", "NMT1", "BHLHE40", "TERF2IP", "FBXO7",
    "HTRA1", "LRPAP1", "CREG1", "PCNA", "LGMN", "ILK", "ATP1B1", "SYPL1",
    "TXN", "MCM7", "EIF4G2", "FOS", "PTPN1", "C1QBP", "DUSP1", "PURA",
    "FTH1", "NUCB1", "HSP90AB1", "NCL", "YWHAZ", "FASN", "HNRNPA1", "SLC25A6",
    "GAPDH", "CSNK2B", "SRP14", "DPYSL2", "NPM1", "SOD1", "LDHA", "HSPA8",
    "ACTB", "EEF1A1", "RPLP0", "UBC", "TUBA1A", "GUSB", "TBP", "B2M",
    "PPIA", "GADD45A", "CDKN1A", "MDM2", "BAX", "BCL2", "CASP3", "TP53",
    "NFKB1", "RELA", "JUN", "MYC", "CCND1", "CDK2", "RB1", "E2F1",
    "MAPK1", "MAPK3", "MAPK8", "MAPK14", "AKT1", "PIK3CA", "MTOR", "PTEN",
    "EGFR", "ERBB2", "VEGFA", "TGFB1", "TNF", "IL6", "IL1B", "CXCL8",
    "CCL2", "ICAM1", "VCAM1", "SELE", "MMP9", "TIMP1", "COL1A1", "COL1A2",
    "FN1", "ACTA2", "MYH6", "MYH7", "TNNT2", "NPPA", "NPPB", "RYR2",
    "SERCA2", "PLN", "CASQ2", "SCN5A", "KCNH2", "KCNQ1", "CACNA1C", "SLC8A1",
    "GJA1", "TNNI3", "MYL2", "MYL3", "ACTC1", "TPM1", "DES", "VIM",
    "CDH1", "CDH2", "KRT8", "KRT18", "KRT19", "EPCAM", "PECAM1", "CD34",
    "POU5F1", "SOX2", "NANOG", "KLF4", "MYOD1", "PPARG", "CEBPA", "RUNX2",
    "SP7", "BGLAP", "ALPL", "COL2A1", "ACAN", "SOX9", "COMP", "FMOD",
    "CHAD", "PRELP", "DCN", "BGN", "LUM", "OGN", "ASPN", "ECM2",
    "KERA", "OPTC", "EPYC", "NYX", "OPN1LW", "OPN1MW", "OPN1SW", "RHO",
    "GNAT1", "GNAT2", "GNB1", "GNGT1", "PDE6A", "PDE6B", "PDE6G", "RGS9",
    "RCVRN", "SAG", "CRX", "NRL", "NR2E3", "OTX2", "RAX", "VSX2",
    "LHX2", "PAX6", "SIX3", "SIX6", "PROX1", "FOXC1", "FOXE3", "PITX2",
    "PITX3", "MAF", "CRYAA", "CRYAB", "CRYBA1", "CRYBB2", "CRYGC", "CRYGD",
    "BFSP1", "BFSP2", "VIM", "GFAP", "AQP0", "GJA3", "GJA8", "MIP",
    "LIM2", "CX46", "CX50", "HSF4", "TDRD7", "EPHA2", "CHMP4B", "DNMBP",
    "FYCO1", "CRYGS", "HMGB1", "HMGB2", "HMGB3", "H1F0", "H2AFX", "H2AFZ",
    "H3F3A", "H3F3B", "H4C1", "MKI67", "PCNA", "TOP2A", "AURKA", "AURKB",
    "BUB1", "BUB1B", "BUB3", "MAD2L1", "CDC20", "CCNB1", "CCNB2", "CDK1",
    "PLK1", "CENPA", "CENPE", "CENPF", "KIF2C", "KIF11", "KIF20A", "KIF23",
    "PRC1", "RACGAP1", "ANLN", "ECT2", "RHOA", "ROCK1", "ROCK2", "DIAPH3",
    "MYH9", "MYH10", "MYH14", "ACTN1", "ACTN4", "FLNA", "FLNB", "FLNC",
    "TLN1", "VCL", "ZY X", "PXN", "ILK", "PARVA", "PARVB", "PARVG",
    "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGAV", "ITGB1",
    "ITGB3", "ITGB4", "ITGB5", "ITGB6", "LAMA1", "LAMA2", "LAMA3", "LAMA4",
    "LAMA5", "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2", "NID1", "NID2",
    "HSPG2", "AGRN", "DAG1", "DMD", "UTRN", "SGCA", "SGCB", "SGCG",
    "SGCD", "SGCE", "SGCZ", "SSPN", "DYSF", "CAV3", "TRIM32", "FKRP",
    "FKTN", "POMT1", "POMT2", "POMGNT1", "LARGE1", "B3GNT1", "ISPD", "GMPPB",
    "DPM1", "DPM2", "DPM3", "DOLK", "ALG13", "ALG14", "B4GAT1", "POMK",
    "TMEM5", "RXYLT1", "FKTN", "FKRP", "LARGE2", "DAG1",
    # Note: This is a partial list; the actual gene order comes from the CSV
]
