"""
Scoring module for the xpert-fda-v1 workflow.

Implements the three MIRI efficacy scoring pipelines (WCS+RGES, Zhang/ssCMap,
pathway RGES) and their integration. This is the refactored successor of the
scripts frozen in miri_reposition/xpert-fda-v1_history/scripts/. The scoring
logic is unchanged from the archived pipeline; only the paths are workflow-
root-aware.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from tqdm import tqdm

from utils import load_drug_deg


# ── Core enrichment statistic (GSEA-style, shared by all methods) ──────

def enrichment_score(drug_sorted_deg, hit_mask, weight=1.0):
    """Weighted Kolmogorov-Smirnov enrichment statistic."""
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


# ── Data I/O ───────────────────────────────────────────────────────────

def _gene_to_idx(gene_order):
    """Map L1000 gene symbols to their column index (shared by all loaders)."""
    return {g: i for i, g in enumerate(gene_order)}


def _zscore(series):
    """Standardize a pandas Series across rows."""
    return (series - series.mean()) / series.std()


def load_gene_set(json_path, gene_order):
    with open(json_path) as f:
        data = json.load(f)
    gene_to_idx = _gene_to_idx(gene_order)
    up_idx = np.array([gene_to_idx[g] for g in data.get("up_genes", []) if g in gene_to_idx], dtype=int)
    down_idx = np.array([gene_to_idx[g] for g in data.get("down_genes", []) if g in gene_to_idx], dtype=int)
    return up_idx, down_idx


def load_rank_reference(json_path, gene_order):
    with open(json_path) as f:
        data = json.load(f)
    gene_to_idx = _gene_to_idx(gene_order)
    signed_ranks = np.zeros(len(gene_order), dtype=float)
    for entry in data.get("genes", []):
        if entry["gene"] in gene_to_idx:
            signed_ranks[gene_to_idx[entry["gene"]]] = entry["signed_rank"]
    return signed_ranks


def load_pathway_gene_sets(json_path, gene_order, min_genes=6):
    with open(json_path) as f:
        data = json.load(f)
    gene_to_idx = _gene_to_idx(gene_order)
    pathways = []
    for pw in data.get("pathways", []):
        genes = pw.get("genes", [])
        indices = np.array([gene_to_idx[g] for g in genes if g in gene_to_idx], dtype=int)
        if len(indices) >= min_genes:
            pathways.append({
                "name": pw["name"],
                "direction": pw.get("direction", "UP"),
                "gene_indices": indices,
            })
    return pathways


def load_critical_names(pathways_json):
    """Names of the curated critical MIRI pathways (shared by scoring + filtering)."""
    with open(pathways_json) as f:
        data = json.load(f)
    return [p["name"] for p in data["pathways"] if p.get("critical")]


def signed_rank_by_abs_value(values):
    """Assign signed ranks by absolute value (descending)."""
    n = len(values)
    abs_rank = rankdata(-np.abs(values), method="average")
    signed = np.sign(values) * (n - abs_rank + 1)
    return np.where(values == 0, 0, signed)


# ── Score 1: WCS + RGES ───────────────────────────────────────────────

def compute_wcs(drug_deg, up_idx, down_idx, weight=1.0):
    n_genes = len(drug_deg)
    sort_idx = np.argsort(-drug_deg)
    sorted_deg = drug_deg[sort_idx]
    up_mask = np.zeros(n_genes, dtype=bool); up_mask[up_idx] = True
    down_mask = np.zeros(n_genes, dtype=bool); down_mask[down_idx] = True
    ES_up = enrichment_score(sorted_deg, up_mask[sort_idx], weight)
    ES_down = enrichment_score(sorted_deg, down_mask[sort_idx], weight)
    if np.sign(ES_up) != np.sign(ES_down):
        WCS = (ES_up - ES_down) / 2.0
    else:
        WCS = 0.0
    return {"WCS": WCS, "ES_up": ES_up, "ES_down": ES_down}


def score_all_wcs_rges(deg_csv, gene_set, zscore=False):
    meta, gene_cols, deg_matrix = load_drug_deg(deg_csv, zscore=zscore)
    up_idx, down_idx = load_gene_set(gene_set, gene_cols)
    rows = []
    for i in tqdm(range(deg_matrix.shape[0]), desc="WCS+RGES"):
        r = compute_wcs(deg_matrix[i], up_idx, down_idx)
        rows.append({
            "drug_name": meta.iloc[i]["drug_name"], "smiles": meta.iloc[i]["smiles"],
            "drug_idx": int(meta.iloc[i]["drug_idx"]),
            "WCS": r["WCS"], "ES_up": r["ES_up"], "ES_down": r["ES_down"],
            "RGES": abs(r["ES_up"]) - abs(r["ES_down"]),
        })
    s = pd.DataFrame(rows)
    s["WCS_zs"] = _zscore(s["WCS"])
    s["RGES_zs"] = _zscore(s["RGES"])
    return s


# ── Score 2: Zhang / ssCMap ───────────────────────────────────────────

def compute_zhang_ranked(drug_deg, disease_signed_ranks):
    n = len(drug_deg)
    R = signed_rank_by_abs_value(drug_deg)
    mask = disease_signed_ranks != 0
    m = mask.sum()
    if m == 0:
        return {"c": np.nan, "p_approx": 1.0}
    C = np.sum(R[mask] * disease_signed_ranks[mask])
    R_abs = np.abs(R[mask]); s_abs = np.abs(disease_signed_ranks[mask])
    C_max = np.sum(np.sort(R_abs)[::-1] * np.sort(s_abs)[::-1])
    if C_max == 0:
        C_max = 1.0
    c = C / C_max
    var_R = np.var(R); sum_s_sq = np.sum(disease_signed_ranks[mask] ** 2)
    var_C = sum_s_sq * var_R / n if n > 0 else 1.0
    p = 2 * norm.sf(abs(C) / np.sqrt(var_C)) if var_C > 0 else 1.0
    return {"c": c, "p_approx": p}


def score_all_zhang(deg_csv, rank_ref, zscore=False):
    meta, gene_cols, deg_matrix = load_drug_deg(deg_csv, zscore=zscore)
    disease_ranks = load_rank_reference(rank_ref, gene_cols)
    rows = []
    for i in tqdm(range(deg_matrix.shape[0]), desc="Zhang"):
        r = compute_zhang_ranked(deg_matrix[i], disease_ranks)
        rows.append({
            "drug_name": meta.iloc[i]["drug_name"], "smiles": meta.iloc[i]["smiles"],
            "drug_idx": int(meta.iloc[i]["drug_idx"]),
            "Zhang_c": r["c"], "Zhang_p_approx": r["p_approx"],
        })
    s = pd.DataFrame(rows)
    s["Zhang_zs"] = _zscore(s["Zhang_c"])
    return s


# ── Score 3: Pathway RGES ─────────────────────────────────────────────

def compute_pathway_rges(drug_deg, gene_indices, direction, weight=1.0):
    n_genes = len(drug_deg)
    sort_idx = np.argsort(-drug_deg)
    sorted_deg = drug_deg[sort_idx]
    hit_mask = np.zeros(n_genes, dtype=bool)
    hit_mask[gene_indices] = True
    ES = enrichment_score(sorted_deg, hit_mask[sort_idx], weight)
    efficacy = -ES if direction == "UP" else ES
    return {"ES": ES, "efficacy": efficacy}


def score_all_pathways(deg_csv, pathways_json, zscore=False):
    meta, gene_cols, deg_matrix = load_drug_deg(deg_csv, zscore=zscore)
    pathways = load_pathway_gene_sets(pathways_json, gene_cols)
    with open(pathways_json) as f:
        pw_meta = json.load(f)
    critical = [p["name"] for p in pw_meta["pathways"] if p.get("critical")]

    rows = []
    for i in tqdm(range(deg_matrix.shape[0]), desc="Pathways"):
        row = {"drug_name": meta.iloc[i]["drug_name"], "smiles": meta.iloc[i]["smiles"],
               "drug_idx": int(meta.iloc[i]["drug_idx"])}
        total, n_eff = 0.0, 0
        for pw in pathways:
            r = compute_pathway_rges(deg_matrix[i], pw["gene_indices"], pw["direction"])
            total += r["efficacy"]
            if r["efficacy"] > 0.1:
                n_eff += 1
            row[f"{pw['name']}_efficacy"] = r["efficacy"]
        row["total_efficacy"] = total
        row["mean_efficacy"] = total / len(pathways) if pathways else 0.0
        row["n_effective"] = n_eff
        row["pathway_reversal"] = n_eff / len(pathways) if pathways else 0.0
        rows.append(row)
    s = pd.DataFrame(rows)
    s["pathway_reversal_zs"] = _zscore(s["pathway_reversal"])
    for name in critical:
        col = f"{name}_efficacy"
        if col in s.columns:
            mu, sigma = s[col].mean(), s[col].std()
            s[f"{col}_zs"] = (s[col] - mu) / sigma if sigma > 1e-10 else 0.0
    return s, critical


# ── Integration (the 6-step screening pipeline's step 1) ──────────────

def merge_components(wcs, zhang, pw, critical):
    """Merge the three scorers' outputs into one per-drug components table."""
    scores = wcs[["drug_name", "smiles", "drug_idx", "WCS", "RGES",
                  "WCS_zs", "RGES_zs"]].merge(
        zhang[["drug_name", "Zhang_c", "Zhang_p_approx", "Zhang_zs"]],
        on="drug_name", how="inner").merge(
        pw[["drug_name", "mean_efficacy", "n_effective", "pathway_reversal",
            "pathway_reversal_zs"]], on="drug_name", how="inner")

    pw_extra_cols = ["drug_name"]
    for n in critical:
        pw_extra_cols.append(f"{n}_efficacy")
        if f"{n}_efficacy_zs" in pw.columns:
            pw_extra_cols.append(f"{n}_efficacy_zs")
    return scores.merge(pw[pw_extra_cols], on="drug_name", how="inner")


def add_integrated(scores, critical):
    """Compute the integrated score from the z-component columns.

    integrated = -0.3*RGES_zs - 0.2*Zhang_zs + 0.3*pathway_reversal_zs
                 + 0.025*sum(critical efficacy zs)
    """
    scores = scores.copy()
    crit_sum = sum(scores.get(f"{n}_efficacy_zs", pd.Series(0, index=scores.index))
                   for n in critical)
    scores["integrated"] = (
        -0.3 * scores["RGES_zs"] - 0.2 * scores["Zhang_zs"]
        + 0.3 * scores["pathway_reversal_zs"] + 0.025 * crit_sum
    )
    scores["integrated_zs"] = ((scores["integrated"] - scores["integrated"].mean())
                               / scores["integrated"].std())
    return scores


def integrate_scores(wcs, zhang, pw, critical, zscore=True):
    """Combine the four z-components into an integrated score."""
    scores = add_integrated(merge_components(wcs, zhang, pw, critical), critical)
    return scores.sort_values("integrated", ascending=False).reset_index(drop=True)


def run_scoring(deg_csv, gene_set, rank_ref, pathways_json, zscore=True):
    """Run all three scorers and return the integrated 1-scores table."""
    wcs = score_all_wcs_rges(deg_csv, gene_set, zscore=zscore)
    zhang = score_all_zhang(deg_csv, rank_ref, zscore=zscore)
    pw, critical = score_all_pathways(deg_csv, pathways_json, zscore=zscore)
    return integrate_scores(wcs, zhang, pw, critical, zscore=zscore)
