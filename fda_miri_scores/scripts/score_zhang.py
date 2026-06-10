"""
Score 2: Zhang Score (ssCMap) — Rank-based drug efficacy scoring.

Two reference modes:

  Type A (--gene-set): Binary up/down gene sets (GSE160516_MIRI.json).
      Disease signed ranks approximated from drug |logFC| within the disease set.

  Type B (--rank-ref): Full ranked gene list with disease logFC magnitudes
      (MIRI_rank.json). Disease signed ranks are pre-computed from the actual
      disease experiment, preserving true effect magnitudes. This is the
      preferred mode for Zhang Score.

Method: Zhang & Gant 2008, BMC Bioinformatics. DOI: 10.1186/1471-2105-9-258

Negative c → drug reverses disease signature.

Usage:
    # Type A (up/down sets):
    python scripts/score_zhang.py \\
        --deg-csv xpert_fda-cm_deg.csv \\
        --gene-set reference-genesets/GSE160516_MIRI.json \\
        --output scores_zhang.csv

    # Type B (ranked list, preferred):
    python scripts/score_zhang.py \\
        --deg-csv xpert_fda-cm_deg.csv \\
        --rank-ref reference-genesets/MIRI_rank.json \\
        --output scores_zhang.csv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (signed_rank_by_abs_value, load_drug_deg,
                   load_gene_set, load_rank_reference)


# ── Zhang Score (Type B: pre-computed disease signed ranks) ─────────

def compute_zhang_score_ranked(drug_deg, disease_signed_ranks):
    """
    Zhang Connectivity Score c ∈ [-1, 1] using pre-computed disease ranks.

    This is the preferred Type B implementation:
      - R(g): signed rank of drug gene by |drug_deg| (computed from drug data)
      - s(g): signed rank of disease gene by |disease_logFC| (from MIRI_rank.json)

    Negative c → drug reverses disease signature.

    Parameters
    ----------
    drug_deg : np.ndarray, shape (N,)
    disease_signed_ranks : np.ndarray, shape (N,)
        Pre-computed disease signed ranks. 0 for non-disease genes.

    Returns
    -------
    dict: c, C, C_max, p_approx, n_genes_used
    """
    n = len(drug_deg)

    # Drug gene signed ranks
    R = signed_rank_by_abs_value(drug_deg)

    # Disease gene mask
    disease_mask = disease_signed_ranks != 0
    m = disease_mask.sum()
    if m == 0:
        return {"c": np.nan, "C": 0, "C_max": 1, "p_approx": 1.0, "n_genes": 0}

    # Connection: Σ R(g) · s(g) over all disease genes
    C = np.sum(R[disease_mask] * disease_signed_ranks[disease_mask])

    # Maximum possible: perfect alignment of R with s
    # R ranges from +n to -n, s ranges from +m to -m
    # C_max = Σ |R_best_i| · |s_i| where R_best_i sorts |R| descending
    R_abs = np.abs(R[disease_mask])
    s_abs = np.abs(disease_signed_ranks[disease_mask])
    # Sort both descending — this is the best possible (most negative) alignment
    C_max = np.sum(np.sort(R_abs)[::-1] * np.sort(s_abs)[::-1])
    if C_max == 0:
        C_max = 1.0

    c = C / C_max

    # p-value via normal approximation
    var_R = np.var(R)
    sum_s_sq = np.sum(disease_signed_ranks[disease_mask] ** 2)
    var_C = sum_s_sq * var_R / n if n > 0 else 1.0
    if var_C > 0:
        z = abs(C) / np.sqrt(var_C)
        p_approx = 2 * norm.sf(z)
    else:
        z = 0.0
        p_approx = 1.0

    return {"c": c, "C": C, "C_max": C_max, "p_approx": p_approx, "n_genes": m}


# ── Zhang Score (Type A: up/down gene sets — legacy) ────────────────

def compute_zhang_score_updown(drug_deg, up_idx, down_idx):
    """
    Zhang Connectivity Score using up/down gene sets (Type A reference).

    Disease signed ranks are approximated from |drug_deg| within the disease set.
    Use compute_zhang_score_ranked() for proper disease ranks.
    """
    n = len(drug_deg)
    R = signed_rank_by_abs_value(drug_deg)

    all_disease_idx = np.concatenate([up_idx, down_idx])
    m = len(all_disease_idx)
    if m == 0:
        return {"c": np.nan, "C": 0, "C_max": 1, "p_approx": 1.0, "n_genes": 0}

    disease_abs_vals = np.abs(drug_deg[all_disease_idx])
    disease_abs_rank = rankdata(-disease_abs_vals, method="average")
    disease_signed_rank = np.zeros(n)
    for i, idx in enumerate(all_disease_idx):
        direction = 1 if idx in up_idx else -1
        disease_signed_rank[idx] = direction * (m - disease_abs_rank[i] + 1)

    C = np.sum(R[all_disease_idx] * disease_signed_rank[all_disease_idx])
    C_max = np.sum([(n - i + 1) * (m - i + 1) for i in range(1, m + 1)])
    c = C / C_max if C_max > 0 else 0.0

    var_R = np.var(R)
    sum_s_sq = np.sum(disease_signed_rank[all_disease_idx] ** 2)
    var_C = sum_s_sq * var_R / n if n > 0 else 1.0
    if var_C > 0:
        z = abs(C) / np.sqrt(var_C)
        p_approx = 2 * norm.sf(z)
    else:
        p_approx = 1.0

    return {"c": c, "C": C, "C_max": C_max, "p_approx": p_approx, "n_genes": m}


# ── Batch scoring ───────────────────────────────────────────────────

def score_all_zhang(deg_csv_path, gene_set_path=None, rank_ref_path=None,
                    output_path=None, zscore=False):
    """Batch Zhang Score for all FDA drugs.

    Parameters
    ----------
    gene_set_path : str, optional
        Type A: path to up/down gene set JSON.
    rank_ref_path : str, optional
        Type B: path to ranked reference JSON (MIRI_rank.json). Preferred.
    """
    print(f"Loading DEG from {deg_csv_path} ...")
    if zscore:
        print("  Applying per-gene z-score normalization...")
    meta, gene_cols, deg_matrix = load_drug_deg(deg_csv_path, zscore=zscore)
    n_drugs, n_genes = deg_matrix.shape
    print(f"  {n_drugs} drugs × {n_genes} genes")

    # Determine mode
    if rank_ref_path:
        # ── Type B: pre-computed disease signed ranks ──
        print(f"Loading ranked reference from {rank_ref_path} ...")
        disease_ranks, ref_meta = load_rank_reference(rank_ref_path, gene_cols)
        print(f"  {ref_meta['name']}: {ref_meta['n_up_mapped']} UP + "
              f"{ref_meta['n_down_mapped']} DOWN (mapped to L1000)")

        results = []
        for i in tqdm(range(n_drugs), desc="Scoring Zhang (rank)"):
            r = compute_zhang_score_ranked(deg_matrix[i], disease_ranks)
            results.append({
                "drug_name": meta.iloc[i]["drug_name"],
                "smiles": meta.iloc[i]["smiles"],
                "drug_idx": int(meta.iloc[i]["drug_idx"]),
                "Zhang_c": r["c"],
                "Zhang_C": r["C"],
                "Zhang_C_max": r["C_max"],
                "Zhang_p_approx": r["p_approx"],
                "n_disease_genes": r["n_genes"],
            })

    elif gene_set_path:
        # ── Type A: up/down gene sets ──
        print(f"Loading gene set from {gene_set_path} ...")
        up_idx, down_idx, gs_meta = load_gene_set(gene_set_path, gene_cols)
        print(f"  {gs_meta['name']}: {len(up_idx)} UP, {len(down_idx)} DOWN "
              f"(mapped to L1000)")

        results = []
        for i in tqdm(range(n_drugs), desc="Scoring Zhang (up/down)"):
            r = compute_zhang_score_updown(deg_matrix[i], up_idx, down_idx)
            results.append({
                "drug_name": meta.iloc[i]["drug_name"],
                "smiles": meta.iloc[i]["smiles"],
                "drug_idx": int(meta.iloc[i]["drug_idx"]),
                "Zhang_c": r["c"],
                "Zhang_C": r["C"],
                "Zhang_C_max": r["C_max"],
                "Zhang_p_approx": r["p_approx"],
                "n_disease_genes": r["n_genes"],
            })
    else:
        raise ValueError("Either --gene-set (Type A) or --rank-ref (Type B) required")

    scores = pd.DataFrame(results)
    scores["rank_Zhang"] = scores["Zhang_c"].rank()
    scores = scores.sort_values("Zhang_c")

    if output_path:
        scores.to_csv(output_path, index=False)
        print(f"Saved {len(scores)} drugs to {output_path}")

    n_strong = scores["Zhang_c"] < -0.3
    print(f"\n  Strong candidates (Zhang c < -0.3): {n_strong.sum()}/{len(scores)}")
    print(f"  Top 10:\n{scores.head(10)[['drug_name', 'Zhang_c', 'Zhang_p_approx']].to_string(index=False)}")

    return scores


# ── CLI ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Zhang Score (ssCMap) drug efficacy scoring (rank-based)")
    parser.add_argument("--deg-csv", required=True,
                        help="Path to xpert_fda-cm_deg.csv")
    parser.add_argument("--gene-set", default=None,
                        help="Type A: Path to up/down gene set JSON")
    parser.add_argument("--rank-ref", default=None,
                        help="Type B: Path to MIRI_rank.json (preferred)")
    parser.add_argument("--output", default="scores/scores_zhang.csv",
                        help="Output CSV path")
    parser.add_argument("--zscore", action="store_true",
                        help="Z-score normalize per gene across drugs")
    args = parser.parse_args()

    if not args.gene_set and not args.rank_ref:
        parser.error("Either --gene-set or --rank-ref is required")

    score_all_zhang(
        deg_csv_path=args.deg_csv,
        gene_set_path=args.gene_set,
        rank_ref_path=args.rank_ref,
        output_path=args.output,
        zscore=args.zscore,
    )


if __name__ == "__main__":
    main()
