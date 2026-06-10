"""
Score 1: WCS + RGES — Gene-set-based drug efficacy scoring.

Reference gene sets: MIRI up/down gene sets (Type A reference).

Methods:
  - WCS (Weighted Connectivity Score)  — Subramanian et al. 2017 Cell
  - RGES (Reverse Gene Expression Score) — Chen et al. 2017 Nat. Commun.

Both share the weighted KS enrichment statistic from utils.enrichment_score().

Usage:
    python scripts/score_wcs_rges.py \
        --deg-csv xpert_fda-cm_deg.csv \
        --gene-set reference-genesets/GSE160516_MIRI.json \
        --output scores_wcs_rges.csv \
        --weight 1.0
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# Allow running from xpert_fda-miri_scores/ directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import enrichment_score, load_drug_deg, load_gene_set


# ── WCS ────────────────────────────────────────────────────────────

def compute_wcs(drug_deg, up_idx, down_idx, weight=1.0):
    """
    Weighted Connectivity Score.

    Parameters
    ----------
    drug_deg : np.ndarray, shape (978,)
        Predicted drug logFC for 978 landmark genes.
    up_idx, down_idx : np.ndarray
        Indices of MIRI up/down-regulated genes in the 978-gene array.
    weight : float
        GSEA weight exponent.

    Returns
    -------
    dict with keys: WCS, ES_up, ES_down
    """
    n_genes = len(drug_deg)
    sort_idx = np.argsort(-drug_deg)  # descending logFC
    sorted_deg = drug_deg[sort_idx]

    up_mask = np.zeros(n_genes, dtype=bool)
    up_mask[up_idx] = True
    down_mask = np.zeros(n_genes, dtype=bool)
    down_mask[down_idx] = True

    ES_up = enrichment_score(sorted_deg, up_mask[sort_idx], weight)
    ES_down = enrichment_score(sorted_deg, down_mask[sort_idx], weight)

    if np.sign(ES_up) != np.sign(ES_down):
        WCS = (ES_up - ES_down) / 2.0
    else:
        WCS = 0.0

    return {"WCS": WCS, "ES_up": ES_up, "ES_down": ES_down}


# ── RGES ───────────────────────────────────────────────────────────

def compute_rges(drug_deg, up_idx, down_idx, weight=1.0):
    """
    Reverse Gene Expression Score.

    RGES = |ES_up| - |ES_down|

    Negative RGES → drug reverses disease signature (therapeutic signal).
    """
    wcs_result = compute_wcs(drug_deg, up_idx, down_idx, weight)
    RGES = abs(wcs_result["ES_up"]) - abs(wcs_result["ES_down"])
    return {"RGES": RGES, **wcs_result}


# ── Batch scoring ───────────────────────────────────────────────────

def score_all_wcs_rges(deg_csv_path, gene_set_path, weight=1.0,
                       output_path=None, zscore=False):
    """Batch WCS + RGES scoring for all FDA drugs."""
    print(f"Loading DEG from {deg_csv_path} ...")
    if zscore:
        print("  Applying per-gene z-score normalization...")
    meta, gene_cols, deg_matrix = load_drug_deg(deg_csv_path, zscore=zscore)
    n_drugs, n_genes = deg_matrix.shape
    print(f"  {n_drugs} drugs × {n_genes} genes")

    print(f"Loading gene set from {gene_set_path} ...")
    up_idx, down_idx, gs_meta = load_gene_set(gene_set_path, gene_cols)
    print(f"  {gs_meta['name']}: {len(up_idx)} UP, {len(down_idx)} DOWN "
          f"(mapped to L1000)")

    results = []
    for i in tqdm(range(n_drugs), desc="Scoring WCS+RGES"):
        drug_deg = deg_matrix[i]
        r = compute_rges(drug_deg, up_idx, down_idx, weight)
        results.append({
            "drug_name": meta.iloc[i]["drug_name"],
            "smiles": meta.iloc[i]["smiles"],
            "drug_idx": int(meta.iloc[i]["drug_idx"]),
            "WCS": r["WCS"],
            "ES_up": r["ES_up"],
            "ES_down": r["ES_down"],
            "RGES": r["RGES"],
        })

    scores = pd.DataFrame(results)
    # Rank: lower WCS/RGES = better (more negative = stronger reversal)
    scores["rank_WCS"] = scores["WCS"].rank()
    scores["rank_RGES"] = scores["RGES"].rank()
    scores["avg_rank"] = (scores["rank_WCS"] + scores["rank_RGES"]) / 2.0
    scores = scores.sort_values("avg_rank")

    if output_path:
        scores.to_csv(output_path, index=False)
        print(f"Saved {len(scores)} drugs to {output_path}")

    # Summary
    n_strong = (scores["WCS"] < -0.2) & (scores["RGES"] < -0.2)
    print(f"\n  Strong candidates (WCS < -0.2 AND RGES < -0.2): {n_strong.sum()}/{len(scores)}")
    print(f"  Top 10:\n{scores.head(10)[['drug_name', 'WCS', 'RGES', 'avg_rank']].to_string(index=False)}")

    return scores


# ── CLI ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="WCS + RGES drug efficacy scoring (gene-set based)")
    parser.add_argument("--deg-csv", required=True,
                        help="Path to xpert_fda-cm_deg.csv")
    parser.add_argument("--gene-set", required=True,
                        help="Path to reference gene set JSON")
    parser.add_argument("--output", default="scores/scores_wcs_rges.csv",
                        help="Output CSV path")
    parser.add_argument("--weight", type=float, default=1.0,
                        help="GSEA weight exponent (0=KS, 1=standard)")
    parser.add_argument("--zscore", action="store_true",
                        help="Z-score normalize per gene across drugs")
    args = parser.parse_args()

    score_all_wcs_rges(
        deg_csv_path=args.deg_csv,
        gene_set_path=args.gene_set,
        weight=args.weight,
        output_path=args.output,
        zscore=args.zscore,
    )


if __name__ == "__main__":
    main()
