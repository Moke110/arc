"""
Score 3: Pathway-Level RGES — Per-pathway drug efficacy scoring.

Reference: MIRI pathway gene sets with direction labels (Type C reference).

For each MIRI-relevant pathway with a known disease direction (UP/DOWN),
computes RGES independently, then aggregates into a pathway-level efficacy
profile per drug. This provides mechanistic interpretability: which pathways
does each drug reverse, and how strongly?

Usage:
    python scripts/score_pathway_rges.py \
        --deg-csv xpert_fda-cm_deg.csv \
        --pathways reference-genesets/MIRI_pathways.json \
        --output scores_pathway_rges.csv \
        --weight 1.0
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import enrichment_score, load_drug_deg, load_pathway_gene_sets


# ── Per-pathway RGES ────────────────────────────────────────────────

def compute_pathway_rges(drug_deg, pathway_gene_indices, pathway_direction,
                         weight=1.0):
    """
    Compute RGES for a single pathway.

    Pathway direction interpretation:
      - "UP" in disease: reversal means drug should DOWN-regulate these genes
                         → ES should be negative (enriched at bottom)
                         → efficacy = -ES (positive = good reversal)
      - "DOWN" in disease: reversal means drug should UP-regulate these genes
                           → ES should be positive (enriched at top)
                           → efficacy = ES (positive = good reversal)

    Parameters
    ----------
    drug_deg : np.ndarray, shape (N,)
    pathway_gene_indices : np.ndarray (int)
        Indices of pathway genes in the N-gene array.
    pathway_direction : str
        "UP" or "DOWN" — direction of this pathway in MIRI.
    weight : float
        GSEA weight exponent.

    Returns
    -------
    dict: pathway_name, direction, ES, efficacy (RGES for this pathway),
          n_genes (mapped)
    """
    n_genes = len(drug_deg)
    sort_idx = np.argsort(-drug_deg)
    sorted_deg = drug_deg[sort_idx]

    hit_mask = np.zeros(n_genes, dtype=bool)
    hit_mask[pathway_gene_indices] = True
    ES = enrichment_score(sorted_deg, hit_mask[sort_idx], weight)

    # Efficacy: how well does drug reverse this pathway?
    if pathway_direction == "UP":
        # Disease UP → drug should DOWN-regulate → ES < 0 = good
        # RGES-style: |ES_down| - ... but here we only have one gene set
        # We use -ES: positive efficacy means good reversal
        efficacy = -ES
    else:  # DOWN
        # Disease DOWN → drug should UP-regulate → ES > 0 = good
        efficacy = ES

    return {
        "ES": ES,
        "efficacy": efficacy,
        "n_genes": len(pathway_gene_indices),
    }


# ── Batch scoring ───────────────────────────────────────────────────

def score_all_pathways(deg_csv_path, pathways_json_path, weight=1.0,
                       output_path=None, zscore=False):
    """Batch pathway-level RGES for all FDA drugs."""
    print(f"Loading DEG from {deg_csv_path} ...")
    if zscore:
        print("  Applying per-gene z-score normalization...")
    meta, gene_cols, deg_matrix = load_drug_deg(deg_csv_path, zscore=zscore)
    n_drugs, n_genes = deg_matrix.shape
    print(f"  {n_drugs} drugs × {n_genes} genes")

    print(f"Loading pathways from {pathways_json_path} ...")
    pathways, pw_meta = load_pathway_gene_sets(pathways_json_path, gene_cols)
    print(f"  {pw_meta['name']}: {pw_meta['n_pathways_valid']}/{pw_meta['n_pathways_total']} "
          f"pathways with ≥3 mapped genes")
    for pw in pathways:
        print(f"    {pw['name']} ({pw['direction']}): {pw['n_mapped']}/{pw['n_total']} genes")

    # Per-drug, per-pathway scores
    all_drug_scores = []

    for i in tqdm(range(n_drugs), desc="Scoring pathways"):
        drug_deg = deg_matrix[i]
        drug_name = meta.iloc[i]["drug_name"]

        pw_scores = {}
        total_efficacy = 0.0
        n_effective = 0

        for pw in pathways:
            r = compute_pathway_rges(
                drug_deg, pw["gene_indices"], pw["direction"], weight
            )
            pw_scores[pw["name"]] = r

            # Count pathways where drug shows reversal (efficacy > 0)
            if r["efficacy"] > 0.1:
                n_effective += 1
            total_efficacy += r["efficacy"]

        # Build summary row
        row = {
            "drug_name": drug_name,
            "smiles": meta.iloc[i]["smiles"],
            "drug_idx": int(meta.iloc[i]["drug_idx"]),
            "total_efficacy": total_efficacy,
            "mean_efficacy": total_efficacy / len(pathways) if pathways else 0.0,
            "n_pathways": len(pathways),
            "n_effective": n_effective,
            "pathway_reversal": n_effective / len(pathways) if pathways else 0.0,
        }

        # Add per-pathway efficacy and ES
        for pw in pathways:
            row[f"{pw['name']}_ES"] = pw_scores[pw["name"]]["ES"]
            row[f"{pw['name']}_efficacy"] = pw_scores[pw["name"]]["efficacy"]

        all_drug_scores.append(row)

    scores = pd.DataFrame(all_drug_scores)

    # Rank by mean_efficacy (higher = more pathways reversed)
    scores["rank_pathway"] = scores["mean_efficacy"].rank(ascending=False)
    # Z-score normalization across drugs (pathway_reversal = fraction of pathways reversed)
    scores["pathway_reversal_zs"] = ((scores["pathway_reversal"] - scores["pathway_reversal"].mean())
                                     / scores["pathway_reversal"].std())

    # Compute Z-scores for critical pathway efficacies across drugs
    import json
    with open(pathways_json_path) as f:
        _pw_meta = json.load(f)
    _critical_names = [p["name"] for p in _pw_meta["pathways"] if p.get("critical")]
    for _name in _critical_names:
        _col = f"{_name}_efficacy"
        if _col in scores.columns:
            _mu = scores[_col].mean()
            _sigma = scores[_col].std()
            if _sigma > 1e-10:
                scores[f"{_col}_zs"] = (scores[_col] - _mu) / _sigma
            else:
                scores[f"{_col}_zs"] = 0.0

    scores = scores.sort_values("mean_efficacy", ascending=False)

    if output_path:
        scores.to_csv(output_path, index=False)
        print(f"Saved {len(scores)} drugs to {output_path}")

    print(f"\n  Top 10 by mean pathway efficacy:")
    top_cols = ["drug_name", "mean_efficacy", "n_effective", "pathway_reversal"]
    print(scores.head(10)[top_cols].to_string(index=False))

    return scores


# ── CLI ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Pathway-Level RGES drug efficacy scoring")
    parser.add_argument("--deg-csv", required=True,
                        help="Path to xpert_fda-cm_deg.csv")
    parser.add_argument("--pathways", required=True,
                        help="Path to MIRI pathway gene sets JSON")
    parser.add_argument("--output", default="scores/scores_pathway_rges.csv",
                        help="Output CSV path")
    parser.add_argument("--weight", type=float, default=1.0,
                        help="GSEA weight exponent")
    parser.add_argument("--zscore", action="store_true",
                        help="Z-score normalize per gene across drugs")
    args = parser.parse_args()

    score_all_pathways(
        deg_csv_path=args.deg_csv,
        pathways_json_path=args.pathways,
        weight=args.weight,
        output_path=args.output,
        zscore=args.zscore,
    )


if __name__ == "__main__":
    main()
