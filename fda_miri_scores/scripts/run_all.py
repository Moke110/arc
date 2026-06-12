#!/usr/bin/env python3
"""
MIRI Drug Efficacy Scoring — Main Pipeline.

Orchestrates three complementary scoring approaches (3 methods, 4 rank components):

  1. WCS + RGES            — gene-set-based, highest adoption in literature
  2. Zhang Score (ssCMap)  — rank-based, preserves expression magnitude
  3. Pathway-level RGES     — mechanistic interpretability

Output:
  fda_miri_scores_integrated.csv — all scores + integrated ranking (mean of 4 ranks).

Usage:
    python scripts/run_all.py \
        --deg-csv xpert_fda-cm_deg.csv \
        --gene-set reference-genesets/GSE160516_MIRI.json \
        --pathways reference-genesets/MIRI_pathways.json \
        --output-dir ./

Or with defaults (from xpert_fda-miri_scores/ directory):
    python scripts/run_all.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to path
_PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PROJECT_DIR / "scripts"))

from score_wcs_rges import score_all_wcs_rges
from score_zhang import score_all_zhang
from score_pathway_rges import score_all_pathways


# ── Integration ─────────────────────────────────────────────────────

def integrate_scores(wcs_df, zhang_df, pathway_df,
                     output_path=None):
    """
    Merge scores from all three methods and compute integrated ranking.

    Ranking strategy:
      - rank_WCS: lower is better (more negative = stronger reversal)
      - rank_RGES: lower is better
      - rank_Zhang: lower is better (ascending c, most negative first)
      - rank_pathway: DESCENDING (higher mean_efficacy = better)

    Final rank = mean of 4 individual method ranks.
    """
    # Merge on drug identifiers
    merged = wcs_df[["drug_name", "smiles", "drug_idx",
                      "WCS", "RGES", "rank_WCS", "rank_RGES"]].merge(
        zhang_df[["drug_name", "Zhang_c", "Zhang_p_approx", "rank_Zhang"]],
        on="drug_name", how="inner"
    ).merge(
        pathway_df[["drug_name", "mean_efficacy", "n_effective",
                     "pathway_reversal", "rank_pathway"]],
        on="drug_name", how="inner"
    )

    # Integrated rank: mean of 4 individual ranks
    merged["rank_integrated"] = (
        merged["rank_WCS"] + merged["rank_RGES"] +
        merged["rank_Zhang"] + merged["rank_pathway"]
    ) / 4.0

    merged = merged.sort_values("rank_integrated")

    if output_path:
        merged.to_csv(output_path, index=False)
        print(f"Integrated scores saved to {output_path}")

    return merged


# ── Main ────────────────────────────────────────────────────────────

def run_pipeline(deg_csv, gene_set_json, pathways_json, output_dir="scores",
                  zscore=False, rank_ref_json=None):
    """Run all three scoring methods and integrate results (4 rank components)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MIRI Drug Efficacy Scoring Pipeline")
    print("=" * 70)
    print(f"  DEG data:      {deg_csv}")
    print(f"  Gene set:      {gene_set_json}")
    print(f"  Rank ref:      {rank_ref_json}")
    print(f"  Pathways:      {pathways_json}")
    print(f"  Z-score:       {zscore}")
    print(f"  Output dir:    {output_dir}")
    print(f"\n  Methods: WCS + RGES + Zhang + Pathway RGES")
    print()

    # ── Score 1: WCS + RGES ──────────────────────────────────────
    print("\n" + "─" * 50)
    print("SCORE 1: WCS + RGES (gene-set based)")
    print("─" * 50)
    wcs_out = output_dir / "scores_wcs_rges.csv"
    wcs_scores = score_all_wcs_rges(
        deg_csv_path=deg_csv,
        gene_set_path=gene_set_json,
        weight=1.0,
        output_path=wcs_out,
        zscore=zscore,
    )

    # ── Score 2: Zhang Score ─────────────────────────────────────
    print("\n" + "─" * 50)
    print("SCORE 2: Zhang Score (rank-based)")
    print("─" * 50)
    zhang_out = output_dir / "scores_zhang.csv"
    zhang_scores = score_all_zhang(
        deg_csv_path=deg_csv,
        gene_set_path=gene_set_json if gene_set_json and Path(gene_set_json).exists() else None,
        rank_ref_path=rank_ref_json if rank_ref_json and Path(rank_ref_json).exists() else None,
        output_path=zhang_out,
        zscore=zscore,
    )

    # ── Score 3: Pathway-Level RGES ──────────────────────────────
    print("\n" + "─" * 50)
    print("SCORE 3: Pathway-Level RGES")
    print("─" * 50)
    pw_out = output_dir / "scores_pathway_rges.csv"
    pw_scores = score_all_pathways(
        deg_csv_path=deg_csv,
        pathways_json_path=pathways_json,
        weight=1.0,
        output_path=pw_out,
        zscore=zscore,
    )

    # ── Integration ─────────────────────────────────────────────
    print("\n" + "─" * 50)
    print("INTEGRATION: Combined Ranking (mean of 4 ranks)")
    print("─" * 50)
    integrated_out = output_dir / "fda_miri_scores_integrated.csv"
    integrated = integrate_scores(wcs_scores, zhang_scores, pw_scores,
                                  output_path=integrated_out)

    # ── Final Summary ────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("FINAL: Top 20 MIRI Drug Candidates (Integrated Rank)")
    print("=" * 70)
    cols = ["drug_name", "smiles", "WCS", "RGES", "Zhang_c",
            "mean_efficacy", "rank_integrated"]
    available = [c for c in cols if c in integrated.columns]
    print(integrated.head(20)[available].to_string(index=False))

    # Consistency check: drugs ranked top-50 by all methods
    top50_wcs = set(integrated.nsmallest(50, "rank_WCS")["drug_name"])
    top50_rges = set(integrated.nsmallest(50, "rank_RGES")["drug_name"])
    top50_zhang = set(integrated.nsmallest(50, "rank_Zhang")["drug_name"])
    top50_pw = set(integrated.nlargest(50, "mean_efficacy")["drug_name"])

    consensus = top50_wcs & top50_rges & top50_zhang & top50_pw
    print(f"\n  Consensus top-50 (all 4 rank components agree): {len(consensus)} drugs")
    if consensus:
        print(f"  {sorted(consensus)[:10]}...")

    # Method agreement (Spearman between ranks)
    from scipy.stats import spearmanr
    print("\n  Method rank agreement (Spearman ρ):")
    pairs = [("rank_WCS", "rank_RGES"), ("rank_WCS", "rank_Zhang"),
             ("rank_RGES", "rank_Zhang"),
             ("rank_WCS", "rank_pathway"), ("rank_RGES", "rank_pathway")]
    for a, b in pairs:
        rho, p = spearmanr(integrated[a], integrated[b])
        print(f"    {a} vs {b}: ρ = {rho:.3f}")

    return integrated


# ── CLI ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="MIRI Drug Efficacy Scoring — Full Pipeline")
    parser.add_argument("--deg-csv",
                        default=str(_PROJECT_DIR / "xpert_fda-cm_deg.csv"),
                        help="Path to xpert_fda-cm_deg.csv")
    parser.add_argument("--gene-set",
                        default=str(_PROJECT_DIR / "reference-genesets" /
                                    "MIRI_up-down.json"),
                        help="Path to up/down gene set JSON (Type A)")
    parser.add_argument("--rank-ref",
                        default=str(_PROJECT_DIR / "reference-genesets" /
                                    "MIRI_rank.json"),
                        help="Path to MIRI_rank.json (Type B, for Zhang)")
    parser.add_argument("--pathways",
                        default=str(_PROJECT_DIR / "reference-genesets" /
                                    "MIRI_pathways.json"),
                        help="Path to pathway gene sets JSON")
    parser.add_argument("--output-dir", "-o",
                        default=str(_PROJECT_DIR / "scores"),
                        help="Output directory for score CSV files")
    parser.add_argument("--zscore", action="store_true",
                        help="Z-score normalize per gene across drugs "
                             "(recommended: improves drug discrimination)")
    args = parser.parse_args()

    run_pipeline(
        deg_csv=args.deg_csv,
        gene_set_json=args.gene_set,
        pathways_json=args.pathways,
        output_dir=args.output_dir,
        zscore=args.zscore,
        rank_ref_json=args.rank_ref,
    )


if __name__ == "__main__":
    main()
