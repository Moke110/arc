#!/usr/bin/env python3
"""
Step 08 — MIRI reference gene sets + pathways.

Consumes the per-GSE DEG lists in disease_data/{GSE}/local/deg/ and produces the
integrated MIRI reference DEG list used by the efficacy scoring step:

  08_miri_phenotype/output/MIRI_up-down.json — union UP/DOWN gene set (majority vote)
  08_miri_phenotype/output/MIRI_rank.json    — ranked reference for Zhang scoring

input/MIRI_pathways.json is static (curated) and is not rewritten by this step.
This is the regenerable, cross-dataset reference; the per-GSE lists under
disease_data/ are shared inputs reused by future workflows.

Usage:
    python 08_miri_phenotype/miri_phenotype.py
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, DISEASE_DATA_ROOT, REFERENCE_ROOT, XPERT_METADATA  # noqa: E402

REF_DIR = WORKFLOW_ROOT / "08_miri_phenotype" / "output"
PRIMARY_GSE = "GSE160516"          # mouse MIRI (largest L1000 overlap)

# Display order + labels for the union metadata (matches the frozen reference).
GSE_ORDER = ["GSE160516", "GSE108940", "GSE137293"]
GSE_LABELS = {
    "GSE160516": "mouse MIRI",
    "GSE108940": "rat MIRI",
    "GSE137293": "pig scar vs border",
}
PRIMARY_PLATFORM = "Affymetrix Clariom S"


GSE_DEG_LISTS = {
    "GSE160516": DISEASE_DATA_ROOT / "GSE160516" / "local" / "deg" / "GSE160516_MIRI.json",
    # The formal v2 screen predates the shared source correction in 7d1d0dcc.
    "GSE108940": REFERENCE_ROOT / "GSE108940_MIRI_v2.json",
    "GSE137293": DISEASE_DATA_ROOT / "GSE137293" / "local" / "deg" / "GSE137293_scar_vs_border.json",
}


def load_per_gse_deg_list(gse: str) -> dict:
    """Load the exact source for one series in the accepted v2 screen."""
    with GSE_DEG_LISTS[gse].open() as f:
        return json.load(f)


def build_union(per_gse: dict) -> dict:
    """Union of UP/DOWN genes across species with majority-vote conflict resolution."""
    up_votes, down_votes = Counter(), Counter()
    source_up, source_down = {}, {}
    for gse, gs in per_gse.items():
        for g in gs["up_genes"]:
            up_votes[g] += 1
            source_up.setdefault(g, []).append(gse)
        for g in gs["down_genes"]:
            down_votes[g] += 1
            source_down.setdefault(g, []).append(gse)

    conflicts = set(up_votes) & set(down_votes)
    resolved_up, resolved_down = set(), set()
    details = []
    for g in sorted(conflicts):
        if up_votes[g] > down_votes[g]:
            resolved_up.add(g)
        elif down_votes[g] > up_votes[g]:
            resolved_down.add(g)
        details.append(
            f"{g}: UP={up_votes[g]} ({source_up.get(g, [])}), "
            f"DOWN={down_votes[g]} ({source_down.get(g, [])}) → "
            f"{'UP' if up_votes[g] > down_votes[g] else 'DOWN' if down_votes[g] > up_votes[g] else 'REMOVED'}"
        )

    resolved_up |= (set(up_votes) - conflicts)
    resolved_down |= (set(down_votes) - conflicts)

    def support(g, direction):
        return sum(1 for gse, gs in per_gse.items()
                   if g in gs["down_genes" if direction == "DOWN" else "up_genes"])

    up_sorted = sorted(resolved_up, key=lambda g: (-support(g, "UP"), g))
    down_sorted = sorted(resolved_down, key=lambda g: (-support(g, "DOWN"), g))

    contributions = ", ".join(
        f"{gse} ({GSE_LABELS.get(gse, gse)}, "
        f"{len(per_gse[gse]['up_genes'])}+{len(per_gse[gse]['down_genes'])})"
        for gse in GSE_ORDER if gse in per_gse
    )
    gene_set = {
        "name": "MIRI_up-down",
        "species": "Multi-species union (Mus musculus + Rattus norvegicus + Sus scrofa)",
        "source_geo": " + ".join(gse for gse in GSE_ORDER if gse in per_gse),
        "comparison": "MIRI/IR vs Sham/Control (union across species)",
        "method": "union",
        "description": (
            f"Union of MIRI dysregulated genes from {len(per_gse)} species mapped to human L1000. "
            "Genes are included if they appear as UP or DOWN in ANY dataset. "
            "Conflicts (gene UP in one dataset, DOWN in another) resolved by majority vote; "
            "ties removed. "
            f"Dataset contributions: {contributions}."
        ),
        "n_datasets": len(per_gse),
        "n_conflicts": len(conflicts),
        "conflict_resolution": "majority_vote",
        "_conflict_details": details,
        "up_genes": up_sorted,
        "down_genes": down_sorted,
        "n_up_l1000": len(up_sorted),
        "n_down_l1000": len(down_sorted),
        "up_by_support": {
            "3_datasets": [g for g in up_sorted if support(g, "UP") == 3],
            "2_datasets": [g for g in up_sorted if support(g, "UP") == 2],
            "1_dataset": [g for g in up_sorted if support(g, "UP") == 1],
        },
        "down_by_support": {
            "3_datasets": [g for g in down_sorted if support(g, "DOWN") == 3],
            "2_datasets": [g for g in down_sorted if support(g, "DOWN") == 2],
            "1_dataset": [g for g in down_sorted if support(g, "DOWN") == 1],
        },
    }
    return gene_set


def build_rank_reference(per_gse: dict, l1000_meta: pd.DataFrame) -> dict:
    """Ranked reference for Zhang scoring from the primary (mouse) GSE.

    Uses the full DEG results of the primary GSE, filtered to human L1000
    landmarks, with signed ranks by |disease logFC|.
    """
    deg = pd.read_csv(DISEASE_DATA_ROOT / PRIMARY_GSE / "source"
                      / f"{PRIMARY_GSE}_deg_results.csv")
    up_l1000 = pd.read_csv(DISEASE_DATA_ROOT / PRIMARY_GSE / "local" / "deg"
                           / f"{PRIMARY_GSE}_up_l1000.csv")
    down_l1000 = pd.read_csv(DISEASE_DATA_ROOT / PRIMARY_GSE / "local" / "deg"
                             / f"{PRIMARY_GSE}_down_l1000.csv")

    landmark = set(l1000_meta["gene_symbol"])
    mouse_to_human, mouse_logfc = {}, {}
    for df in (up_l1000, down_l1000):
        for _, row in df.iterrows():
            mouse_to_human[row.iloc[0]] = row["human_gene"]
            mouse_logfc[row.iloc[0]] = row["logFC"]

    mapped = []
    for _, row in deg.iterrows():
        mouse_gene = row["Gene"]
        if mouse_gene not in mouse_to_human or mouse_to_human[mouse_gene] not in landmark:
            continue
        mapped.append({
            "gene": mouse_to_human[mouse_gene],
            "mouse_gene": mouse_gene,
            "logFC": row["logFC"],
            "P_Value": row["P.Value"],
            "adj_P_Val": row["adj.P.Val"],
        })
    if not mapped:
        raise RuntimeError(f"No {PRIMARY_GSE} genes map to L1000!")

    mapped.sort(key=lambda g: abs(g["logFC"]), reverse=True)
    up_genes = [g for g in mapped if g["logFC"] > 0]
    down_genes = [g for g in mapped if g["logFC"] < 0]

    for group in (up_genes, down_genes):
        if not group:
            continue
        abs_vals = np.array([abs(g["logFC"]) for g in group])
        ranks = rankdata(-abs_vals, method="average")
        sign = 1 if group is up_genes else -1
        for g, r in zip(group, ranks):
            g["signed_rank"] = int(sign * (len(group) - r + 1))
            g["direction"] = "UP" if sign > 0 else "DOWN"

    all_ranked = sorted(up_genes + down_genes, key=lambda g: abs(g["logFC"]), reverse=True)
    return {
        "name": "MIRI_rank",
        "type": "ranked_gene_list",
        "species": "Mus musculus",
        "source_geo": PRIMARY_GSE,
        "platform": PRIMARY_PLATFORM,
        "comparison": "MIRI (pooled IR) vs Sham",
        "description": (
            f"Full {PRIMARY_GSE} MIRI DEG results filtered to human L1000 978 landmark genes. "
            "Genes ranked by |disease logFC| (largest effect first). "
            "Each gene has a signed rank: UP genes get +m (largest |logFC|) to +1 (smallest), "
            "DOWN genes get -m (largest |logFC|) to -1 (smallest). "
            "Used as Type B reference for Zhang Score (ssCMap) computation, "
            "preserving disease effect magnitude information."
        ),
        "n_genes_total": len(all_ranked),
        "n_up": len(up_genes),
        "n_down": len(down_genes),
        "genes": all_ranked,
    }


def main():
    parser = argparse.ArgumentParser(description="Build the integrated MIRI reference DEG list")
    parser.add_argument("--l1000-meta",
                        default=str(XPERT_METADATA / "l1000_gene_info_978.csv"),
                        help="L1000 landmark gene metadata CSV")
    args = parser.parse_args()

    l1000_meta = pd.read_csv(args.l1000_meta)
    per_gse = {gse: load_per_gse_deg_list(gse) for gse in GSE_ORDER}
    for gse, genes in per_gse.items():
        print(f"  {gse}: {len(genes['up_genes'])} UP + {len(genes['down_genes'])} DOWN")

    REF_DIR.mkdir(parents=True, exist_ok=True)

    union = build_union(per_gse)
    with open(REF_DIR / "MIRI_up-down.json", "w") as f:
        json.dump(union, f, indent=2, ensure_ascii=False)
    print(f"\nUnion UP: {union['n_up_l1000']}, DOWN: {union['n_down_l1000']}, "
          f"conflicts: {union['n_conflicts']}")
    print(f"  Saved: {REF_DIR / 'MIRI_up-down.json'}")

    rank = build_rank_reference(per_gse, l1000_meta)
    with open(REF_DIR / "MIRI_rank.json", "w") as f:
        json.dump(rank, f, indent=2, ensure_ascii=False)
    print(f"\nRank reference: {rank['n_genes_total']} genes ({rank['n_up']} UP + {rank['n_down']} DOWN)")
    print(f"  Saved: {REF_DIR / 'MIRI_rank.json'}")


if __name__ == "__main__":
    main()
