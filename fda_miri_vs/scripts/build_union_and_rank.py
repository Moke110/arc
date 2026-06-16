#!/usr/bin/env python3
"""
Build union-based MIRI reference gene set and ranked reference for Zhang scoring.

1. MIRI_up-down.json — union of UP/DOWN genes from all 3 species (resolve conflicts)
2. Delete MIRI_consensus.json  — replaced by union approach
3. MIRI_rank.json — GSE160516 genes ranked by |disease logFC|, with signed ranks
"""

import json
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import rankdata

PROJECT_DIR = Path(__file__).resolve().parents[1]
REF_DIR = PROJECT_DIR / "reference-genesets"
GEO_DIR = PROJECT_DIR / "geo"

# ── L1000 gene set ──────────────────────────────────────────────────

def load_l1000_set():
    df = pd.read_csv("/home/chang/arc/xpert/metadata/l1000_gene_info_978.csv")
    return set(df["gene_symbol"].dropna().unique())


# ═══════════════════════════════════════════════════════════════════
# Task 1: Create union-based MIRI up/down gene set
# ═══════════════════════════════════════════════════════════════════

def build_union(l1000_set):
    """Union of UP/DOWN genes from all datasets, with conflict resolution."""

    datasets = {
        "GSE160516": {
            "up": REF_DIR / "GSE160516_MIRI.json",
            "down": REF_DIR / "GSE160516_MIRI.json",
        },
        "GSE108940": {
            "up": REF_DIR / "GSE108940_MIRI.json",
            "down": REF_DIR / "GSE108940_MIRI.json",
        },
        "GSE137293": {
            "up": REF_DIR / "GSE137293_scar_vs_border.json",
            "down": REF_DIR / "GSE137293_scar_vs_border.json",
        },
    }

    # Collect votes: gene → {UP: count, DOWN: count}
    up_votes = Counter()
    down_votes = Counter()
    dataset_up = {}
    dataset_down = {}
    source_up = {}  # gene → which datasets contributed UP
    source_down = {}

    for name, paths in datasets.items():
        try:
            with open(paths["up"]) as f:
                gs = json.load(f)
        except FileNotFoundError:
            print(f"  Skipping {name} — JSON not found")
            continue

        for g in gs["up_genes"]:
            up_votes[g] += 1
            source_up.setdefault(g, []).append(name)
        for g in gs["down_genes"]:
            down_votes[g] += 1
            source_down.setdefault(g, []).append(name)

        dataset_up[name] = set(gs["up_genes"])
        dataset_down[name] = set(gs["down_genes"])
        print(f"  {name}: {len(gs['up_genes'])} UP, {len(gs['down_genes'])} DOWN")

    # All unique genes
    all_up_candidates = set(up_votes.keys())
    all_down_candidates = set(down_votes.keys())

    # Resolve conflicts: gene voted UP in some datasets, DOWN in others
    conflicts = all_up_candidates & all_down_candidates
    resolved_up = set()
    resolved_down = set()
    conflict_details = []

    for g in sorted(conflicts):
        nu = up_votes[g]
        nd = down_votes[g]
        if nu > nd:
            resolved_up.add(g)
        elif nd > nu:
            resolved_down.add(g)
        else:
            # Tie — remove from both
            pass
        conflict_details.append(
            f"{g}: UP={nu} ({source_up.get(g, [])}), DOWN={nd} ({source_down.get(g, [])}) → "
            f"{'UP' if nu > nd else 'DOWN' if nd > nu else 'REMOVED'}"
        )

    # Non-conflicting genes
    resolved_up |= (all_up_candidates - conflicts)
    resolved_down |= (all_down_candidates - conflicts)

    # Count dataset support
    def count_support(gene, direction):
        if direction == "UP":
            return sum(1 for name, genes in dataset_up.items() if gene in genes)
        else:
            return sum(1 for name, genes in dataset_down.items() if gene in genes)

    up_with_support = [(g, count_support(g, "UP"), source_up.get(g, [])) for g in resolved_up]
    down_with_support = [(g, count_support(g, "DOWN"), source_down.get(g, [])) for g in resolved_down]

    # Sort: higher support first, then alphabetically
    up_with_support.sort(key=lambda x: (-x[1], x[0]))
    down_with_support.sort(key=lambda x: (-x[1], x[0]))

    up_genes = [g for g, _, _ in up_with_support]
    down_genes = [g for g, _, _ in down_with_support]

    # Build JSON
    gene_set = {
        "name": "MIRI_up-down",
        "species": "Multi-species union (Mus musculus + Rattus norvegicus + Sus scrofa)",
        "source_geo": "GSE160516 + GSE108940 + GSE137293",
        "comparison": "MIRI/IR vs Sham/Control (union across species)",
        "method": "union",
        "description": (
            f"Union of MIRI dysregulated genes from 3 species mapped to human L1000. "
            f"Genes are included if they appear as UP or DOWN in ANY dataset. "
            f"Conflicts (gene UP in one dataset, DOWN in another) resolved by majority vote; "
            f"ties removed. "
            f"Dataset contributions: "
            f"GSE160516 (mouse MIRI, 75+32), "
            f"GSE108940 (rat MIRI, 4+6), "
            f"GSE137293 (pig scar vs border, 2+1)."
        ),
        "n_datasets": 3,
        "n_conflicts": len(conflicts),
        "conflict_resolution": "majority_vote",
        "_conflict_details": conflict_details,
        "up_genes": up_genes,
        "down_genes": down_genes,
        "n_up_l1000": len(up_genes),
        "n_down_l1000": len(down_genes),
        "up_by_support": {
            "3_datasets": [g for g, s, _ in up_with_support if s == 3],
            "2_datasets": [g for g, s, _ in up_with_support if s == 2],
            "1_dataset": [g for g, s, _ in up_with_support if s == 1],
        },
        "down_by_support": {
            "3_datasets": [g for g, s, _ in down_with_support if s == 3],
            "2_datasets": [g for g, s, _ in down_with_support if s == 2],
            "1_dataset": [g for g, s, _ in down_with_support if s == 1],
        },
    }

    out_path = REF_DIR / "MIRI_up-down.json"
    with open(out_path, 'w') as f:
        json.dump(gene_set, f, indent=2, ensure_ascii=False)

    print(f"\n  Union UP: {len(up_genes)} (3-ds={sum(1 for _,s,_ in up_with_support if s==3)}, "
          f"2-ds={sum(1 for _,s,_ in up_with_support if s==2)}, "
          f"1-ds={sum(1 for _,s,_ in up_with_support if s==1)})")
    print(f"  Union DOWN: {len(down_genes)} (3-ds={sum(1 for _,s,_ in down_with_support if s==3)}, "
          f"2-ds={sum(1 for _,s,_ in down_with_support if s==2)}, "
          f"1-ds={sum(1 for _,s,_ in down_with_support if s==1)})")
    print(f"  Conflicts resolved: {len(conflicts)}")
    for d in conflict_details:
        print(f"    {d}")
    print(f"  Saved: {out_path}")

    return gene_set


# ═══════════════════════════════════════════════════════════════════
# Task 3: Create MIRI_rank.json (Type B reference for Zhang scoring)
# ═══════════════════════════════════════════════════════════════════

def build_rank_reference(l1000_set):
    """
    Build ranked reference from GSE160516 full DEG results.

    For Zhang Score (ssCMap):
      - Disease genes are assigned signed ranks s(g) by |disease_logFC|
      - Up-regulated genes → positive rank (+m to +1, largest |logFC| = +m)
      - Down-regulated genes → negative rank (-1 to -m, largest |logFC| = -m)

    Output format: each gene has gene symbol, logFC, direction, signed_rank.
    Genes sorted by |logFC| descending.
    """
    # Load GSE160516 DEG results
    deg = pd.read_csv(GEO_DIR / "GSE160516_deg_results.csv")

    # Also load pre-mapped L1000 files to get human symbols
    up_l1000 = pd.read_csv(GEO_DIR / "GSE160516_up_l1000.csv")
    down_l1000 = pd.read_csv(GEO_DIR / "GSE160516_down_l1000.csv")

    # Build mouse→human mapping (for ALL DEG genes, not just sig ones)
    mouse_to_human = {}
    mouse_to_logfc_up = {}
    mouse_to_logfc_down = {}

    for _, row in up_l1000.iterrows():
        mouse_to_human[row["mouse_gene"]] = row["human_gene"]
        mouse_to_logfc_up[row["mouse_gene"]] = row["logFC"]

    for _, row in down_l1000.iterrows():
        mouse_to_human[row["mouse_gene"]] = row["human_gene"]
        mouse_to_logfc_down[row["mouse_gene"]] = row["logFC"]

    # Filter DEGs: only genes that MAP to human AND are in L1000
    # (This gives us the genes we can actually score against)
    deg_mapped = []
    for _, row in deg.iterrows():
        mouse_gene = row["Gene"]
        if mouse_gene not in mouse_to_human:
            continue
        human_gene = mouse_to_human[mouse_gene]
        if human_gene not in l1000_set:
            continue
        deg_mapped.append({
            "gene": human_gene,
            "mouse_gene": mouse_gene,
            "logFC": row["logFC"],
            "P_Value": row["P.Value"],
            "adj_P_Val": row["adj.P.Val"],
        })

    print(f"\n  GSE160516 DEGs: {len(deg)} total")
    print(f"  Mapped to human + L1000: {len(deg_mapped)}")

    if len(deg_mapped) == 0:
        print("  ERROR: No GSE160516 genes map to L1000!")
        return None

    # Sort by |logFC| descending
    deg_mapped.sort(key=lambda x: abs(x["logFC"]), reverse=True)

    # Assign signed ranks within each direction group
    up_genes = [g for g in deg_mapped if g["logFC"] > 0]
    down_genes = [g for g in deg_mapped if g["logFC"] < 0]

    # Rank UP genes by |logFC| descending → rank 1 = largest |logFC|
    up_abs = np.array([abs(g["logFC"]) for g in up_genes])
    up_ranks = rankdata(-up_abs, method="average")  # 1 = largest |logFC|
    up_signed_ranks = len(up_genes) - up_ranks + 1  # +m = largest |logFC|, +1 = smallest

    down_abs = np.array([abs(g["logFC"]) for g in down_genes])
    down_ranks = rankdata(-down_abs, method="average")
    down_signed_ranks = -(len(down_genes) - down_ranks + 1)  # -m = largest |logFC|

    for i, g in enumerate(up_genes):
        g["signed_rank"] = int(up_signed_ranks[i])
        g["direction"] = "UP"

    for i, g in enumerate(down_genes):
        g["signed_rank"] = int(down_signed_ranks[i])
        g["direction"] = "DOWN"

    # Combine and sort by |logFC| descending
    all_ranked = up_genes + down_genes
    all_ranked.sort(key=lambda x: abs(x["logFC"]), reverse=True)

    rank_set = {
        "name": "MIRI_rank",
        "type": "ranked_gene_list",
        "species": "Mus musculus",
        "source_geo": "GSE160516",
        "platform": "Affymetrix Clariom S",
        "comparison": "MIRI (pooled IR) vs Sham",
        "description": (
            "Full GSE160516 MIRI DEG results filtered to human L1000 978 landmark genes. "
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

    out_path = REF_DIR / "MIRI_rank.json"
    with open(out_path, 'w') as f:
        json.dump(rank_set, f, indent=2, ensure_ascii=False)

    print(f"  UP: {len(up_genes)} (ranked +{int(max(up_signed_ranks))} to +{int(min(up_signed_ranks))})")
    print(f"  DOWN: {len(down_genes)} (ranked {int(max(down_signed_ranks))} to {int(min(down_signed_ranks))})")
    print(f"  Top 10 by |logFC|:")
    for g in all_ranked[:10]:
        print(f"    {g['gene']:10s}  dir={g['direction']:4s}  logFC={g['logFC']:+.4f}  "
              f"rank={g['signed_rank']:+d}")
    print(f"  Saved: {out_path}")

    return rank_set


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    l1000_set = load_l1000_set()
    print(f"L1000 landmark genes: {len(l1000_set)}")

    # ── Task 1: Union-based gene set ──
    print("\n═══ Task 1: MIRI_up-down.json (union) ═══")
    union_gs = build_union(l1000_set)

    # ── Task 2: Delete consensus ──
    print("\n═══ Task 2: Delete MIRI_consensus.json ═══")
    consensus_path = REF_DIR / "MIRI_consensus.json"
    if consensus_path.exists():
        consensus_path.unlink()
        print(f"  Deleted: {consensus_path}")
    else:
        print(f"  Already absent: {consensus_path}")

    # ── Task 3: Ranked reference ──
    print("\n═══ Task 3: MIRI_rank.json (Zhang reference) ═══")
    rank_set = build_rank_reference(l1000_set)

    print("\n✓ Done.")


if __name__ == "__main__":
    main()
