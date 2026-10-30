#!/usr/bin/env python3
"""
Step 02 — Baseline transcriptome (authoritative gene mapping).

Builds the consensus cardiomyocyte baseline transcriptome from DTOXS control
wells (CTRL, CTRL-DMSO, CTRL-EtOH): computes log2(CPM+1) and averages across
all control wells into a single consensus baseline. Gene matching uses the
workflow-authoritative DToxS->L1000 mapping (02_cm_baseline/output/
dtoxs_to_l1000_gene_map.csv), which resolves the 18 legacy-symbol landmarks
(e.g. WASHC5 <- KIAA0196, ELP1 <- IKBKAP, JPT2 <- HN1L) by Entrez ID so their
real control expression is not lost. Idempotent (pass --force to rebuild).

Output: 02_cm_baseline/output/cm_ctl_landmark.csv
  Columns: gene_symbol, entrez_id, mean_log2_cpm1, std_log2_cpm1, n_controls
  Rows: 978 landmark genes in L1000 gene_idx order

Usage:
    python 02_cm_baseline/cm_baseline.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, WORKSPACE_ROOT, XPERT_ROOT  # noqa: E402
import gene_mapping  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "02_cm_baseline"
OUTPUT = STEP_DIR / "output" / "cm_ctl_landmark.csv"
GENE_MAP = STEP_DIR / "output" / "dtoxs_to_l1000_gene_map.csv"
DTOXS_ROOT = WORKSPACE_ROOT / "pert_data" / "LINCS-LDS-1593" / "local" / "v1-frozen"

# Control state labels in DTOXS experiment designs
CONTROL_STATES = {"CTRL", "CTRL-EtOH", "CTRL-DMSO"}

# Experiment batches: (name, design_tsv, counts_tsv)
EXPERIMENTS = [
    ("Exp1", "Sequencing Experiment Number 1/PromoCell-DGE-Experiment-Design-20150409.tsv",
             "Sequencing Experiment Number 1/PromoCell-DGE-Read-Counts-20150409.tsv"),
    ("Exp2", "Sequencing Experiment Number 2/PromoCell-DGE-Experiment-Design-20150503.tsv",
             "Sequencing Experiment Number 2/PromoCell-DGE-Read-Counts-20150503.tsv"),
    ("Exp4", "Sequencing Experiment Number 4/PromoCell-DGE-Experiment-Design-20151120.tsv",
             "Sequencing Experiment Number 4/PromoCell-DGE-Read-Counts-20151120.tsv"),
    ("Exp5_Cor", "Sequencing Experiment Number 5/Cor/PromoCell-DGE-Experiment-Design-20161108-Cor.tsv",
                 "Sequencing Experiment Number 5/Cor/PromoCell-DGE-Read-Counts-20161108-Cor.tsv"),
    ("Exp5_SF", "Sequencing Experiment Number 5/Spec-Fixed/PromoCell-DGE-Experiment-Design-20161108-Spec-Fixed.tsv",
                "Sequencing Experiment Number 5/Spec-Fixed/PromoCell-DGE-Read-Counts-20161108-Spec-Fixed.tsv"),
    ("Exp5_SV1", "Sequencing Experiment Number 5/Spec-Var-1/PromoCell-DGE-Experiment-Design-20161108-Spec-Var-1.tsv",
                 "Sequencing Experiment Number 5/Spec-Var-1/PromoCell-DGE-Read-Counts-20161108-Spec-Var-1.tsv"),
    ("Exp5_SV2", "Sequencing Experiment Number 5/Spec-Var-2/PromoCell-DGE-Experiment-Design-20161108-Spec-Var-2.tsv",
                 "Sequencing Experiment Number 5/Spec-Var-2/PromoCell-DGE-Read-Counts-20161108-Spec-Var-2.tsv"),
    ("Exp5_SV3", "Sequencing Experiment Number 5/Spec-Var-3/PromoCell-DGE-Experiment-Design-20161108-Spec-Var-3.tsv",
                 "Sequencing Experiment Number 5/Spec-Var-3/PromoCell-DGE-Read-Counts-20161108-Spec-Var-3.tsv"),
]


def load_authoritative_gene_map():
    """Return {dtoxs_source_symbol -> l1000_gene_symbol} from the workflow mapping."""
    if not GENE_MAP.exists():
        print(f"Building authoritative gene map -> {GENE_MAP}")
        raw, _ = gene_mapping.load_dtoxs_raw_counts()
        sym2entrez = gene_mapping.build_symbol_to_entrez(set(raw.index))
        mp = gene_mapping.build_mapping(raw, [], sym2entrez)
        GENE_MAP.parent.mkdir(parents=True, exist_ok=True)
        mp.to_csv(GENE_MAP, index=False)
    else:
        mp = pd.read_csv(GENE_MAP)
    valid = mp[mp["dtoxs_source_symbol"].notna()
               & (mp["dtoxs_source_symbol"] != "")]
    return dict(zip(valid["dtoxs_source_symbol"], valid["l1000_gene_symbol"]))


def build_baseline():
    level3_root = DTOXS_ROOT / "Data" / "Level3"
    dtoxs2l1000 = load_authoritative_gene_map()

    # Load L1000 gene info for ordering
    l1000 = pd.read_csv(XPERT_ROOT / "metadata" / "l1000_gene_info_978.csv")
    l1000_order = l1000.sort_values("gene_idx")
    landmark_set = set(l1000_order['gene_symbol'])

    # Collect control CPM data
    all_parts = []
    total_ctrl = 0

    for name, design_rel, counts_rel in EXPERIMENTS:
        design_path = level3_root / design_rel
        counts_path = level3_root / counts_rel

        if not design_path.exists() or not counts_path.exists():
            print(f"  Skipping {name}: files not found")
            continue

        design = pd.read_csv(design_path, sep='\t')
        ctrl_wells = sorted(design[design['State'].isin(CONTROL_STATES)]['Well'].tolist())
        if not ctrl_wells:
            continue

        data = pd.read_csv(counts_path, sep='\t', index_col=0)
        available = [w for w in ctrl_wells if w in data.columns]
        if not available:
            continue

        ctrl = data[available].astype(float)
        totals = ctrl.sum(axis=0)
        valid = totals > 0
        ctrl = ctrl.loc[:, valid]
        totals = totals[valid]
        cpm = ctrl.div(totals, axis=1) * 1e6
        all_parts.append(cpm)
        total_ctrl += cpm.shape[1]
        print(f"  {name}: {len(available)} control wells ({ctrl.shape[1]} with non-zero counts)")

    if not all_parts:
        print("ERROR: No control data found!")
        sys.exit(1)

    # Concatenate all control CPMs (align by gene symbol)
    all_cpm = pd.concat(all_parts, axis=1)
    print(f"\nTotal control samples: {all_cpm.shape[1]}")

    # Map DTOXS gene symbols to L1000 landmark genes (authoritative, Entrez-verified)
    dtoxs_genes = all_cpm.index.tolist()
    l1000_mapped_genes = []
    for g in dtoxs_genes:
        mapped = dtoxs2l1000.get(g, g)
        l1000_mapped_genes.append(mapped if mapped in landmark_set else None)

    mask = [g is not None for g in l1000_mapped_genes]
    l1000_cpm = all_cpm.loc[mask].copy()
    l1000_cpm.index = [l1000_mapped_genes[i] for i in range(len(l1000_mapped_genes))
                       if l1000_mapped_genes[i] is not None]

    # Average across wells for genes that appear multiple times
    l1000_cpm = l1000_cpm.groupby(l1000_cpm.index).mean()
    print(f"Genes mapped to L1000 landmarks: {l1000_cpm.shape[0]}/978")

    # log2(CPM+1) transformation
    l1000_log2 = np.log2(l1000_cpm + 1)
    means = l1000_log2.mean(axis=1)
    stds = l1000_log2.std(axis=1)

    # Build result in L1000 gene_idx order
    result = pd.DataFrame({
        'gene_symbol': l1000_order['gene_symbol'].values,
        'entrez_id': l1000_order['gene_id'].astype(str).values,
    })
    result['mean_log2_cpm1'] = result['gene_symbol'].map(means)
    result['std_log2_cpm1'] = result['gene_symbol'].map(stds)
    result['n_controls'] = total_ctrl

    n_missing = result['mean_log2_cpm1'].isna().sum()
    if n_missing > 0:
        print(f"Warning: {n_missing} landmark genes not found in DTOXS control data "
              f"(will be set to 0)")
        result['mean_log2_cpm1'] = result['mean_log2_cpm1'].fillna(0.0)
        result['std_log2_cpm1'] = result['std_log2_cpm1'].fillna(0.0)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, index=False)

    print(f"\nBaseline summary:")
    print(f"  Mean expression: {result['mean_log2_cpm1'].mean():.3f}")
    print(f"  Range: [{result['mean_log2_cpm1'].min():.3f}, {result['mean_log2_cpm1'].max():.3f}]")
    print(f"  Genes with mean < 1: {(result['mean_log2_cpm1'] < 1).sum()}")
    print(f"  Genes exactly 0 (no DTOXS control signal): {(result['mean_log2_cpm1'] == 0).sum()}")
    print(f"  Saved to: {OUTPUT}")


def main():
    parser = argparse.ArgumentParser(description="Build the workflow-owned baseline transcriptome")
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if the cached baseline exists")
    args = parser.parse_args()

    if OUTPUT.exists() and not args.force:
        print(f"Cached baseline present: {OUTPUT} — skipping")
        return

    build_baseline()
    print(f"Step 02 complete — baseline cached at {OUTPUT}")


if __name__ == "__main__":
    main()
