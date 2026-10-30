#!/usr/bin/env python3
"""
Step 02 — Baseline transcriptome (authoritative gene mapping).

Builds the consensus cardiomyocyte baseline transcriptome from the DToxS-primary
control wells (CTRL, CTRL-DMSO, CTRL-EtOH) of the released LINCS LDS-1593
well-level package (pert_data/LINCS-LDS-1593/source/Data/Level3): computes
log2(CPM+1) and averages across all control wells into a single consensus
baseline. Gene matching uses the workflow-authoritative DToxS-primary -> L1000
mapping (03_cm_finetune_data/output/gene_maps/dtoxs_primary_to_l1000_gene_map.csv),
which resolves the 18 legacy-symbol landmarks (e.g. WASHC5 <- KIAA0196,
ELP1 <- IKBKAP, JPT2 <- HN1L) so their real control expression is not lost.
The baseline is the DToxS-primary control consensus — identical to v1's
baseline, which used the same underlying control wells. (The GSE146095/
GSE146096 raw GEO deposit was deleted 2026-08-18; the LINCS package is its
strict superset and carries the same control wells.) Idempotent (pass --force
to rebuild).

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
from utils import WORKFLOW_ROOT, XPERT_ROOT, DTOXS_LINCS_LEVEL3  # noqa: E402
import gene_mapping  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "02_cm_baseline"
OUTPUT = STEP_DIR / "output" / "cm_ctl_landmark.csv"
# The DToxS-primary gene map is produced by step 03 (shared with the union).
GENE_MAP = (WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "gene_maps"
            / "dtoxs_primary_to_l1000_gene_map.csv")

# Control state labels in the DToxS-primary experiment designs
CONTROL_STATES = {"CTRL", "CTRL-EtOH", "CTRL-DMSO"}


def load_authoritative_gene_map():
    """Return {dtoxs_source_symbol -> l1000_gene_symbol} from the DToxS-primary map."""
    if not GENE_MAP.exists():
        print(f"Building DToxS-primary gene map -> {GENE_MAP}")
        mp = gene_mapping.build_dataset_map("dtoxs_primary")
        gene_mapping.write_dataset_map("dtoxs_primary", mp)
    else:
        mp = pd.read_csv(GENE_MAP)
    valid = mp[mp["source_symbol"].notna() & (mp["source_symbol"] != "")]
    return dict(zip(valid["source_symbol"], valid["l1000_gene_symbol"]))


def build_baseline():
    dtoxs2l1000 = load_authoritative_gene_map()

    # L1000 gene info for ordering
    l1000 = pd.read_csv(XPERT_ROOT / "metadata" / "l1000_gene_info_978.csv")
    l1000_order = l1000.sort_values("gene_idx")
    landmark_set = set(l1000_order['gene_symbol'])

    # Per-experiment control wells over the released LINCS LDS-1593 package,
    # exactly like v1 (per-experiment well positions — the Cor/Spec batches
    # share control wells and each physical well is counted once per batch).
    import glob

    level3_root = DTOXS_LINCS_LEVEL3
    all_parts = []
    total_ctrl = 0
    for counts_rel in sorted(glob.glob(
            str(level3_root / "**" / "PromoCell-DGE-Read-Counts-*.tsv"),
            recursive=True)):
        design_rel = counts_rel.replace("PromoCell-DGE-Read-Counts",
                                        "PromoCell-DGE-Experiment-Design")
        if not Path(design_rel).exists():
            continue
        design = pd.read_csv(design_rel, sep='\t')
        ctrl_wells = sorted(design[design['State'].isin(CONTROL_STATES)]
                            ['Well'].tolist())
        if not ctrl_wells:
            continue
        data = pd.read_csv(counts_rel, sep='\t', index_col=0)
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
        print(f"  {Path(counts_rel).parent.name}: "
              f"{len(available)} control wells "
              f"({ctrl.shape[1]} with non-zero counts)")

    if not all_parts:
        print("ERROR: no DToxS-primary control data found under "
              f"{level3_root}")
        sys.exit(1)
    all_cpm = pd.concat(all_parts, axis=1)
    print(f"\nTotal control samples: {all_cpm.shape[1]}")

    # Map DToxS gene symbols to L1000 landmark genes (authoritative)
    dtoxs_genes = all_cpm.index.tolist()
    l1000_mapped_genes = []
    for g in dtoxs_genes:
        mapped = dtoxs2l1000.get(g, g)
        l1000_mapped_genes.append(mapped if mapped in landmark_set else None)

    mask = [g is not None for g in l1000_mapped_genes]
    l1000_cpm = all_cpm.loc[mask].copy()
    l1000_cpm.index = [l1000_mapped_genes[i] for i in range(len(l1000_mapped_genes))
                       if l1000_mapped_genes[i] is not None]

    l1000_cpm = l1000_cpm.groupby(l1000_cpm.index).mean()
    print(f"Genes mapped to L1000 landmarks: {l1000_cpm.shape[0]}/978")

    l1000_log2 = np.log2(l1000_cpm + 1)
    means = l1000_log2.mean(axis=1)
    stds = l1000_log2.std(axis=1)

    result = pd.DataFrame({
        'gene_symbol': l1000_order['gene_symbol'].values,
        'entrez_id': l1000_order['gene_id'].astype(str).values,
    })
    result['mean_log2_cpm1'] = result['gene_symbol'].map(means)
    result['std_log2_cpm1'] = result['gene_symbol'].map(stds)
    result['n_controls'] = total_ctrl

    n_missing = result['mean_log2_cpm1'].isna().sum()
    if n_missing > 0:
        print(f"Warning: {n_missing} landmark genes not found in DToxS control data "
              f"(will be set to 0)")
        result['mean_log2_cpm1'] = result['mean_log2_cpm1'].fillna(0.0)
        result['std_log2_cpm1'] = result['std_log2_cpm1'].fillna(0.0)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, index=False)

    print(f"\nBaseline summary:")
    print(f"  Mean expression: {result['mean_log2_cpm1'].mean():.3f}")
    print(f"  Range: [{result['mean_log2_cpm1'].min():.3f}, {result['mean_log2_cpm1'].max():.3f}]")
    print(f"  Genes with mean < 1: {(result['mean_log2_cpm1'] < 1).sum()}")
    print(f"  Genes exactly 0 (no DToxS control signal): {(result['mean_log2_cpm1'] == 0).sum()}")
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
