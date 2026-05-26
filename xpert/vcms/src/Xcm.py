"""
Build unified control cardiomyocyte transcriptome X_cm by aggregating all
matched-pair control profiles from dtoxs.h5ad obsm['X_ctl'].
"""
import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import anndata


def main():
    project_root = Path(__file__).resolve().parent.parent.parent
    h5ad_path = project_root / "finetune" / "data" / "dtoxs.h5ad"
    output_path = project_root / "vcms" / "Xcm.csv"

    print(f"Loading {h5ad_path}...")
    adata = anndata.read_h5ad(h5ad_path)

    X_ctl = adata.obsm["X_ctl"]
    print(f"X_ctl shape: {X_ctl.shape}, type: {type(X_ctl).__name__}")

    X_ctl_dense = X_ctl.toarray() if hasattr(X_ctl, "toarray") else X_ctl
    X_cm = X_ctl_dense.mean(axis=0)  # (978,)

    print(f"X_cm shape: {X_cm.shape}, range: [{X_cm.min():.4f}, {X_cm.max():.4f}]")

    gene_symbols = list(adata.var["gene_symbol"])
    assert len(gene_symbols) == 978, f"Expected 978 genes, got {len(gene_symbols)}"

    df = pd.DataFrame([X_cm], columns=gene_symbols)
    df.to_csv(output_path, index=False)
    print(f"Saved X_cm to {output_path}")
    print(f"  Columns: {len(df.columns)}, Values: non-null, range [{X_cm.min():.4f}, {X_cm.max():.4f}]")


if __name__ == "__main__":
    main()
