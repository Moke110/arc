"""
Post-process fda_cm_perturb.csv to produce:
  - fda_cm_deg.csv: pure deg (logFC) = y_pred - ctl_mean
  - fda_cm_degrank.csv: gene ranks per drug, highest deg first
"""
import numpy as np
import pandas as pd
from scipy import sparse
import h5py

FDA_DIR = "vcms/fda"
PERTURB_CSV = f"{FDA_DIR}/fda_cm_perturb.csv"

# Load ctl_mean from DTOXS data
print("Loading DTOXS control mean...")
with h5py.File("finetune/data/dtoxs.h5ad", "r") as f:
    g = f["obsm/X_ctl"]
    X_ctl = sparse.csr_matrix(
        (g["data"][:], g["indices"][:], g["indptr"][:]),
        shape=(int(g.attrs.get("shape", [521, 978])[0]), int(g.attrs.get("shape", [521, 978])[1]))
    ).toarray()
    gene_symbols = f["var/gene_symbol"][:]
ctl_mean = X_ctl.mean(axis=0).astype(np.float32)
gene_names = [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in gene_symbols]
print(f"ctl_mean shape: {ctl_mean.shape}, genes: {len(gene_names)}")

# Load predictions
print(f"Loading {PERTURB_CSV}...")
df_pred = pd.read_csv(PERTURB_CSV)
smiles = df_pred["canonical_smiles"].values
expr_cols = [c for c in df_pred.columns if c != "canonical_smiles"]
y_pred = df_pred[expr_cols].values.astype(np.float32)
print(f"y_pred shape: {y_pred.shape}")

# Verify gene order matches
assert len(expr_cols) == len(gene_names), f"Column mismatch: {len(expr_cols)} vs {len(gene_names)}"
assert list(expr_cols) == list(gene_names), "Gene name/order mismatch!"

# Compute deg
deg = y_pred - ctl_mean[np.newaxis, :]
print(f"deg shape: {deg.shape}")
print(f"deg range: [{deg.min():.4f}, {deg.max():.4f}], mean: {deg.mean():.4f}")

# Save deg
df_deg = pd.DataFrame(deg, columns=gene_names)
df_deg.insert(0, "canonical_smiles", smiles)
df_deg.to_csv(f"{FDA_DIR}/fda_cm_deg.csv", index=False)
print(f"Saved {FDA_DIR}/fda_cm_deg.csv")

# Compute rank per drug (highest deg = rank 1, ties averaged)
print("Computing gene ranks (1 = most up-regulated)...")
ranks = np.zeros_like(deg, dtype=np.float32)
for i in range(deg.shape[0]):
    # argsort descending: highest deg first
    order = np.argsort(-deg[i])
    ranks[i, order] = np.arange(1, len(order) + 1, dtype=np.float32)

df_rank = pd.DataFrame(ranks, columns=gene_names)
df_rank = df_rank.astype(int)
df_rank.insert(0, "canonical_smiles", smiles)
df_rank.to_csv(f"{FDA_DIR}/fda_cm_degrank.csv", index=False)
print(f"Saved {FDA_DIR}/fda_cm_degrank.csv")

# Quick summary
print("\nDone! Summary for first 5 drugs:")
for i in range(min(5, len(smiles))):
    top5_idx = np.argsort(-deg[i])[:5]
    top5_genes = [f"{gene_names[j]}({deg[i,j]:.3f})" for j in top5_idx]
    print(f"  {smiles[i][:60]}...")
    print(f"    Top5 up: {', '.join(top5_genes)}")
