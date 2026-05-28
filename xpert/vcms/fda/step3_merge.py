"""
Step 3: Merge FDA drug embeddings into existing all_drugs_*.npy files.
Reports drug counts before and after merging.
"""
import numpy as np
import pandas as pd
from collections import OrderedDict

FDA_DIR = "vcms/fda"
EMBED_DIR = "drug_embeddings"

# Load FDA SMILES
print("=" * 60)
print("Loading FDA drug SMILES...")
fda_df = pd.read_csv(f"{FDA_DIR}/fda_drugs.csv")
fda_smiles_list = fda_df["canonical_smiles"].dropna().str.strip().tolist()
print(f"FDA drugs: {len(fda_smiles_list)}")

# Load existing SMILES dict
old_smi_dict = np.load(f"{EMBED_DIR}/all_drugs_idx2smi_8981.npy", allow_pickle=True).item()
n_existing = len(old_smi_dict)
print(f"Existing drugs: {n_existing}")

# Create SMILES set for fast lookup (canonical SMILES)
old_smiles_set = set(old_smi_dict.values())
print(f"Unique existing SMILES: {len(old_smiles_set)}")

# Find new drugs (not in existing set)
new_indices = []
existing_indices = []
for i, smi in enumerate(fda_smiles_list):
    if smi not in old_smiles_set:
        new_indices.append(i)
    else:
        existing_indices.append(i)

print(f"\nFDA drugs already in existing: {len(existing_indices)}")
print(f"FDA drugs NEW (to be added): {len(new_indices)}")

# === Merge SMILES ===
print("\n" + "=" * 60)
print("Merging SMILES...")
new_smi_dict = OrderedDict()
for k, v in old_smi_dict.items():
    new_smi_dict[k] = v

next_idx = max(new_smi_dict.keys()) + 1
fda_idx_map = {}  # Map fda_internal_idx -> new_global_idx
for fda_i in new_indices:
    fda_idx_map[fda_i] = next_idx
    new_smi_dict[next_idx] = fda_smiles_list[fda_i]
    next_idx += 1

np.save(f"{EMBED_DIR}/all_drugs_idx2smi.npy", new_smi_dict)
print(f"SMILES: {n_existing} -> {len(new_smi_dict)} (+{len(new_indices)})")
print(f"Saved {EMBED_DIR}/all_drugs_idx2smi.npy")

# === Merge HG ===
print("\n" + "=" * 60)
print("Merging HG embeddings...")
old_hg = np.load(f"{EMBED_DIR}/all_drugs_idx2HG.npy", allow_pickle=True)
print(f"Old HG shape: {old_hg.shape}")
# Load FDA HG (dict format)
fda_hg_dict = np.load(f"{FDA_DIR}/fda_HG.npy", allow_pickle=True).item()
fda_hg_new = np.stack([fda_hg_dict[i] for i in new_indices])
print(f"New HG to add: {fda_hg_new.shape}")
new_hg = np.vstack([old_hg, fda_hg_new])
np.save(f"{EMBED_DIR}/all_drugs_idx2HG.npy", new_hg.astype(np.float32))
print(f"HG: {old_hg.shape[0]} -> {new_hg.shape[0]} (+{len(new_indices)})")
print(f"Saved {EMBED_DIR}/all_drugs_idx2HG.npy")

# === Merge KPGT ===
print("\n" + "=" * 60)
print("Merging KPGT features...")
old_kpgt_dict = np.load(f"{EMBED_DIR}/all_drugs_idx2KPGT.npy", allow_pickle=True).item()
old_kpgt_n = len(old_kpgt_dict)
new_kpgt_dict = OrderedDict()
for k, v in old_kpgt_dict.items():
    new_kpgt_dict[k] = v
fda_kpgt_dict = np.load(f"{FDA_DIR}/fda_KPGT.npy", allow_pickle=True).item()
nk = max(new_kpgt_dict.keys()) + 1
for fda_i in new_indices:
    new_kpgt_dict[nk] = fda_kpgt_dict[fda_i]
    nk += 1
np.save(f"{EMBED_DIR}/all_drugs_idx2KPGT.npy", new_kpgt_dict)
print(f"KPGT: {old_kpgt_n} -> {len(new_kpgt_dict)} (+{len(new_indices)})")
print(f"Saved {EMBED_DIR}/all_drugs_idx2KPGT.npy")

# === Merge Morgan ===
print("\n" + "=" * 60)
print("Merging Morgan fingerprints...")
old_morgan_dict = np.load(f"{EMBED_DIR}/all_drugs_idx2morgan.npy", allow_pickle=True).item()
old_morgan_n = len(old_morgan_dict)
new_morgan_dict = OrderedDict()
for k, v in old_morgan_dict.items():
    new_morgan_dict[k] = v
fda_morgan_dict = np.load(f"{FDA_DIR}/fda_morgan.npy", allow_pickle=True).item()
nk = max(new_morgan_dict.keys()) + 1
for fda_i in new_indices:
    new_morgan_dict[nk] = fda_morgan_dict[fda_i]
    nk += 1
np.save(f"{EMBED_DIR}/all_drugs_idx2morgan.npy", new_morgan_dict)
print(f"Morgan: {old_morgan_n} -> {len(new_morgan_dict)} (+{len(new_indices)})")
print(f"Saved {EMBED_DIR}/all_drugs_idx2morgan.npy")

# === Merge UniMol ===
print("\n" + "=" * 60)
print("Merging UniMol features...")
old_unimol = np.load(f"{EMBED_DIR}/all_drugs_unimol_arr.npy", allow_pickle=True)
print(f"Old UniMol shape: {old_unimol.shape}")
fda_unimol = np.load(f"{FDA_DIR}/fda_unimol.npy", allow_pickle=True)
fda_unimol_new = fda_unimol[list(new_indices)]  # select only new drugs
print(f"New UniMol to add: {fda_unimol_new.shape}")
new_unimol = np.concatenate([old_unimol, fda_unimol_new], axis=0)
np.save(f"{EMBED_DIR}/all_drugs_unimol_arr.npy", new_unimol.astype(np.float64))
print(f"UniMol: {old_unimol.shape[0]} -> {new_unimol.shape[0]} (+{len(new_indices)})")
print(f"Saved {EMBED_DIR}/all_drugs_unimol_arr.npy")

# Summary
print("\n" + "=" * 60)
print("MERGE SUMMARY")
print("=" * 60)
print(f"Existing drugs:        {n_existing}")
print(f"FDA drugs total:       {len(fda_smiles_list)}")
print(f"Already in dataset:    {len(existing_indices)}")
print(f"New drugs added:       {len(new_indices)}")
print(f"Total after merge:     {n_existing + len(new_indices)}")
print(f"\nAll files updated in {EMBED_DIR}/")
