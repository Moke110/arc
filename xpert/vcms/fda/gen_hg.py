"""
Generate HG embeddings for FDA drugs.
Approach: Learn linear projection from KPGT(2304) to HG_input(512)
using the 8981 shared drugs, then run HG model inference.
Run in xpert conda env.
"""
import numpy as np
import torch
import sys; sys.path.insert(0, "/home/chang/arc/hg")
from src.model_hg import HeteroGNN
from torch_geometric.data import HeteroData
from tqdm import tqdm

# Paths
FDA_DIR = "vcms/fda"
EMBED_DIR = "drug_embeddings"
HG_DATA = "/home/chang/arc/hg/data"
HG_MODEL = "/home/chang/arc/hg/model/20260518_102121/global_fold_early_stop.pth"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")

# Load existing KPGT and HG features for 8981 drugs
print("Loading shared drug features...")
kpgt_orig = np.load(f"{EMBED_DIR}/all_drugs_idx2KPGT.npy", allow_pickle=True).item()
kpgt_orig_arr = np.array([kpgt_orig[i] for i in range(len(kpgt_orig))], dtype=np.float32)
print(f"  KPGT orig: {kpgt_orig_arr.shape}")

drug_node_feat = np.load(f"{HG_DATA}/all_drug_node_feat_8981.npy").astype(np.float32)
print(f"  HG node feat: {drug_node_feat.shape}")

hg_orig = np.load(f"{EMBED_DIR}/all_drugs_idx2HG.npy", allow_pickle=True)
print(f"  HG embedding: {hg_orig.shape}")

# Load FDA KPGT features
print("Loading FDA KPGT features...")
kpgt_fda = np.load(f"{FDA_DIR}/fda_KPGT.npy", allow_pickle=True).item()
kpgt_fda_arr = np.array([kpgt_fda[i] for i in range(len(kpgt_fda))], dtype=np.float32)
print(f"  FDA KPGT: {kpgt_fda_arr.shape}")

# Learn linear projection: KPGT (2304) -> HG node feat (512)
print("Learning KPGT -> HG-node-feat projection...")
# Ridge regression via SVD
X = kpgt_orig_arr  # (8981, 2304)
Y = drug_node_feat  # (8981, 512)
lam = 10.0
U, S, Vt = np.linalg.svd(X, full_matrices=False)
S_reg = S / (S ** 2 + lam)
W_proj = Vt.T @ np.diag(S_reg) @ U.T @ Y  # (2304, 512)
print(f"  Projection matrix: {W_proj.shape}")

# Apply projection to FDA KPGT
fda_hg_input = kpgt_fda_arr @ W_proj  # (3127, 512)
print(f"  FDA HG input: {fda_hg_input.shape}")

# Load HG model
print("Loading HG model...")
model = HeteroGNN(hidden_dim=256, out_dim=256, num_layers=3,
                  include_edges=["PPI", "DTI", "DDS"]).to(device)
ckpt = torch.load(HG_MODEL, map_location=device, weights_only=False)
model.load_state_dict(ckpt["model_state_dict"])
model.eval()

# Load gene features and edge indices from the existing HG graph
print("Building HeteroData for inference...")
gene_feat = torch.tensor(
    np.load(f"{HG_DATA}/all_gene_node_feat_19392.npy").astype(np.float32),
    dtype=torch.float, device=device
)

# Load all edge indices
ppi = np.load(f"{HG_DATA}/PPI_edge_index_all_901260_pairs_weighted.npy", allow_pickle=True).item()
dti = np.load(f"{HG_DATA}/DTI_edge_index_all_12890_pairs.npy", allow_pickle=True).item()
dds = np.load(f"{HG_DATA}/DDS_edge_index_all_287834_pairs_weighted.npy", allow_pickle=True).item()

# We need to add new drugs to the HG graph
# The new drugs share gene nodes with the original graph
# but have no DTI or DDS edges to existing drugs
# However, DTI edges connect drugs to genes; we don't have DTI for new drugs
# We add the new drugs as isolated drug nodes

# Combine original drug features + new FDA drug features
drug_feat_combined = np.vstack([drug_node_feat, fda_hg_input])  # (8981+3127, 512)
drug_feat_tensor = torch.tensor(drug_feat_combined, dtype=torch.float, device=device)
fda_hg_input_tensor = torch.tensor(fda_hg_input, dtype=torch.float, device=device)

# Create HeteroData with all drugs (old + new) and all genes
# Old edges remain valid since they reference drug indices 0..8980
data = HeteroData()
data["gene"].x = gene_feat
data["drug"].x = drug_feat_tensor

# Add edges (drug indices 0..8980 are original drugs, 8981..12108 are FDA)
data["gene", "PPI", "gene"].edge_index = torch.tensor(ppi["edge_index"], dtype=torch.long, device=device)
data["drug", "DTI", "gene"].edge_index = torch.tensor(dti["edge_index"], dtype=torch.long, device=device)
data["gene", "DTI", "drug"].edge_index = torch.tensor(
    np.stack([dti["edge_index"][1], dti["edge_index"][0]]), dtype=torch.long, device=device
)
data["drug", "DDS", "drug"].edge_index = torch.tensor(dds["edge_index"], dtype=torch.long, device=device)

# Run inference
print(f"Running HG inference on {drug_feat_combined.shape[0]} drugs ({len(fda_hg_input)} new)...")
with torch.no_grad():
    x_dict = model(data.x_dict, data.edge_index_dict)
    all_hg_embeddings = x_dict["drug"].cpu().numpy()  # (12108, 256)

# Extract FDA drug embeddings (last n_fda entries)
fda_hg_embeddings = all_hg_embeddings[-len(fda_hg_input):]  # (3127, 256)
print(f"FDA HG embeddings shape: {fda_hg_embeddings.shape}")

# Save as dict format matching existing (indexed 0..3126)
hg_dict = {i: fda_hg_embeddings[i] for i in range(len(fda_hg_embeddings))}
np.save(f"{FDA_DIR}/fda_HG.npy", hg_dict)
print(f"Saved {FDA_DIR}/fda_HG.npy")

# Also save the combined HG embeddings (old + new) for merging
combined_hg = all_hg_embeddings  # (12108, 256)
np.save(f"{FDA_DIR}/fda_HG_all.npy", combined_hg)
print(f"Saved {FDA_DIR}/fda_HG_all.npy (combined, shape={combined_hg.shape})")

print("\nDone! HG embeddings generated successfully.")
