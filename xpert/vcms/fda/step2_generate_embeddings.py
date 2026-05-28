"""
Step 2: Generate all 4 embedding types for FDA drugs.
Run from xpert env for morgan, unimol, hg.
KPGT requires the KPGT conda env (separate script).
"""
import numpy as np
import pandas as pd
import os
import sys
import csv
import pickle
from tqdm import tqdm

FDA_DIR = "vcms/fda"
EMBED_DIR = "drug_embeddings"

def load_fda_smiles():
    """Load FDA drug SMILES, return dict {index: canonical_smiles}"""
    df = pd.read_csv(f"{FDA_DIR}/fda_drugs.csv")
    smi_dict = {}
    for i, row in df.iterrows():
        smiles = str(row['canonical_smiles']).strip()
        if smiles and smiles != 'nan':
            smi_dict[i] = smiles
    print(f"Loaded {len(smi_dict)} FDA drugs with valid SMILES")
    return smi_dict

def generate_morgan(smi_dict, fp_size=1024, radius=2):
    """Generate Morgan fingerprints for FDA drugs."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    morgan_dict = {}
    for idx, smiles in tqdm(smi_dict.items(), desc="Morgan fingerprints"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=fp_size)
            arr = np.zeros(fp_size, dtype=np.float32)
            for bit in fp.GetOnBits():
                arr[bit] = 1.0
            morgan_dict[idx] = arr
        else:
            morgan_dict[idx] = np.zeros(fp_size, dtype=np.float32)
            print(f"WARNING: Invalid SMILES at idx {idx}: {smiles}")

    return morgan_dict

def generate_unimol(smi_dict, max_atoms=122):
    """Generate UniMol features for FDA drugs."""
    from unimol_tools import UniMolRepr

    smiles_list = list(smi_dict.values())
    print("Running UniMol inference...")
    clf = UniMolRepr(data_type='molecule', remove_hs=False)
    unimol_repr = clf.get_repr(smiles_list, return_atomic_reprs=True)
    assert len(unimol_repr['cls_repr']) == len(smiles_list)

    # Reuse the atom2id mapping from xpert
    atom2id = {'N': 2, 'Se': 3, 'I': 4, 'As': 5, 'Pt': 6, 'S': 7, 'Sn': 8,
               'Mg': 9, 'P': 10, 'Cl': 11, 'Au': 12, 'F': 13, 'Co': 14, 'Br': 15,
               'Si': 16, 'O': 17, 'H': 18, 'B': 19, 'C': 20, 'Ca': 21, 'Hg': 22,
               'Li': 23, "[UNK]": 24, "Na": 25, 'hg_molecule': 0, 'molecule': 1}
    n_features = 514  # 512 unimol + 2 extra (atom_symbol_id, input_mask)

    unimol_feat = {}
    for i, k in enumerate(tqdm(smi_dict.keys(), desc="UniMol wrapping")):
        molecule_feat = unimol_repr['cls_repr'][i]
        atom_feat = unimol_repr['atomic_reprs'][i]
        atom_symbols = unimol_repr['atomic_symbol'][i]

        atom_feat = np.vstack((np.array([molecule_feat, molecule_feat]), atom_feat))
        atom_symbols = [*['hg_molecule', 'molecule'], *atom_symbols]
        atom_symbols_ids = [atom2id.get(atom, 24) for atom in atom_symbols]

        l = len(atom_feat)
        if l < max_atoms:
            atom_feat_pad = np.pad(atom_feat, ((0, max_atoms - l), (0, 0)), 'constant', constant_values=0)
            atom_symbols_ids_pad = np.pad(atom_symbols_ids, (0, max_atoms - l), 'constant', constant_values=0)
            input_mask = ([1] * l) + ([0] * (max_atoms - l))
        else:
            atom_feat_pad = atom_feat[:max_atoms]
            atom_symbols_ids_pad = atom_symbols_ids[:max_atoms]
            input_mask = [1] * max_atoms

        tmp = np.hstack((np.array(atom_symbols_ids_pad)[:, None], atom_feat_pad))
        unimol_feat[k] = np.hstack((np.array(input_mask)[:, None], tmp))

    # Convert to 3D array (n_drugs, max_atoms, n_features)
    n_drugs = len(unimol_feat)
    unimol_arr = np.zeros((n_drugs, max_atoms, n_features), dtype=np.float64)
    for k, v in unimol_feat.items():
        unimol_arr[k] = v

    return unimol_arr

def generate_hg(kpgt_features_fda, model_path, hg_data_path):
    """
    Generate HG embeddings for FDA drugs using:
    1. PCA projection: KPGT (2304) → 512 (trained on shared drugs)
    2. HG model inference (isolated nodes)
    """
    import torch
    import sys
    sys.path.insert(0, '/home/chang/arc/hg')
    from src.model_hg import HeteroGNN
    from torch_geometric.data import HeteroData

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # Load original data to learn mapping
    kpgt_orig = np.load(f"{EMBED_DIR}/all_drugs_idx2KPGT.npy", allow_pickle=True).item()
    kpgt_orig_arr = np.array([kpgt_orig[i] for i in range(len(kpgt_orig))])

    drug_node_feat = np.load(f"{hg_data_path}/all_drug_node_feat_8981.npy")

    # Learn linear projection: KPGT → HG node feat (ridge regression via SVD)
    print("Learning KPGT → HG-node-feat projection...")
    # X = KPGT, Y = HG_node_feat
    # W = (X^T X + λI)^(-1) X^T Y
    X = kpgt_orig_arr
    Y = drug_node_feat
    lam = 1.0
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    # Ridge: W = V * diag(S/(S^2+lam)) * U^T * Y
    S_reg = S / (S**2 + lam)
    W = Vt.T @ np.diag(S_reg) @ U.T @ Y
    print(f"Projection matrix: {W.shape}")

    # Apply projection to FDA KPGT features
    fda_kpgt_arr = np.array([kpgt_features_fda[i] for i in range(len(kpgt_features_fda))])
    fda_hg_input = fda_kpgt_arr @ W

    # Load HG model
    print("Loading HG model...")
    model = HeteroGNN(hidden_dim=256, out_dim=256, num_layers=3,
                      include_edges=['PPI', 'DTI', 'DDS']).to(device)

    ckpt = torch.load(model_path, map_location=device, weights_only=False)
    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()

    # Create HeteroData with new drugs as isolated nodes + existing gene nodes
    gene_feat = torch.tensor(np.load(f"{hg_data_path}/all_gene_node_feat_19392.npy"),
                             dtype=torch.float, device=device)

    # Load edge indices
    ppi = np.load(f"{hg_data_path}/PPI_edge_index_all_901260_pairs_weighted.npy",
                  allow_pickle=True).item()
    dti = np.load(f"{hg_data_path}/DTI_edge_index_all_12890_pairs.npy",
                  allow_pickle=True).item()
    dds = np.load(f"{hg_data_path}/DDS_edge_index_all_287834_pairs_weighted.npy",
                  allow_pickle=True).item()

    fda_hg_input_tensor = torch.tensor(fda_hg_input, dtype=torch.float, device=device)

    data = HeteroData()
    data['gene'].x = gene_feat
    data['drug'].x = fda_hg_input_tensor
    data['gene', 'PPI', 'gene'].edge_index = torch.tensor(ppi['edge_index'], dtype=torch.long, device=device)
    data['drug', 'DTI', 'gene'].edge_index = torch.tensor(dti['edge_index'], dtype=torch.long, device=device)
    data['gene', 'DTI', 'drug'].edge_index = torch.tensor(np.stack([dti['edge_index'][1], dti['edge_index'][0]]), dtype=torch.long, device=device)
    data['drug', 'DDS', 'drug'].edge_index = torch.tensor(dds['edge_index'], dtype=torch.long, device=device)

    print(f"Running HG inference for {fda_hg_input_tensor.shape[0]} drugs...")
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        hg_embeddings = x_dict['drug'].cpu().numpy()

    # Convert to dict format matching existing
    hg_dict = {i: hg_embeddings[i] for i in range(len(hg_embeddings))}
    return hg_dict

def main():
    os.makedirs(FDA_DIR, exist_ok=True)

    smi_dict = load_fda_smiles()

    # === Morgan fingerprints ===
    print("\n" + "="*60)
    print("Generating Morgan fingerprints...")
    morgan_dict = generate_morgan(smi_dict)
    np.save(f"{FDA_DIR}/fda_morgan.npy", morgan_dict)
    print(f"Saved fda_morgan.npy: {len(morgan_dict)} drugs, dim={morgan_dict[0].shape}")

    # === UniMol features ===
    print("\n" + "="*60)
    print("Generating UniMol features...")
    unimol_arr = generate_unimol(smi_dict)
    np.save(f"{FDA_DIR}/fda_unimol.npy", unimol_arr)
    print(f"Saved fda_unimol.npy: shape={unimol_arr.shape}")

    print("\n" + "="*60)
    print("Done! Now run step2_kpgt.sh in KPGT conda env for KPGT features.")
    print("Then run step2_hg.sh in xpert conda env for HG embeddings.")

if __name__ == "__main__":
    main()
