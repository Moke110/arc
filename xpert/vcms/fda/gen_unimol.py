"""Generate UniMol features for FDA drugs. Run in xpert conda env."""
import numpy as np
import pandas as pd
from tqdm import tqdm
from unimol_tools import UniMolRepr

# Load FDA SMILES
df = pd.read_csv("vcms/fda/fda_drugs.csv")
smi_dict = {}
for i, row in df.iterrows():
    s = str(row["canonical_smiles"]).strip()
    if s and s != "nan":
        smi_dict[i] = s
print(f"Loaded {len(smi_dict)} FDA drugs")

smiles_list = list(smi_dict.values())
print("Running UniMol inference...")
clf = UniMolRepr(data_type="molecule", remove_hs=False)
unimol_repr = clf.get_repr(smiles_list, return_atomic_reprs=True)
assert len(unimol_repr["cls_repr"]) == len(smiles_list)

atom2id = {
    "N": 2, "Se": 3, "I": 4, "As": 5, "Pt": 6, "S": 7, "Sn": 8,
    "Mg": 9, "P": 10, "Cl": 11, "Au": 12, "F": 13, "Co": 14, "Br": 15,
    "Si": 16, "O": 17, "H": 18, "B": 19, "C": 20, "Ca": 21, "Hg": 22,
    "Li": 23, "[UNK]": 24, "Na": 25, "hg_molecule": 0, "molecule": 1,
}
max_atoms = 122
n_features = 514
n_drugs = len(smi_dict)
unimol_arr = np.zeros((n_drugs, max_atoms, n_features), dtype=np.float64)

for i, k in enumerate(tqdm(smi_dict.keys(), desc="UniMol wrap")):
    mf = unimol_repr["cls_repr"][i]
    af = unimol_repr["atomic_reprs"][i]
    asym = unimol_repr["atomic_symbol"][i]

    af = np.vstack((np.array([mf, mf]), af))
    asym = ["hg_molecule", "molecule"] + list(asym)
    asym_ids = [atom2id.get(a, 24) for a in asym]

    L = len(af)
    if L < max_atoms:
        af_pad = np.pad(af, ((0, max_atoms - L), (0, 0)), constant_values=0)
        asym_ids_pad = np.pad(asym_ids, (0, max_atoms - L), constant_values=0)
        mask = [1] * L + [0] * (max_atoms - L)
    else:
        af_pad = af[:max_atoms]
        asym_ids_pad = asym_ids[:max_atoms]
        mask = [1] * max_atoms

    tmp = np.hstack((np.array(asym_ids_pad)[:, None], af_pad))
    unimol_arr[k] = np.hstack((np.array(mask)[:, None], tmp))

print(f"UniMol shape: {unimol_arr.shape}")
np.save("vcms/fda/fda_unimol.npy", unimol_arr)
print("Saved vcms/fda/fda_unimol.npy")
