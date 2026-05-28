"""Extract KPGT features for FDA drugs. Run in KPGT conda env."""
import sys; sys.path.insert(0, "/home/chang/arc/kpgt")
import torch, numpy as np, os
from torch.utils.data import DataLoader
from src.utils import set_random_seed
from src.data.featurizer import Vocab, N_ATOM_TYPES, N_BOND_TYPES
from src.data.finetune_dataset import MoleculeDataset
from src.data.collator import Collator_tune
from src.model.light import LiGhTPredictor as LiGhT
from src.model_config import config_dict

set_random_seed(22, 1)
config = config_dict["base"]
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")

vocab = Vocab(N_ATOM_TYPES, N_BOND_TYPES)
collator = Collator_tune(config["path_length"])
dp = "/home/chang/arc/kpgt/feature_extract"
mol_dataset = MoleculeDataset(root_path=dp, dataset="FDA_drugs", dataset_type=None)
loader = DataLoader(mol_dataset, batch_size=32, shuffle=False, num_workers=4, drop_last=False, collate_fn=collator)

model = LiGhT(
    d_node_feats=config["d_node_feats"], d_edge_feats=config["d_edge_feats"],
    d_g_feats=config["d_g_feats"], d_hpath_ratio=config["d_hpath_ratio"],
    n_mol_layers=config["n_mol_layers"], path_length=config["path_length"],
    n_heads=config["n_heads"], n_ffn_dense_layers=config["n_ffn_dense_layers"],
    input_drop=0, attn_drop=0, feat_drop=0, n_node_types=vocab.vocab_size
).to(device)

mp_full = "/home/chang/arc/kpgt/models/pretrained/base/base.pth"
state = torch.load(mp_full, map_location=device)
model.load_state_dict({k.replace("module.", ""): v for k, v in state.items()})
print("Model loaded!")

fps_list = []
for batch_idx, batched_data in enumerate(loader):
    (_, g, ecfp, md, labels) = batched_data
    ecfp, md, g = ecfp.to(device), md.to(device), g.to(device)
    fps = model.generate_fps(g, ecfp, md)
    fps_list.extend(fps.detach().cpu().numpy().tolist())
    if batch_idx % 20 == 0:
        print(f"  Batch {batch_idx}/{len(loader)}")

out = np.array(fps_list)
print(f"KPGT features shape: {out.shape}")
np.savez_compressed(f"{dp}/FDA_drugs/kpgt_base.npz", fps=out)

fda_kpgt_dict = {i: out[i] for i in range(len(out))}
np.save("/home/chang/arc/xpert/vcms/fda/fda_KPGT.npy", fda_kpgt_dict)
print(f"Saved fda_KPGT.npy: {len(fda_kpgt_dict)} drugs x {out.shape[1]} dim")
