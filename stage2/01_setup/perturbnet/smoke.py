#!/usr/bin/env python
"""PerturbNet setup smoke: load the released cINN + ChemicalVAE and run a forward.

The released checkpoints are the third-party reproduction assets in
`pert_model/perturbnet/assets/`. Full end-to-end inference needs the original
Python 3.7 / TensorFlow 1.15 cell-VAE; this smoke proves the checkpoint and the
in-repo flow architecture are compatible on CPU.
"""
import json
import sys
import types
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
MODEL = ROOT / "pert_model" / "perturbnet"
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

# The flow modules import torchvision at module load even though this path does
# not use it; the pinned CPU env intentionally has no torchvision.
tv = types.ModuleType("torchvision")
tv.models = types.SimpleNamespace()
tv.transforms = types.SimpleNamespace()
sys.modules.setdefault("torchvision", tv)

sys.path.insert(0, str(MODEL))
from perturbnet.net2net.modules.flow.flatflow import ConditionalFlatCouplingFlow  # noqa: E402
from perturbnet.chemicalvae.chemicalVAE import ChemicalVAE  # noqa: E402


def load_cinn(path: Path, cond_dim: int):
    ckpt = torch.load(path, map_location="cpu", weights_only=False)
    flow = ConditionalFlatCouplingFlow(
        in_channels=10,
        conditioning_dim=cond_dim,
        embedding_dim=10,
        hidden_dim=1024,
        hidden_depth=2,
        n_flows=20,
        conditioning_option="none",
        activation="none",
        conditioner_use_bn=True,
    )
    flow.load_state_dict({k[5:]: v for k, v in ckpt.items() if k.startswith("flow.")})
    chem = ChemicalVAE(n_char=35, z_dim=196, max_len=120)
    chem.load_state_dict({k[10:]: v for k, v in ckpt.items() if k.startswith("model_con.")})
    flow.eval()
    chem.eval()
    return flow, chem, len(ckpt)


def main():
    lincs = MODEL / "assets" / "models_lincs_PerturbNet_holdout1_cinn_iso_clean_100ep_100ep_model_params.pt"
    sciplex = MODEL / "assets" / "models_sciplex_PerturbNet_holdout_1_cinn_100ep_model_params.pt"
    chemvae_path = MODEL / "assets" / "pretrained_model_chemicalVAE_model_params_525.pt"

    flow_l, chem_l, n_l = load_cinn(lincs, 196)
    flow_s, chem_s, n_s = load_cinn(sciplex, 200)

    chem_pretrained = ChemicalVAE(n_char=35, z_dim=196, max_len=120)
    chem_pretrained.load_state_dict(torch.load(chemvae_path, map_location="cpu", weights_only=False))
    chem_pretrained.eval()

    x = torch.randn(2, 10, 1, 1)
    cond_l = torch.randn(2, 196)
    cond_s = torch.randn(2, 200)
    with torch.no_grad():
        z_l, ld_l = flow_l(x, cond_l)
        z_s, ld_s = flow_s(x, cond_s)
        smiles = torch.zeros(2, 120, 35)
        smiles[:, :, 0] = 1.0
        _, mu, logvar, z_chem = chem_pretrained(smiles)

    np.save(OUT / "pred_latent_lincs.npy", z_l.squeeze(-1).squeeze(-1).numpy())
    np.save(OUT / "pred_latent_sciplex.npy", z_s.squeeze(-1).squeeze(-1).numpy())
    meta = {
        "model": "PerturbNet",
        "lincs_tensors": n_l,
        "sciplex_tensors": n_s,
        "lincs_flow_output": list(z_l.shape),
        "sciplex_flow_output": list(z_s.shape),
        "chemicalvae_z": list(z_chem.shape),
        "chemicalvae_mu": list(mu.shape),
        "chemicalvae_params": sum(p.numel() for p in chem_pretrained.parameters()),
        "note": "checkpoint load + latent forward; full inference needs Python 3.7/TF1.15",
    }
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
