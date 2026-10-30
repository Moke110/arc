#!/usr/bin/env python
"""Squidiff setup smoke: 2,000-gene drug-structure-conditioned forward pass.

The drug condition is the upstream encoding: an rFCFP (Morgan feature)
fingerprint computed from a real SMILES string and scaled by log10(dose+1).
No checkpoint is released, so the network itself is randomly initialised.
"""

import sys

# Keep `pert_model/*/source/` byte-identical to its download (ADR 0030): this
# smoke imports model code from the vendored tree, so byte-code must not land
# there. Set before any other import; the run.sh wrapper exports the same flag.
sys.dont_write_bytecode = True

import json
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent


def find_repo_root(start: Path) -> Path:
    """Walk up to the directory containing .git (ADR 0030 §4)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


ROOT = find_repo_root(HERE)
MODEL = ROOT / "pert_model" / "Squidiff" / "source"
OUT = HERE / "output"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MODEL))

from Squidiff.script_util import create_model  # noqa: E402


def rfcfp_drug_dose(smiles, dose, n_bits=1024):
    """Mirror Squidiff's Drug_dose_encoder for comb_num=1."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, 2, useFeatures=True, nBits=n_bits
    ).ToBitString()
    return np.array(list(fp), dtype=np.float32) * np.log10(dose + 1.0)


def main():
    torch.manual_seed(0)
    n_genes = 2000
    model = create_model(
        n_genes,
        num_layers=2,
        output_dim=n_genes,
        use_encoder=True,
        use_drug_structure=True,
        drug_dimension=1024,
    )
    model.eval()

    # Dexamethasone at two doses: the drug input differs only by the
    # log10(dose+1) scale, which is how Squidiff encodes dose.
    smiles = "CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C"
    doses = [1.0, 100.0]
    drug_vectors = np.stack([rfcfp_drug_dose(smiles, d) for d in doses])
    assert np.count_nonzero(drug_vectors[0]) > 0, "empty rFCFP fingerprint"
    assert not np.allclose(drug_vectors[0], drug_vectors[1]), "dose scaling lost"

    x = torch.randn(2, n_genes)
    t = torch.tensor([1, 2])
    out = model(
        x,
        t,
        x_start=torch.randn(2, n_genes),
        drug_dose=torch.from_numpy(drug_vectors),
        control_feature=torch.randn(2, n_genes),
    )
    np.save(OUT / "pred.npy", out.detach().cpu().numpy())
    meta = {
        "model": "Squidiff",
        "n_genes": n_genes,
        "params": sum(p.numel() for p in model.parameters()),
        "pred_shape": list(out.shape),
        "weights": "random initialization (no released checkpoint)",
        "drug_smiles": smiles,
        "doses": doses,
        "drug_input": "rFCFP (Morgan feature fingerprint, 1024 bits) x log10(dose+1)",
        "note": "real SMILES/dose encoding; network weights are random",
    }
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
