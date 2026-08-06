#!/usr/bin/env python3
"""
FDA DDS edge recomputation with the L1000 pretraining recipe.

v2-owned artifact (user story 05 / ticket 24). Replaces the calibrated-threshold
DDS edge file (``drug_emb_model/hg/data/DDS_edge_index_fda_drugs.npy``) with
edges computed using the pretraining recipe: Morgan (radius-2, 1024-bit)
Tanimoto similarity strictly above 0.5. Edge weight = Tanimoto similarity.
Drug indices are global (0-11103).

The Morgan fingerprints reproduce the L1000 DDS edge weights exactly (verified:
``GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)`` equals the training
fingerprint file shipped by the article's figshare release), so the FDA-side
edges are computed in the same fingerprint space as the training side.

The library's Morgan fingerprint / Tanimoto helpers are imported read-only;
nothing under ``drug_emb_model/hg`` is modified and the L1000 DDS edges are
never touched.

Output: 02_drug_embeddings/output/DDS_edge_index_fda_tanimoto.npy
  dict{edge_index: (2, N) global indices, edge_attr: (N,) Tanimoto weights,
       stats: {...}}

Usage:
    python 02_drug_embeddings/compute_dds_fda_tanimoto.py
"""

import logging
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, HG_ROOT  # noqa: E402

sys.path.insert(0, str(HG_ROOT))
from src.compute_dds_fda import compute_fda_dds_edges  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

N_TRAIN = 8981
N_TOTAL = 11104
THRESHOLD = 0.5
MORGAN_BITS = 1024
RADIUS = 2

SMI_MAPPING = WORKFLOW_ROOT / "input" / "all_drugs_idx2smi.npy"
OUT = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "DDS_edge_index_fda_tanimoto.npy"


def morgan_from_smiles(smiles):
    """1024-bit radius-2 Morgan bit vector (dense float) for one SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(MORGAN_BITS, dtype=np.float32)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=MORGAN_BITS)
    arr = np.zeros(MORGAN_BITS, dtype=np.float32)
    arr[list(fp.GetOnBits())] = 1.0
    return arr


def load_all_morgan():
    idx2smi = np.load(SMI_MAPPING, allow_pickle=True).item()
    mats = np.zeros((N_TOTAL, MORGAN_BITS), dtype=np.float32)
    for i in range(N_TOTAL):
        mats[i] = morgan_from_smiles(idx2smi[i])
    return mats[:N_TRAIN].astype(np.float64), mats[N_TRAIN:].astype(np.float64)


def main():
    print(f"Computing Morgan fingerprints for {N_TOTAL} drugs "
          f"(radius {RADIUS}, {MORGAN_BITS} bits)...")
    morgan_train, morgan_fda = load_all_morgan()
    print(f"  train {morgan_train.shape}, fda {morgan_fda.shape}")

    result = compute_fda_dds_edges(
        morgan_train, morgan_fda,
        threshold=THRESHOLD,
        n_train=N_TRAIN,
    )
    stats = result["stats"]
    print(f"FDA DDS edges (Tanimoto > {THRESHOLD}): "
          f"{stats['n_edges']} total "
          f"({stats['n_fda_train']} FDA<->train, {stats['n_fda_fda']} FDA<->FDA)")

    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    np.save(OUT, result, allow_pickle=True)
    print(f"Saved v2-owned FDA DDS edges to {OUT}")
    print(f"Stats: {stats}")


if __name__ == "__main__":
    main()
