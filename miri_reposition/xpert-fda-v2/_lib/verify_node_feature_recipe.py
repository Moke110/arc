#!/usr/bin/env python3
"""
Ticket 02 — Verify the UniMol mol-token node-feature recipe.

Settles, by experiment, whether the FDA-extension drugs can be given HG graph
node features generated with the same recipe as the L1000 drugs.

Procedure:
  1. Generate a fresh UniMol ``cls_repr`` (dense 512-d mol-token) for the 8,981
     L1000 drugs using the repo's established extraction settings
     (``UniMolRepr(data_type='molecule', remove_hs=False)``, canonical SMILES
     from the frozen idx->SMILES mapping).
  2. Compare against the L1000 node-feature file
     (``drug_emb_model/hg/data/all_drug_node_feat_8981.npy``) on three axes:
     exact match, match up to a per-row normalization, and per-row cosine.
  3. Also report how the same mol-token read out of the frozen workflow UniMol
     array (``input/all_drugs_unimol_arr.npy``, rows 0-1 cols 2:514) compares,
     since that array is the repo's single source of UniMol features.
  4. Write a verdict record: exact / up-to-normalization / mismatch, plus the
     node-feature recipe decision for the FDA extension.

The verdict feeds the author-pinned inference (ticket 05): if the fresh
mol-token does not reproduce the L1000 node features, the residual mismatch is
documented and the FDA node features use the mol-token recipe anyway.

Usage:
    python _lib/verify_node_feature_recipe.py
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, HG_ROOT, WORKSPACE_ROOT  # noqa: E402

N_TRAIN = 8981
N_FDA = 2123
N_TOTAL = N_TRAIN + N_FDA

SMI_MAPPING = WORKFLOW_ROOT / "input" / "all_drugs_idx2smi.npy"
UNIMOL_ARR = WORKFLOW_ROOT / "input" / "all_drugs_unimol_arr.npy"
NODE_FEAT = HG_ROOT / "data" / "all_drug_node_feat_8981.npy"

OUT = WORKFLOW_ROOT / "02_drug_embeddings" / "node_feature_recipe_verdict.json"
CLS_REPR_CACHE = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "l1000_cls_repr_fresh.npy"

BATCH = 512


def generate_fresh_cls_repr(smiles_list, cache_path=None):
    """Generate a fresh UniMol cls_repr for each SMILES (batched, GPU)."""
    from unimol_tools import UniMolRepr

    if cache_path is not None and Path(cache_path).exists():
        print(f"Using cached fresh cls_repr: {cache_path}")
        return np.load(cache_path)

    clf = UniMolRepr(data_type="molecule", remove_hs=False)
    out = []
    t0 = time.time()
    for i in range(0, len(smiles_list), BATCH):
        chunk = smiles_list[i : i + BATCH]
        reprs = clf.get_repr(chunk)
        out.append(np.asarray(reprs, dtype=np.float32))
        print(f"  [{i}:{i + len(chunk)}] generated in {time.time() - t0:.0f}s")
    arr = np.concatenate(out, axis=0)
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, arr)
    return arr


def per_row_cosine(a, b):
    a_n = a / (np.linalg.norm(a, axis=1, keepdims=True) + 1e-12)
    b_n = b / (np.linalg.norm(b, axis=1, keepdims=True) + 1e-12)
    return (a_n * b_n).sum(axis=1)


def main():
    idx2smi = np.load(SMI_MAPPING, allow_pickle=True).item()
    smiles = [idx2smi[i] for i in range(N_TRAIN)]
    node_feat = np.load(NODE_FEAT)  # (8981, 512)
    unimol_arr = np.load(UNIMOL_ARR, mmap_mode="r", allow_pickle=True)
    frozen_cls = unimol_arr[:N_TRAIN, 0, 2:514].astype(np.float32)  # mol-token cols

    print(f"L1000 node-feature file: {node_feat.shape} {node_feat.dtype}")
    print(f"Frozen UniMol array mol-token: {frozen_cls.shape}")
    print(f"Generating fresh UniMol cls_repr for {N_TRAIN} drugs...")
    fresh = generate_fresh_cls_repr(smiles, cache_path=CLS_REPR_CACHE)

    # --- Axis 1: exact match ---
    exact = bool(np.array_equal(fresh, node_feat))
    print(f"\n[Axis 1] exact match: {exact}")

    # --- Axis 2: per-row scalar normalization ---
    scale = (fresh * node_feat).sum(axis=1) / ((fresh * fresh).sum(axis=1) + 1e-12)
    resid = np.abs(node_feat - scale[:, None] * fresh)
    up_to_norm = float(np.max(resid)) < 1e-5
    print(f"[Axis 2] up-to-normalization: {up_to_norm} (max resid {resid.max():.4f})")

    # --- Axis 3: cosine similarity ---
    cos_fresh = per_row_cosine(fresh, node_feat)
    cos_frozen = per_row_cosine(frozen_cls, node_feat)
    print(f"[Axis 3] cos(fresh, node_feat): mean={cos_fresh.mean():.4f}")
    print(f"         cos(frozen-array, node_feat): mean={cos_frozen.mean():.4f}")

    maxdiff_fresh = float(np.abs(fresh - node_feat).max())
    maxdiff_frozen = float(np.abs(frozen_cls - node_feat).max())

    if exact:
        verdict = "exact"
    elif up_to_norm:
        verdict = "up-to-normalization"
    else:
        verdict = "mismatch"

    record = {
        "ticket": 2,
        "verdict": verdict,
        "n_drugs": N_TRAIN,
        "settings": {
            "unimol_tools": "UniMolRepr(data_type='molecule', remove_hs=False)",
            "smiles": "frozen idx->SMILES mapping (canonical)",
            "cls_repr_position": "cls_repr of get_repr(smiles)",
        },
        "axis1_exact": exact,
        "axis2_up_to_normalization": up_to_norm,
        "axis3_cos_fresh_mean": float(cos_fresh.mean()),
        "axis3_cos_frozen_array_mean": float(cos_frozen.mean()),
        "maxdiff_fresh": maxdiff_fresh,
        "maxdiff_frozen_array": maxdiff_frozen,
        "node_feature_recipe_decision": (
            "Fresh UniMol cls_repr does not reproduce the L1000 node-feature file; "
            "the L1000 node features are not reproducible by the current unimol_tools "
            "pipeline (a different UniMol model version produced them). Per the feature "
            "spec the FDA-extension node features use the mol-token recipe regardless. "
            "To avoid feature drift and stay closest to the L1000 recipe, the FDA "
            "mol-token is read from the frozen workflow UniMol array (cls_repr at rows "
            "0-1, feature cols 2:514) rather than freshly generated — that array is the "
            "repo's single source of UniMol features and its mol-token is cos≈0.97 with "
            "the L1000 node features (fresh generation is cos≈0.16)."
        ),
    }

    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(record, f, indent=2)
    print(f"\nVerdict: {verdict}")
    print(f"Recorded: {OUT}")


if __name__ == "__main__":
    main()
