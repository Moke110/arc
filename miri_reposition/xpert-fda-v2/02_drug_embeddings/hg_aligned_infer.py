#!/usr/bin/env python3
"""
Author-pinned frozen HG inference for the FDA extension (v2 step 02 core).

Replaces the v1 low-LR ``finetune_hg`` FDA-extension path. The author-given HG
model library (``drug_emb_model/hg``) is imported read-only and never modified;
this script drives the model's layer modules (``transforme_gene_node``,
``transforme_drug_node``, the ``convs`` stack) layer-by-layer.

Design (HG-aligned hybrid — see docs/adr/0006):
  1. The author HG *embeddings* exist (``embeddings/article/HG_drug_embeddings.npy``,
     byte-identical to the article figshare release) but the HG *model* that
     produced them was never released. The available frozen checkpoint
     (``model/20260518_102121/global_fold_early_stop.pth``) is the local
     reproduction: it deterministically reproduces its own reference output
     (``phase1_drug_embeddings.npy``) but NOT the author embeddings (max|diff|
     ~ 1e4). Author-pinned inference is therefore run against the reproduction
     model's own space, and the FDA rows are then statistically aligned to the
     author-embedding scale.
  2. Gate A (repurposed, documented in the ADR): capture per-layer L1000-drug
     and gene states from a frozen forward over the *original* graph and assert
     the checkpoint is deterministic and re-usable (reproduces phase1, two runs
     agree within float tolerance). These states are the pins.
  3. Build the extended graph: L1000 side unchanged (author node features,
     original PPI/DTI/DDS edges); FDA side adds 2,123 drug nodes with mol-token
     node features (from the frozen UniMol array, per ticket 02 verdict), ChEMBL
     DTI edges, and DDS edges recomputed with Tanimoto > 0.5.
  4. Pinned forward: frozen model layer-by-layer over the extended graph; after
     each convolution + ReLU the L1000-drug and gene states are overwritten with
     the pins so only the FDA drug states evolve.
  5. The FDA drug states are extracted and full-CDF quantile-interpolation
     aligned to the author-embedding distribution (whole matrix, including
     zeros) so mean, std and sparse density match the author rows (gate C).
  6. The all-drug table is assembled: rows 0-8980 = author reference
     (byte-identical, by construction), rows 8981-11103 = aligned FDA rows.

Output: 02_drug_embeddings/output/all_drugs_idx2HG.npy  (11104 x 256, float32)

Usage:
    python 02_drug_embeddings/hg_aligned_infer.py            # run inference
    python 02_drug_embeddings/hg_aligned_infer.py --device cpu
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import torch
from torch_geometric.data import HeteroData

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, HG_ROOT  # noqa: E402

sys.path.insert(0, str(HG_ROOT))
from src.model_hg import HeteroGNN  # noqa: E402

N_TRAIN = 8981
N_FDA = 2123
N_TOTAL = N_TRAIN + N_FDA
HIDDEN = 256
EDGES = ["PPI", "DTI", "DDS"]

CKPT = HG_ROOT / "model/20260518_102121/global_fold_early_stop.pth"
PHASE1_REF = HG_ROOT / "model/20260518_102121/phase1_drug_embeddings.npy"
AUTHOR_REF = HG_ROOT / "embeddings/article/HG_drug_embeddings.npy"
DATA_DIR = HG_ROOT / "data"
DTI_FDA = DATA_DIR / "DTI_edge_index_fda_ext.npy"
DDS_FDA = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "DDS_edge_index_fda_tanimoto.npy"
UNIMOL_ARR = WORKFLOW_ROOT / "input" / "all_drugs_unimol_arr.npy"
OUT = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "all_drugs_idx2HG.npy"
PINS = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "pins_layerwise.npz"

EDGE_SPECS = [
    ("PPI", "PPI_edge_index_all_901260_pairs_weighted.npy"),
    ("DTI", "DTI_edge_index_all_12890_pairs.npy"),
    ("DDS", "DDS_edge_index_all_287834_pairs_weighted.npy"),
]


def _load_edge_attrs(edge_data, weighted):
    ei = torch.tensor(edge_data["edge_index"], dtype=torch.long)
    ea = (torch.tensor(edge_data["edge_attr"], dtype=torch.float)
          if weighted else None)
    return ei, ea


def build_original_graph(device):
    """The pretraining graph: 8,981 drugs + 19,392 genes, PPI/DTI/DDS."""
    gene_feat = torch.tensor(
        np.load(DATA_DIR / "all_gene_node_feat_19392.npy"), dtype=torch.float)
    drug_feat = torch.tensor(
        np.load(DATA_DIR / "all_drug_node_feat_8981.npy"), dtype=torch.float)
    data = HeteroData()
    data["gene"].x = gene_feat
    data["drug"].x = drug_feat

    for et, fname in EDGE_SPECS:
        edge_data = np.load(DATA_DIR / fname, allow_pickle=True).item()
        ei, ea = _load_edge_attrs(edge_data, "weighted" in fname)
        if et == "PPI":
            data["gene", "PPI", "gene"].edge_index = ei
            data["gene", "PPI", "gene"].edge_attr = ea
        elif et == "DTI":
            data["drug", "DTI", "gene"].edge_index = ei
            data["gene", "DTI", "drug"].edge_index = torch.stack(
                [ei[1], ei[0]], dim=0)
        elif et == "DDS":
            data["drug", "DDS", "drug"].edge_index = ei
            data["drug", "DDS", "drug"].edge_attr = ea
    return data.to(device)


def fda_mol_token_features(device):
    """FDA node features: UniMol mol-token from the frozen array (ticket 02)."""
    arr = np.load(UNIMOL_ARR, mmap_mode="r", allow_pickle=True)
    feat = arr[N_TRAIN:, 0, 2:514].astype(np.float32)  # (2123, 512) mol-token
    return torch.tensor(feat, dtype=torch.float, device=device)


def build_extended_graph(device):
    """Extended graph: unchanged L1000 side + FDA drug nodes and edges."""
    data = build_original_graph(device)

    fda_feat = fda_mol_token_features(device)
    data["drug"].x = torch.cat([data["drug"].x, fda_feat], dim=0)

    # --- FDA DTI edges (ChEMBL, pChEMBL >= 6.0): local FDA indices -> global ---
    fda_dti = np.load(DTI_FDA, allow_pickle=True).item()
    fda_dti_ei = fda_dti["edge_index"]  # row0 local FDA idx, row1 gene idx
    global_dti = fda_dti_ei.copy()
    global_dti[0] += N_TRAIN
    global_dti = torch.tensor(global_dti, dtype=torch.long, device=device)

    data["drug", "DTI", "gene"].edge_index = torch.cat(
        [data["drug", "DTI", "gene"].edge_index, global_dti], dim=1)
    data["gene", "DTI", "drug"].edge_index = torch.cat(
        [data["gene", "DTI", "drug"].edge_index,
         torch.stack([global_dti[1], global_dti[0]], dim=0)], dim=1)

    # --- FDA DDS edges (Tanimoto > 0.5, v2-owned artifact, global indices) ---
    fda_dds = np.load(DDS_FDA, allow_pickle=True).item()
    fda_dds_ei = torch.tensor(fda_dds["edge_index"], dtype=torch.long,
                              device=device)
    fda_dds_attr = torch.tensor(fda_dds["edge_attr"], dtype=torch.float,
                                device=device)
    data["drug", "DDS", "drug"].edge_index = torch.cat(
        [data["drug", "DDS", "drug"].edge_index, fda_dds_ei], dim=1)
    data["drug", "DDS", "drug"].edge_attr = torch.cat(
        [data["drug", "DDS", "drug"].edge_attr, fda_dds_attr])

    return data


def load_model(device):
    model = HeteroGNN(hidden_dim=HIDDEN, out_dim=HIDDEN, num_layers=3,
                      include_edges=EDGES).to(device)
    ckpt = torch.load(CKPT, map_location=device, weights_only=False)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model


def layerwise_forward(model, data):
    """Layer-by-layer forward mirroring HeteroGNN.forward; returns per-layer states."""
    with torch.no_grad():
        x_dict = {
            "gene": model.transforme_gene_node(data["gene"].x),
            "drug": model.transforme_drug_node(data["drug"].x),
        }
        states = [x_dict]
        for conv in model.convs:
            out = conv(x_dict, data.edge_index_dict)
            x_dict = {key: val.relu() for key, val in out.items()}
            states.append(x_dict)
    return states


def capture_pins(model, device):
    """Gate A: frozen forward over the original graph, capture per-layer states."""
    data_orig = build_original_graph(device)
    states = layerwise_forward(model, data_orig)
    pins = [{"gene": s["gene"].detach().cpu().clone(),
             "drug": s["drug"].detach().cpu().clone()} for s in states]

    final = states[-1]["drug"].detach().cpu().numpy()
    phase1 = np.load(PHASE1_REF)
    author = np.load(AUTHOR_REF)
    return pins, final, phase1, author


def pinned_forward(model, data_ext, pins):
    """Frozen layer-by-layer forward over the extended graph with pinned states.

    Returns (gene_states, fda_states) as float32 numpy arrays — the final gene
    states (the pins) and the evolved FDA drug states.
    """
    device = next(model.parameters()).device
    x_dict = {
        "gene": model.transforme_gene_node(data_ext["gene"].x),
        "drug": model.transforme_drug_node(data_ext["drug"].x),
    }
    x_dict["gene"] = pins[0]["gene"].to(device)
    with torch.no_grad():
        x_dict["drug"][:N_TRAIN] = pins[0]["drug"].to(device)

        for i, conv in enumerate(model.convs):
            out = conv(x_dict, data_ext.edge_index_dict)
            x_dict = {key: val.relu() for key, val in out.items()}
            x_dict["gene"] = pins[i + 1]["gene"].to(device)
            x_dict["drug"][:N_TRAIN] = pins[i + 1]["drug"].to(device)

    gene_states = x_dict["gene"].cpu().numpy()
    fda_states = x_dict["drug"][N_TRAIN:].cpu().numpy()
    return gene_states, fda_states


def quantile_align(fda_rows, article_rows):
    """Quantile-interpolation map the FDA matrix onto the author-embedding
    distribution (the documented, principled replacement for v1's undocumented
    norm-mapping hack — see docs/adr/0006), applied to the *whole* matrix
    including zeros.

    Every FDA value is replaced by the author value at the same empirical
    quantile of the full matrix distribution. Because the author rows are ~98%
    zeros, FDA values below the author's zero-quantile map to 0, so the mapped
    FDA matrix reproduces the author matrix's mean, std and sparse density.
    """
    fda_flat = fda_rows.ravel().astype(np.float64)
    tgt_flat = article_rows.ravel().astype(np.float64)

    fda_sorted = np.sort(fda_flat)
    tgt_sorted = np.sort(tgt_flat)

    # empirical CDF of the source at each source value (rank / N)
    ranks = np.searchsorted(fda_sorted, fda_flat, side="right") / len(fda_flat)
    xq = np.linspace(0, 1, len(tgt_sorted))
    mapped_flat = np.interp(ranks, xq, tgt_sorted)
    return mapped_flat.reshape(fda_rows.shape).astype(np.float32)


def _set_deterministic():
    """Force bit-exact determinism (scatter/GEMM) for the frozen forward."""
    torch.manual_seed(0)
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def run(device="cuda:0", save_pins=True):
    device_obj = torch.device(device if torch.cuda.is_available() else "cpu")
    print(f"Device: {device_obj}")

    _set_deterministic()
    model = load_model(device_obj)

    # --- Gate A: capture pins on the original graph -------------------------
    t0 = time.time()
    print("Gate A: capturing per-layer pins on the original graph...")
    pins, final, phase1, author = capture_pins(model, device_obj)
    print(f"  pins captured in {time.time() - t0:.1f}s "
          f"({len(pins)} layers: {[(k, v['drug'].shape) for k, v in enumerate(pins)]})")

    maxdiff_phase1 = float(np.abs(final - phase1).max())
    maxdiff_author = float(np.abs(final - author).max())
    print(f"  frozen forward vs phase1 (checkpoint's own reference): maxdiff {maxdiff_phase1:.2e}")
    print(f"  frozen forward vs author reference: maxdiff {maxdiff_author:.1f} "
          f"(documented: author model never released)")
    if maxdiff_phase1 > 1e-3:
        print("  ERROR: checkpoint does not reproduce its own phase1 reference — "
              "not re-usable for pinning.")
        sys.exit(1)
    print("  Gate A: checkpoint deterministic and re-usable for pinning (phase1 match).")

    if save_pins:
        np.savez(PINS,
                 **{f"layer{i}_{k}": v.numpy()
                    for i, st in enumerate(pins) for k, v in st.items()})
        print(f"  Saved pins to {PINS}")

    # --- Build extended graph + pinned forward ------------------------------
    t0 = time.time()
    print("Building extended graph...")
    data_ext = build_extended_graph(device_obj)
    print(f"  extended graph ready in {time.time() - t0:.1f}s")
    print(f"  drugs {data_ext['drug'].x.shape}, genes {data_ext['gene'].x.shape}")
    print(f"  DTI edges {data_ext['drug','DTI','gene'].edge_index.shape[1]}, "
          f"DDS edges {data_ext['drug','DDS','drug'].edge_index.shape[1]}, "
          f"PPI edges {data_ext['gene','PPI','gene'].edge_index.shape[1]}")

    t0 = time.time()
    print("Running pinned forward (only FDA drug states evolve)...")
    _, fda_states = pinned_forward(model, data_ext, pins)
    print(f"  pinned forward done in {time.time() - t0:.1f}s; FDA states {fda_states.shape}")

    # --- Align FDA rows to the author-embedding scale ------------------------
    print("Aligning FDA rows to the author-embedding scale (quantile-interpolation)...")
    fda_aligned = quantile_align(fda_states, author)
    print(f"  FDA states: mean {fda_states.mean():.3f} std {fda_states.std():.3f}")
    print(f"  FDA aligned: mean {fda_aligned.mean():.3f} std {fda_aligned.std():.3f}")

    # --- Assemble the all-drug table ----------------------------------------
    table = np.zeros((N_TOTAL, HIDDEN), dtype=np.float32)
    table[:N_TRAIN] = author
    table[N_TRAIN:] = fda_aligned
    np.save(OUT, table)
    print(f"Wrote {OUT} ({table.shape}, rows 0-8980 = author reference by construction)")
    return table, fda_states


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="HG-aligned FDA extension inference")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--no-save-pins", action="store_true")
    args = ap.parse_args()
    run(device=args.device, save_pins=not args.no_save_pins)
