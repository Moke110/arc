#!/usr/bin/env python3
"""
Gate harness for the v2 HG-aligned FDA extension (the single testing seam).

Mirrors the existing ``verify_parity.py`` seam (docs/adr/0004): verification is
a script + manual output inspection, no unit-test framework. Asserts the four
correctness gates after ``02_drug_embeddings`` produces the all-drug table.

  Gate A — the frozen checkpoint is deterministic and re-usable: a frozen
           forward over the original graph reproduces the checkpoint's own
           reference (``phase1_drug_embeddings.npy``) within float tolerance,
           and two pin-capture runs agree within float tolerance. (Repurposed:
           the author HG model was never released, so reproducing the author
           reference itself is impossible — documented in docs/adr/0006.)
  Gate B — output rows 0-8980 are byte-identical to the author reference.
  Gate C — FDA/train std ratio ~= 1.0 (after quantile alignment; no v1-style
           ad-hoc norm-mapping of L1000 rows).
  Gate D — FDA DTI drug-gene cosine-vs-random ratio >= 2.0.

Gates A and B abort on failure (exit 1). Gates C and D fail loudly but may be
reviewed if borderline (a badly-off C triggers a revisit of the node-feature
recipe, never silent acceptance).

Usage:
    python _lib/verify_gates.py                    # run all gates
    python _lib/verify_gates.py --table <path>     # custom output table
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, HG_ROOT  # noqa: E402

N_TRAIN = 8981
N_TOTAL = 11104
N_GENE = 19392

AUTHOR_REF = HG_ROOT / "embeddings/article/HG_drug_embeddings.npy"
PHASE1_REF = HG_ROOT / "model/20260518_102121/phase1_drug_embeddings.npy"
CKPT = HG_ROOT / "model/20260518_102121/global_fold_early_stop.pth"
DTI_FDA = HG_ROOT / "data" / "DTI_edge_index_fda_ext.npy"
PINS = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "pins_layerwise.npz"
TABLE = WORKFLOW_ROOT / "02_drug_embeddings" / "output" / "all_drugs_idx2HG.npy"

results = []


def check(name, cond, detail=""):
    ok = bool(cond)
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def gate_a():
    print("\nGate A — frozen forward deterministic and re-usable (checkpoint's own reference)")
    import torch

    sys.path.insert(0, str(HG_ROOT))
    sys.path.insert(0, str(WORKFLOW_ROOT / "02_drug_embeddings"))
    from src.model_hg import HeteroGNN
    from hg_aligned_infer import build_original_graph, layerwise_forward

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = HeteroGNN(hidden_dim=256, out_dim=256, num_layers=3,
                      include_edges=["PPI", "DTI", "DDS"]).to(device)
    ckpt = torch.load(CKPT, map_location=device, weights_only=False)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()

    data = build_original_graph(device)
    states1 = layerwise_forward(model, data)
    states2 = layerwise_forward(model, data)
    final = states1[-1]["drug"].detach().cpu().numpy()

    d_self = float(np.abs(states1[-1]["drug"].detach().cpu().numpy()
                          - states2[-1]["drug"].detach().cpu().numpy()).max())
    phase1 = np.load(PHASE1_REF)
    d_phase1 = float(np.abs(final - phase1).max())

    ok = d_phase1 <= 1e-3 and d_self <= 1e-3
    check("frozen forward reproduces phase1 (maxdiff <= 1e-3)", ok,
          f"maxdiff vs phase1 {d_phase1:.2e}, run-to-run {d_self:.2e}")
    return ok


def gate_b(table_path):
    print("\nGate B — rows 0-8980 byte-identical to the author reference")
    table = np.load(table_path, mmap_mode="r")
    author = np.load(AUTHOR_REF)
    ok = table.shape[0] == N_TOTAL and np.array_equal(
        np.asarray(table[:N_TRAIN]), np.asarray(author))
    check("rows 0-8980 byte-identical to author reference", ok,
          f"table {table.shape}, maxdiff {np.abs(np.asarray(table[:N_TRAIN]) - author).max()}")
    return ok


def gate_c(table_path):
    print("\nGate C — FDA/train std ratio ~= 1.0 (no ad-hoc L1000 norm mapping)")
    table = np.load(table_path, mmap_mode="r")
    train_std = float(np.asarray(table[:N_TRAIN]).std())
    fda_std = float(np.asarray(table[N_TRAIN:]).std())
    ratio = fda_std / (train_std + 1e-10)
    ok = abs(ratio - 1.0) <= 0.05
    check("FDA/train std ratio ~= 1.0", ok,
          f"train_std {train_std:.3f}, fda_std {fda_std:.3f}, ratio {ratio:.4f}")
    return ok


def gate_d():
    print("\nGate D — FDA DTI drug-gene cosine vs random ratio >= 2.0")
    import torch

    sys.path.insert(0, str(HG_ROOT))
    sys.path.insert(0, str(WORKFLOW_ROOT / "02_drug_embeddings"))
    from src.model_hg import HeteroGNN
    from hg_aligned_infer import build_extended_graph

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = HeteroGNN(hidden_dim=256, out_dim=256, num_layers=3,
                      include_edges=["PPI", "DTI", "DDS"]).to(device)
    ckpt = torch.load(CKPT, map_location=device, weights_only=False)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()

    pins = np.load(PINS)
    data_ext = build_extended_graph(device)
    with torch.no_grad():
        x = data_ext["gene"].x
        x_dict = {
            "gene": model.transforme_gene_node(data_ext["gene"].x),
            "drug": model.transforme_drug_node(data_ext["drug"].x),
        }
        x_dict["drug"][:N_TRAIN] = torch.tensor(
            pins["layer0_drug"], dtype=torch.float, device=device)
        for i, conv in enumerate(model.convs):
            out = conv(x_dict, data_ext.edge_index_dict)
            x_dict = {key: val.relu() for key, val in out.items()}
            x_dict["gene"] = torch.tensor(
                pins[f"layer{i + 1}_gene"], dtype=torch.float, device=device)
            x_dict["drug"][:N_TRAIN] = torch.tensor(
                pins[f"layer{i + 1}_drug"], dtype=torch.float, device=device)
        gene_emb = x_dict["gene"].cpu().numpy()
        fda_emb = x_dict["drug"][N_TRAIN:].cpu().numpy()

    fda_dti = np.load(DTI_FDA, allow_pickle=True).item()
    drug_local = fda_dti["edge_index"][0]
    gene_idx = fda_dti["edge_index"][1]

    def cos(a, b):
        an = a / (np.linalg.norm(a, axis=1, keepdims=True) + 1e-10)
        bn = b / (np.linalg.norm(b, axis=1, keepdims=True) + 1e-10)
        return (an * bn).sum(axis=1)

    dti_cos = cos(fda_emb[drug_local], gene_emb[gene_idx])
    n_random = min(50000, len(dti_cos) * 3)
    rng = np.random.RandomState(42)
    rand_drugs = rng.randint(0, len(fda_emb), n_random)
    rand_genes = rng.randint(0, len(gene_emb), n_random)
    rand_cos = cos(fda_emb[rand_drugs], gene_emb[rand_genes])
    ratio = float(dti_cos.mean() / (rand_cos.mean() + 1e-10))

    ok = ratio >= 2.0
    check("FDA DTI cosim vs random ratio >= 2.0", ok,
          f"n_edges {len(dti_cos)}, dti_mean {dti_cos.mean():.4f}, "
          f"random_mean {rand_cos.mean():.4f}, ratio {ratio:.2f}x")
    return ok


def main():
    ap = argparse.ArgumentParser(description="v2 HG-aligned gate harness")
    ap.add_argument("--table", default=str(TABLE))
    args = ap.parse_args()
    table_path = Path(args.table)

    import torch
    torch.manual_seed(0)
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    print("=" * 60)
    print("v2 HG-aligned gates (A-D)")
    print("=" * 60)

    g_a = gate_a()
    if not g_a:
        print("\nRESULT: Gate A FAILED — aborting before any extension inference.")
        sys.exit(1)
    g_b = gate_b(table_path)
    if not g_b:
        print("\nRESULT: Gate B FAILED — L1000 rows not byte-identical. Aborting.")
        sys.exit(1)
    g_c = gate_c(table_path)
    g_d = gate_d()

    passed = all(results)
    print("\n" + "=" * 60)
    print(f"RESULT: {'ALL GATES PASSED' if passed else 'A/B passed; C/D need review'}")
    print("=" * 60)
    return 0 if g_a and g_b else 1


if __name__ == "__main__":
    sys.exit(main())
