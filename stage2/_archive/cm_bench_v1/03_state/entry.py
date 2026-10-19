#!/usr/bin/env python3
"""
State (ST-HVG-Tahoe) entrant entry script — runs inside `arc-stage2-v1-state`.

Contract (ticket 11 / _lib/adapter.py):
    python entry.py <train.npz> <test.npz> <out.npz> [--finetune-steps N]

Reads the fold's `Cases` bundles (baseline log2-CPM, drug, dose, Δ target on train)
plus the sidecar JSON, predicts a Δ matrix for the test wells and writes it aligned
to the common gene universe.

Faithful-to-the-checkpoint, with documented approximations (read 03_state/README.md):

  * Space: ST-HVG-Tahoe predicts a 2000-dim HVG vector. The 2000 HVG gene identities
    are not published; they were recovered by exact value-matching against the
    `arcinstitute/State-Tahoe-Filtered` training shards (see
    `hvg_gene_names_recovered.json`). HVG dims that could not be mapped (or whose gene
    is outside the benchmark universe) are dropped.
  * Input: the model was trained on log1p(raw counts); the benchmark only froze
    log2(CPM+1). We feed the benchmark baseline directly (scale approximation).
  * Perturbation: the checkpoint uses a one-hot over 1138 Tahoe drug-dose combos.
    A CM drug is mapped to the closest Tahoe dose for the same base name; CM drugs
    absent from Tahoe fall back to the DMSO control one-hot (=> near-zero Δ). This is
    the model's own cold-drug limitation, reported as-is, never fabricated.
  * Finetune: a short domain-adaptation finetune on the fold's control baselines
    (State's own test-time-finetune path), never from scratch.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent
CONTROL_PERT = "[('DMSO_TF', 0.0, 'uM')]"

import os

STATE_WEIGHTS = Path(os.environ.get(
    "STATE_WEIGHTS", "/home/chang/arc/pert_model/state/weights/ST-HVG-Tahoe"))
CKPT = Path(os.environ.get("STATE_CKPT", STATE_WEIGHTS / "checkpoints" / "final.ckpt"))
PERT_MAP = Path(os.environ.get("STATE_PERT_MAP", STATE_WEIGHTS / "pert_onehot_map.pt"))
HVG_MAP_PATH = Path(os.environ.get("STATE_HVG_MAP", HERE / "hvg_gene_names_recovered.json"))


def load_cases(npz_path):
    z = np.load(npz_path)
    meta = json.loads(Path(npz_path).with_suffix(".json").read_text())
    return {"baseline": z["baseline"], "dose_log": z["dose_log"],
            "delta": z["delta"] if "delta" in z else None, **meta}


def universe_and_hvg_index():
    universe = (WORKFLOW / "output" / "bench" / "universe.txt").read_text().split()
    uni_col = {g: i for i, g in enumerate(universe)}
    hvg_names = json.loads(HVG_MAP_PATH.read_text())["hvg_gene_names"]
    hvg_to_uni = np.full(len(hvg_names), -1, dtype=int)
    for j, g in enumerate(hvg_names):
        if g and g in uni_col:
            hvg_to_uni[j] = uni_col[g]
    return universe, hvg_to_uni


def to_hvg(baseline, hvg_to_uni):
    """(n, G) universe matrix -> (n, 2000) State HVG input (unmapped dims = 0)."""
    X = np.zeros((baseline.shape[0], len(hvg_to_uni)), dtype=np.float32)
    m = hvg_to_uni >= 0
    X[:, m] = baseline[:, hvg_to_uni[m]]
    return X


def norm_name(s):
    s = str(s).lower().strip()
    s = re.sub(r"\s*\(.*?\)", "", s)          # drop salt/hydrate parentheticals
    s = s.replace("_", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s).strip()
    return s


def parse_tahoe_key(k):
    m = re.match(r"\[\('(.+?)',\s*([0-9.eE+-]+),\s*'uM'\)\]", str(k))
    return (norm_name(m.group(1)), float(m.group(2))) if m else None


def build_tahoe_lookup(pert_map):
    lookup = {}
    for k in pert_map.keys():
        pk = parse_tahoe_key(k)
        if pk is None:
            continue
        name, conc = pk
        lookup.setdefault(name, []).append((conc, k))
    return lookup


def pert_vector(drug, dose_log, lookup, pert_map, base_dim):
    cand = lookup.get(norm_name(drug))
    if cand:
        best = min(cand, key=lambda ck: abs(np.log10(max(ck[0], 1e-6)) - dose_log))
        return best[1], pert_map[best[1]].float()
    return None, pert_map[CONTROL_PERT].float()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("train_npz")
    ap.add_argument("test_npz")
    ap.add_argument("out_npz")
    ap.add_argument("--finetune-steps", type=int, default=150)
    ap.add_argument("--finetune-lr", type=float, default=1e-5)
    ap.add_argument("--batch", type=int, default=128)
    args = ap.parse_args()

    torch.manual_seed(42)
    np.random.seed(42)

    from state.tx.models.state_transition import StateTransitionPerturbationModel

    universe, hvg_to_uni = universe_and_hvg_index()
    n_uni = len(universe)
    n_hvg = len(hvg_to_uni)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = StateTransitionPerturbationModel.load_from_checkpoint(str(CKPT))
    model = model.to(device).eval()
    pert_map = torch.load(PERT_MAP, map_location="cpu", weights_only=False)
    lookup = build_tahoe_lookup(pert_map)

    train = load_cases(args.train_npz)
    test = load_cases(args.test_npz)
    Xtr = to_hvg(train["baseline"], hvg_to_uni)
    Xte = to_hvg(test["baseline"], hvg_to_uni)

    # ---- optional domain-adaptation finetune on the fold's control baselines ----
    if args.finetune_steps > 0 and len(Xtr):
        opt = torch.optim.Adam(model.parameters(), lr=args.finetune_lr)
        ctrl = pert_map[CONTROL_PERT].float().to(device)
        Xtr_t = torch.tensor(Xtr, device=device)
        n = Xtr_t.shape[0]
        model.train()
        for step in range(args.finetune_steps):
            idx = torch.randint(0, n, (min(args.batch, n),))
            xb = Xtr_t[idx]
            batch = {"ctrl_cell_emb": xb,
                     "pert_emb": ctrl.unsqueeze(0).repeat(xb.shape[0], 1),
                     "pert_name": [CONTROL_PERT] * xb.shape[0]}
            opt.zero_grad()
            pred = model.forward(batch, padded=False)
            loss = torch.nn.functional.mse_loss(pred, xb)
            loss.backward()
            opt.step()
        model.eval()
        print(f"[state] finetuned {args.finetune_steps} steps on {n} train baselines "
              f"(lr={args.finetune_lr})", file=sys.stderr)

    # ---- predict test wells ----
    Xte_t = torch.tensor(Xte, device=device)
    n_test = Xte_t.shape[0]
    delta_uni = np.zeros((n_test, n_uni), dtype=np.float32)
    seen = 0
    with torch.no_grad():
        for start in range(0, n_test, args.batch):
            end = min(start + args.batch, n_test)
            drugs = test["drug"][start:end]
            doses = test["dose_log"][start:end]
            vecs = []
            for d, dl in zip(drugs, doses):
                key, v = pert_vector(d, float(dl), lookup, pert_map, n_hvg)
                if key is not None:
                    seen += 1
                vecs.append(v)
            pert_emb = torch.stack(vecs).to(device)
            batch = {"ctrl_cell_emb": Xte_t[start:end], "pert_emb": pert_emb,
                     "pert_name": list(drugs)}
            pred = model.predict_step(batch, 0, padded=False)["preds"]
            delta_hvg = (pred - Xte_t[start:end]).cpu().numpy()
            m = hvg_to_uni >= 0
            delta_uni[start:end, hvg_to_uni[m]] = delta_hvg[:, m]
    print(f"[state] test wells={n_test}; mapped HVG dims={int((hvg_to_uni>=0).sum())}; "
          f"seen-Tahoe drug wells={seen}", file=sys.stderr)

    np.savez_compressed(args.out_npz, pred=delta_uni.astype(np.float32))


if __name__ == "__main__":
    main()
