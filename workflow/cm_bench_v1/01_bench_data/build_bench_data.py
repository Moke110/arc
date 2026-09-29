#!/usr/bin/env python3
"""
`01_bench_data` — build the whole-transcriptome CM benchmark.

Walks all three participating CM datasets (`dtoxs_primary`, `dtoxs_ipsc`,
`gse244740`) from raw expression to a single benchmark table:

  raw counts -> log-CPM -> per-plate control baseline -> per-well Δ
             -> common gene universe (intersection, HGNC symbols)
             -> true significant DEG sets (moderated one-sample t, BH-FDR)
             -> cold-drug 5-fold + cold-dose leave-level-out split masks

Every artifact is written under `output/bench/`; gate assertions fail loudly if
any invariant is violated. No 978-landmark artifact is read anywhere.

Run (from anywhere):
    python 01_bench_data/build_bench_data.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
import gates  # noqa: E402
from baseline import per_plate_baselines, per_well_delta  # noqa: E402
from datasets import load_dataset  # noqa: E402
from genes import (align_matrix, collapse_to_symbols,  # noqa: E402
                   common_universe, load_ensembl_map)
from bench_io import (BENCH_DIR, group_key, save_bench, save_cold_dose,  # noqa: E402
                      save_cold_drug, save_dataset_summary, save_true_sets)
from normalize import log_cpm  # noqa: E402
from splits import cold_dose_masks, cold_drug_masks  # noqa: E402
from true_sets import group_true_sets  # noqa: E402
from utils import (DATASETS, DEG_FDR, DEG_ABS_MEAN_GATE,  # noqa: E402
                   DOSE_EPS, GENE_ID_SPACE, INPUT_DIR, MIN_WELLS_PER_GROUP)


def _ingest(dataset, ensembl_map):
    meta, counts, gene_ids = load_dataset(dataset)
    symbols, mat, stats = collapse_to_symbols(
        counts, gene_ids, GENE_ID_SPACE[dataset], ensembl_map)
    norm = log_cpm(mat)
    is_control = meta["is_control"].to_numpy()
    baselines, no_ctrl = per_plate_baselines(meta, norm, is_control)
    treated, delta, base_rows = per_well_delta(meta, norm, baselines, is_control)
    stats.update({"dataset": dataset, "n_wells_total": int(len(meta)),
                  "n_controls": int(is_control.sum()),
                  "n_treated": int(len(treated)),
                  "n_drugs": int(treated["drug"].nunique()),
                  "n_plates_without_controls": len(no_ctrl)})
    return {"meta": treated, "symbols": symbols, "delta": delta,
            "baseline": base_rows, "stats": stats}


def _write_split_manifest(fold_df, dose_df, groups):
    lines = ["# Split manifest", "",
             "## Cold-drug (primary, selection criterion)", ""]
    for f in sorted(fold_df.loc[fold_df["eligible"], "fold"].unique()):
        drugs = fold_df.loc[(fold_df["fold"] == f) & fold_df["eligible"], "drug"]
        lines.append(f"- fold {f}: {len(drugs)} held-out drugs "
                     f"(tested exactly once, zero train/test drug overlap)")
    lines += ["", f"- ineligible drugs (never tested): "
              f"{int((~fold_df['eligible']).sum())}", "",
              "## Cold-dose (secondary, robustness only)", ""]
    for _, r in dose_df[dose_df["held_out"]].iterrows():
        lvl = int(r["dose_level"])
        n = int(len(groups[(groups["dataset"] == "gse244740")
                           & (groups["dose_level"] == lvl)
                           & (groups["n_wells"] >= 2)]))
        lines.append(f"- {r['split']}: held-out dose level Con{lvl}, "
                     f"{n} scoreable gse244740 groups (>=2 wells for the t-test)")
    lines += ["",
              "**Caveat.** cold-dose true sets rest on ~3 wells per (drug, dose) "
              "in gse244740, so they are statistically weaker than the DToxS "
              "cold-drug true sets (>=4, often 8, wells per drug). The cold-dose "
              "axis therefore scores every group with >=2 wells rather than the "
              "4-well cold-drug bar, and is reported as a robustness column only; "
              "it never contaminates the cold-drug selection criterion.", ""]
    (BENCH_DIR / "splits" / "manifest.md").write_text("\n".join(lines) + "\n")


def _require_input_assets():
    missing = [str(p) for p in (INPUT_DIR / "ensembl_to_symbol.tsv",
                                INPUT_DIR / "gse244740_smiles.tsv")
               if not p.exists()]
    if missing:
        raise SystemExit(
            "missing frozen input assets:\n  " + "\n  ".join(missing) +
            "\nRegenerate them once (needs network):\n"
            "  python _lib/fetch_gene_symbol_map.py\n"
            "  python _lib/fetch_gse244740_smiles.py")


def main():
    _require_input_assets()
    ensembl_map = load_ensembl_map()
    per, symbol_sets, summary = {}, {}, []
    for ds in DATASETS:
        per[ds] = _ingest(ds, ensembl_map)
        symbol_sets[ds] = per[ds]["symbols"]
        summary.append(per[ds]["stats"])
        print(f"[ingest] {ds}: {per[ds]['stats']}")

    universe = common_universe(symbol_sets)
    print(f"[universe] {len(universe)} common genes")
    prev_path = BENCH_DIR / "universe.txt"
    if prev_path.exists():
        prev = prev_path.read_text().split()
        gates.gate(prev == universe,
                   "common gene universe changed on re-run (not byte-stable)")

    metas, bases, deltas = [], [], []
    for ds in DATASETS:
        p = per[ds]
        m = p["meta"].copy()
        m["log10_dose"] = np.log10(np.maximum(m["dose_uM"].to_numpy(), DOSE_EPS))
        m["group_key"] = [group_key(d, drug, lvl) for d, drug, lvl in
                          zip(m["dataset"], m["drug"], m["dose_level"])]
        metas.append(m)
        bases.append(align_matrix(p["symbols"], p["baseline"], universe))
        deltas.append(align_matrix(p["symbols"], p["delta"], universe))
    meta = pd.concat(metas, ignore_index=True)
    baseline = np.vstack(bases)
    delta = np.vstack(deltas)
    save_bench(meta, universe, baseline, delta)

    # ── ground-truth true sets per (dataset, drug, dose) ──────────────────
    rows, labels_list, order = [], [], []
    for key, idx in meta.groupby("group_key", sort=True).indices.items():
        g = meta.iloc[idx]
        lab, n_up, n_down = group_true_sets(
            delta[idx], fdr=DEG_FDR, abs_mean_gate=DEG_ABS_MEAN_GATE)
        rows.append({"group_key": key, "dataset": g["dataset"].iloc[0],
                     "drug": g["drug"].iloc[0], "dose_uM": float(g["dose_uM"].iloc[0]),
                     "dose_level": int(g["dose_level"].iloc[0]),
                     "n_wells": int(len(g)),
                     "reliable": bool(len(g) >= MIN_WELLS_PER_GROUP),
                     "n_sig_up": n_up, "n_sig_down": n_down})
        labels_list.append(lab)
        order.append(key)
    groups = pd.DataFrame(rows)
    labels = np.vstack(labels_list) if labels_list else np.zeros((0, len(universe)),
                                                                 dtype=np.int8)
    save_true_sets(groups, labels)

    # ── splits ───────────────────────────────────────────────────────────
    fold_df = cold_drug_masks(groups)
    dose_df = cold_dose_masks(groups)
    save_cold_drug(fold_df)
    save_cold_dose(dose_df)
    _write_split_manifest(fold_df, dose_df, groups)
    save_dataset_summary(pd.DataFrame(summary))

    # ── gates ────────────────────────────────────────────────────────────
    gates.check_universe(universe)
    gates.check_baselines_finite(baseline)
    gates.check_eligibility(groups, fold_df, MIN_WELLS_PER_GROUP)
    gates.check_cold_drug_masks(fold_df, groups)
    gates.check_cold_dose_masks(dose_df, groups)
    gates.check_metric()

    n_elig = int(fold_df["eligible"].sum())
    print(f"[groups] {len(groups)} (dataset,drug,dose) groups; "
          f"{int(groups['reliable'].sum())} reliable (>= {MIN_WELLS_PER_GROUP} wells)")
    print(f"[drugs] {fold_df['drug'].nunique()} union drugs; {n_elig} eligible")
    print(f"[samples] {len(meta)} treated wells over {len(universe)} genes")
    print(f"[out] {BENCH_DIR}")
    print("01_bench_data: OK")


if __name__ == "__main__":
    main()
