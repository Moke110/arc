#!/usr/bin/env python3
"""
`04_evaluate` — score entrants on cold-drug and cold-dose, write the leaderboard.

Loads the `01_bench_data` artifacts and runs each registered adapter (baselines
now; models as they land) through the `predict` seam over the persisted split
masks, then aggregates direction-aware AUPRC into macro (headline), micro and
per-drug-equal views and a written report.

Run (from anywhere):
    python 04_evaluate/evaluate.py
    python 04_evaluate/evaluate.py --adapters no-change mean-effect
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_lib"))
sys.path.insert(0, str(ROOT / "02_baselines"))
import evaluate as ev  # noqa: E402
from baselines import BASELINES  # noqa: E402
from bench_io import BENCH_DIR, load_bench, load_true_sets  # noqa: E402
from utils import OUTPUT_DIR  # noqa: E402

SPLITS_DIR = BENCH_DIR / "splits"


def build_adapters(names):
    from model_registry import MODEL_ADAPTERS  # model entrants (may be empty)

    registry = dict(BASELINES)
    registry.update(MODEL_ADAPTERS)
    if names:
        missing = [n for n in names if n not in registry]
        if missing:
            raise SystemExit(f"unknown adapters: {missing}; known: {sorted(registry)}")
        names = list(names)
    else:
        names = list(registry)
    return {n: registry[n]() for n in names}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adapters", nargs="*", default=None)
    args = ap.parse_args()

    bench = load_bench()
    groups, labels = load_true_sets()
    meta, baseline, delta = bench["meta"], bench["baseline"], bench["delta"]
    fold_df = pd.read_csv(SPLITS_DIR / "cold_drug.csv")
    dose_df = pd.read_csv(SPLITS_DIR / "cold_dose.csv")

    adapters = build_adapters(args.adapters)

    cd_rows, ce_rows = [], []
    per_group = []
    for name, adapter in adapters.items():
        gdf, pooled = ev.evaluate_cold_drug(adapter, meta, baseline, delta,
                                            groups, labels, fold_df)
        s = ev.summarize(gdf, pooled)
        cd_rows.append({"adapter": name, **s})
        per_group.append(gdf)
        gdf2, pooled2 = ev.evaluate_cold_dose(adapter, meta, baseline, delta,
                                              groups, labels, dose_df)
        if pooled2 is None:  # adapter cannot take dose -> N/A
            s2 = {"macro": None, "micro": None, "per_drug_equal": None,
                  "n_groups": 0, "n_datasets": 0, "n_drugs": 0}
        else:
            s2 = ev.summarize(gdf2, pooled2)
        ce_rows.append({"adapter": name, **s2})
        per_group.append(gdf2)
        print(f"[eval] {name}: cold-drug macro={s['macro'] if s['macro'] is None else round(s['macro'],4)}"
              f" cold-dose macro={None if s2['macro'] is None else round(s2['macro'],4)}")

    cd = pd.DataFrame(cd_rows).sort_values("macro", ascending=False,
                                           na_position="last").reset_index(drop=True)
    ce = pd.DataFrame(ce_rows).sort_values("macro", ascending=False,
                                           na_position="last").reset_index(drop=True)
    cd.insert(0, "rank", range(1, len(cd) + 1))
    ce.insert(0, "rank", range(1, len(ce) + 1))
    cd.to_csv(OUTPUT_DIR / "leaderboard_cold_drug.csv", index=False)
    ce.to_csv(OUTPUT_DIR / "leaderboard_cold_dose.csv", index=False)
    pd.concat(per_group, ignore_index=True).to_csv(
        OUTPUT_DIR / "per_group_auprc.csv", index=False)

    _write_report(cd, ce)
    print(f"[out] {OUTPUT_DIR/'leaderboard_cold_drug.csv'}")
    print(f"[out] {OUTPUT_DIR/'report.md'}")
    print("04_evaluate: OK")


def _sha(path):
    import hashlib

    p = Path(path)
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16] if p.exists() else "missing"


def _write_report(cd, ce):
    from model_registry import MODEL_ADAPTERS

    ce_map = {r["adapter"]: r for _, r in ce.iterrows()}
    # winner: highest cold-drug macro, ties broken by cold-dose macro
    def cd_macro(r):
        return -1.0 if pd.isna(r["macro"]) else float(r["macro"])

    def ce_macro(name):
        v = ce_map.get(name, {}).get("macro", None)
        return -1.0 if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)

    ranked = sorted(cd.to_dict("records"),
                    key=lambda r: (cd_macro(r), ce_macro(r["adapter"])), reverse=True)
    winner = ranked[0]["adapter"]

    lines = ["# Stage 2 v1 benchmark report", "",
             "Whole-transcriptome cardiomyocyte drug-perturbation model benchmark "
             "(ADR 0020). Metric: direction-aware AUPRC recovering each "
             "(dataset, drug, dose)'s true significant DEG set; headline = cold-drug "
             "macro AUPRC.", "",
             "## Leaderboard", "",
             "| rank | adapter | cold-drug macro | cold-drug micro | per-drug-equal "
             "| cold-dose macro | n_groups | n_drugs |",
             "|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(ranked, 1):
        cev = ce_map.get(r["adapter"], {})
        cev = cev.get("macro", None)
        cev = "N/A" if cev is None or (isinstance(cev, float) and np.isnan(cev)) \
            else f"{cev:.4f}"
        lines.append(
            f"| {i} | {r['adapter']} | {r['macro']:.4f} | {r['micro']:.4f} | "
            f"{r['per_drug_equal']:.4f} | {cev} | {r['n_groups']} | {r['n_drugs']} |")

    lines += ["", f"**Recommended winner: `{winner}`** — highest cold-drug macro "
              "AUPRC, ties broken by cold-dose macro AUPRC.", "",
              "## Entrants", "",
              f"- Scored adapters: {', '.join(cd['adapter'])}"]
    if MODEL_ADAPTERS:
        lines.append(f"- Model adapters present: {', '.join(MODEL_ADAPTERS)}")
    else:
        lines.append("- Model adapters present: none yet (baselines only); the "
                     "checkpoint audit (`docs/entrant_audit.md`) admits chemCPA and "
                     "State, whose adapters land in tickets 07-11.")

    lines += ["", "## Reproducibility manifest", "",
              "| artifact | sha256[:16] |", "|---|---|"]
    for rel in ["universe.txt", "matrices.npz", "groups.csv", "true_sets.npz",
                "splits/cold_drug.csv", "splits/cold_dose.csv"]:
        lines.append(f"| output/bench/{rel} | `{_sha(BENCH_DIR / rel)}` |")
    lines.append(f"| output/entrant_audit.csv | `{_sha(OUTPUT_DIR / 'entrant_audit.csv')}` |")

    from entrant_audit import AUDIT

    lines += ["", "Entrant environments (ADR 0021):", "",
              "| entrant | status | env | pinned commit |", "|---|---|---|---|",
              "| (baselines) | scored | `arc` | — |"]
    for r in AUDIT:
        env = f"`arc-stage2-v1-{r['model'].lower()}`" if r["status"] == "admitted" else "—"
        lines.append(f"| {r['model']} | {r['status']} | {env} | "
                     f"`{r['pinned_commit'][:10]}` |")

    lines += ["", "## Exit condition", "",
              "Met: a model leaderboard + a recommended winner + a reproducible "
              "benchmark. Stage 2 produces **no MIRI drug ranking**; the "
              "whole-transcriptome MIRI scorer (rebuilt 08/09) is deferred to a "
              "follow-on stage informed by the winner. No 978-landmark artifact is "
              "read anywhere.", ""]
    (OUTPUT_DIR / "report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
