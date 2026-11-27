#!/usr/bin/env python3
"""
_other — proliferation-pathway effect of a v2-screened drug on the
workflow's predicted cardiomyocyte (CM) perturbed transcriptome.

Question asked for this directory: for a drug screened by the v2 workflow
(07_xpert_predict), score every proliferation-related GO pathway for the
direction/strength of the drug's predicted effect on it.

Inputs
------
  * 07_xpert_predict/output/xpert_fda-cm_deg-zs.csv  -- 1,004 drugs x 978 L1000
    landmarks, z-scored per gene across drugs (the workflow's recommended input;
    step 07). Each row is the XPert-predicted DEG (dose_idx 0 / time_idx 2,
    per-plate-normalised union-CM finetune) after gene-level z-scoring.
  * output/cache/go_prolif_terms.json       -- GO proliferation vocabulary.
  * output/cache/go_prolif_annotations.csv  -- human annotations under the
    cell-population-proliferation subtree (QuickGO, cached by go_prolif_fetch.py).

Gene-set scoring
----------------
  gene set of GO term T = human protein annotations with goId == T mapped to the
  978 landmarks by gene symbol; NOT-qualified annotations are dropped. Terms are
  scored on their L1000-landmark overlap (min --min-genes).

  pathway score = mean of the drug's predicted per-gene z over the pathway's
  landmark genes (how far the drug moves the pathway as a block vs the
  1,004-drug mean). A percentile is also reported: rank of the drug's pathway
  mean-z inside the empirical distribution over all 1,004 screened drugs
  (per-drug same pathway mean), so "top decile" etc. is meaningful. A GSEA-style
  running-sum enrichment score (same stat as _lib/scoring.enrichment_score) is
  given for comparability with the pipeline's own pathway metric.

  Terms are tagged into families by GO semantics of the term name:
    process    -- proliferation process itself (promotes proliferation when UP)
    pos_reg    -- positive regulation of <..> proliferation
    neg_reg    -- negative regulation of <..> proliferation (pro-proliferation
                  when DOWN)
    reg        -- generic regulation of <..> proliferation
  For pos_reg/neg_reg/process an explicit polarity (+1/-1/+1) multiplies the
  mean z into a directional effect d (positive = drug pushes the pathway toward
  more proliferation).

Outputs
-------
  output/<drug>_prolif_pathway_scores.csv   per-term table
  output/<drug>_prolif_family.csv           family-level aggregates
  output/<drug>_prolif_landmark_genes.csv   long-form overlap gene sets
  output/<drug>_prolif_report.md            bilingual analysis summary

Usage:
    python _other/mirtazapine_prolif_score.py --drug MIRTAZAPINE
    python _other/mirtazapine_prolif_score.py --drug MIRTAZAPINE --min-genes 5
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

WORKFLOW = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(WORKFLOW / "_lib"))
from scoring import enrichment_score  # noqa: E402  (same ES stat as the pipeline)

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
DEG_ZS = WORKFLOW / "07_xpert_predict" / "output" / "xpert_fda-cm_deg-zs.csv"
TERMS = OUT / "cache" / "go_prolif_terms.json"
ANN = OUT / "cache" / "go_prolif_annotations.csv"


def category_of(name):
    """GO family tag by the term name's regulation semantics."""
    low = name.lower()
    if low.startswith("positive regulation"):
        return "pos_reg"
    if low.startswith("negative regulation"):
        return "neg_reg"
    if "regulation" in low:
        return "reg"
    return "process"


def polarity_of(cat):
    return {"process": 1.0, "pos_reg": 1.0, "neg_reg": -1.0, "reg": np.nan}[cat]


def load_zs(drug):
    df = pd.read_csv(DEG_ZS)
    gene_cols = df.columns[3:].tolist()
    if len(gene_cols) != 978:
        print(f"warning: {len(gene_cols)} gene columns (expected 978)", file=sys.stderr)
    hit = df[df["drug_name"].str.upper() == drug.upper()]
    if hit.empty:
        avail = ", ".join(sorted(df["drug_name"].unique()[:20]))
        raise SystemExit(f"drug '{drug}' not in v2 step-07 table; first drugs: {avail}")
    row = hit.iloc[0]
    z = row[gene_cols].astype(float).values
    return row, gene_cols, z, df, gene_cols


def load_gene_sets():
    terms = json.loads(TERMS.read_text())
    ann = pd.read_csv(ANN, dtype=str)
    ann = ann[ann["qualifier"].fillna("") != "NOT"]          # drop NOT rows
    sets = defaultdict(dict)                                  # go_id -> {symbol}
    for rec in ann.itertuples(index=False):
        go = rec.goId
        sym = rec.symbol
        if go and sym:
            sets[go][sym] = None
    return terms, sets, ann


def family_membership(terms):
    fam = {}
    for g, name in terms.items():
        fam[g] = category_of(name)
    return fam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="MIRTAZAPINE")
    ap.add_argument("--min-genes", type=int, default=6,
                    help="min landmark overlap to score a term (pipeline uses 6)")
    args = ap.parse_args()

    row, gene_cols, z, df_all, _ = load_zs(args.drug)
    drugs = df_all["drug_name"].astype(str).str.upper()
    zmat = df_all[gene_cols].astype(float).values            # drugs x genes

    terms, sets, ann = load_gene_sets()
    fam = family_membership(terms)
    landmarks = set(gene_cols)

    # ---- build per-term records -------------------------------------------
    records = []
    coverage = []
    n_term_skipped_lowcov = 0
    for go, name in terms.items():
        syms_all = sorted(set(sets.get(go, {})) & landmarks)
        cov = {
            "go_id": go, "term": name, "family": fam.get(go, "process"),
            "n_annot_symbols": len(sets.get(go, {})),
            "n_l1000_landmark_overlap": len(syms_all),
        }
        if go not in sets or len(syms_all) < args.min_genes:
            cov["scored"] = False
            coverage.append(cov)
            n_term_skipped_lowcov += 1
            continue
        idx = np.array([gene_cols.index(s) for s in syms_all])
        zs = z[idx]
        mean_z = float(zs.mean())
        drug_means = zmat[:, idx].mean(axis=1)
        pct = float((np.mean(drug_means <= mean_z)) * 100)
        mask = np.zeros(len(gene_cols), bool)
        mask[idx] = True
        es = float(enrichment_score(z, mask, weight=1.0))
        cat = fam.get(go, "process")
        pol = polarity_of(cat)
        cov["scored"] = True
        coverage.append(cov)
        records.append({
            "go_id": go, "term": name, "family": cat, "polarity": pol,
            "n_annot_symbols": len(sets[go]),
            "n_landmark_genes": len(syms_all),
            "drug": args.drug.upper(),
            "mean_z": mean_z,
            "directional_effect": (pol * mean_z if not np.isnan(pol) else np.nan),
            "gsea_ES": es,
            "frac_up": float((zs > 0).mean()),
            "drug_pathway_mean_z_pctile": pct,
            "landmark_genes": ";".join(syms_all),
        })
    if not records:
        raise SystemExit(f"no proliferation terms with >= {args.min_genes} landmarks")
    t = pd.DataFrame(records)

    # reference module gene sets (the 5 canonical proliferation terms)
    ref = [g for g in terms if g in ("GO:0008283", "GO:0042127",
                                     "GO:0008284", "GO:0008285", "GO:0033002")]
    ref_df = t[t["go_id"].isin(ref)][["go_id", "term", "family", "n_landmark_genes"]]

    # ---- family-level aggregates (union of the family's term gene sets) ----
    fam_union = defaultdict(list)
    for go in terms:
        if go not in sets:
            continue
        syms = [s for s in sets[go] if s in landmarks]
        if syms:
            fam_union[fam.get(go, "process")].extend(syms)
    fam_rows = []
    for cat in ["process", "pos_reg", "neg_reg", "reg"]:
        syms = sorted(set(fam_union[cat]))
        if not syms:
            continue
        idx = np.array([gene_cols.index(s) for s in syms])
        zs = z[idx]
        drug_means = zmat[:, idx].mean(axis=1)
        mean_z = float(zs.mean())
        fam_rows.append({
            "family": cat,
            "n_terms": int(sum(fam.get(g) == cat for g in terms)),
            "n_landmark_genes": len(syms),
            "mean_z": mean_z,
            "frac_up": float((zs > 0).mean()),
            "drug_pathway_mean_z_pctile": float((drug_means <= mean_z).mean() * 100),
            "landmark_genes": ";".join(syms),
        })
    fam_df = pd.DataFrame(fam_rows)

    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / f"{args.drug.lower()}_prolif_pathway_scores.csv", index=False)
    fam_df.to_csv(OUT / f"{args.drug.lower()}_prolif_family.csv", index=False)
    pd.DataFrame(coverage).to_csv(
        OUT / f"{args.drug.lower()}_prolif_allterms_coverage.csv", index=False)

    order = t.sort_values("mean_z")
    up = t[t["mean_z"] > 0].sort_values("mean_z", ascending=False)
    down = order[order["mean_z"] < 0].head(10)

    print(f"drug {args.drug.upper()} (drug_idx {int(row['drug_idx'])})")
    print(f"proliferation terms scored: {len(t)} "
          f"(skipped low-coverage: {n_term_skipped_lowcov})")
    print(f"  strongest DOWN (suppressed) block mean-z terms:")
    for _, r in down.iterrows():
        print(f"    {r['go_id']} {r['term'][:60]:60s} n={r['n_landmark_genes']:3d} "
              f"mean_z={r['mean_z']:+.3f} pctile={r['drug_pathway_mean_z_pctile']:5.1f}")
    print(f"  strongest UP (activated) block mean-z terms:")
    for _, r in up.head(10).iterrows():
        print(f"    {r['go_id']} {r['term'][:60]:60s} n={r['n_landmark_genes']:3d} "
              f"mean_z={r['mean_z']:+.3f} pctile={r['drug_pathway_mean_z_pctile']:5.1f}")
    print("family aggregates:")
    print(fam_df[["family", "n_terms", "n_landmark_genes", "mean_z", "frac_up"]].to_string(index=False))
    print("canonical reference terms in table:")
    print(ref_df.to_string(index=False) if not ref_df.empty else "(none with >= min genes)")


if __name__ == "__main__":
    main()
