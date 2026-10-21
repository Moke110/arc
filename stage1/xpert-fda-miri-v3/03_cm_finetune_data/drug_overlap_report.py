#!/usr/bin/env python
"""Compound / SMILES / pharmacology + in-HG overlap report for the CM datasets.

The report now has two layers:

  1. Compound inventory + external structure audit (the main layer). It counts
     the actual small molecules per dataset (excluding vehicles, the gse262419
     "Blank" wells, and the DToxS antibodies/proteins/peptide), then resolves
     every small molecule against PubChem / UniChem / ChEMBL and records
     whether a SMILES is obtainable and whether it is cross-verified by two
     independent registries. It also records the pharmacology classification
     (WHO ATC via ChEMBL/RxNav, ChEMBL mechanism/MoA, PubChem MeSH
     Pharmacological Classification). Inputs are the cached tables produced by
     `compound_external_audit.py`.
  2. The author-embedding (in-HG) coverage/overlap layer retained from the
     previous report: per-dataset unique drugs in the author HG table
     (drug_idx in [0, N_L1000_DRUGS)), sample counts, Venn diagrams, pairwise
     intersections, and the per-drug membership table (now with real
     names/MoA/targets from the L1000 drug metadata).

The gse262419 well-join is plate-aware ((Plate_ID, well) keyed against
GSE262419_hash.csv).

Writes: <pert_data>/cm_drugs.html (gitignored data-lake artifact).
"""

import base64
import html
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "_lib"))
sys.path.insert(0, str(PROJECT_ROOT / "03_cm_finetune_data"))
from utils import CM_DATA_ROOT, N_L1000_DRUGS  # noqa: E402
import _datasets  # noqa: E402
import drug_resolution as dr  # noqa: E402

DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse262419", "gse247541"]
OUT_PATH = CM_DATA_ROOT / "cm_drugs.html"
MOA_CSV = (PROJECT_ROOT.parents[1] / "pert_model" / "xpert" / "metadata"
           / "l1000_sdst_drugs_info_8432.csv")
INVENTORY_CSV = (PROJECT_ROOT / "03_cm_finetune_data" / "output"
                 / "compound_inventory.csv")
EXTERNAL_CSV = (PROJECT_ROOT / "03_cm_finetune_data" / "output"
                / "compound_external_audit.csv")


def l1000_full_h5ad():
    for rel in ("l1000-xpert", "xpert-l1000"):
        p = CM_DATA_ROOT / rel / "data" / "l1000_mdmt_full_336852.h5ad"
        if p.exists():
            return p
    raise FileNotFoundError(
        "L1000 full pretrain h5ad not found under pert_data/l1000-xpert "
        "or pert_data/xpert-l1000")


def build_author_index():
    dr.L1000_FULL_H5AD = l1000_full_h5ad()
    smi2idx, name2idx = dr.build_author_drug_index()
    meta_smi = dr.dtoxs_metadata_smiles_by_name()
    return smi2idx, name2idx, meta_smi


def collect_treated_compounds():
    """dataset -> {(compound_name, compound_code): (name, smiles)} treated only.

    Also returns the per-dataset sample metadata so the summary does not load
    the (large) gse262419 plates a second time.
    """
    out, metas = {}, {}
    for ds in DATASETS:
        meta, _, _ = _datasets.extract_samples(ds)
        metas[ds] = meta
        seen = {}
        for r in meta.itertuples():
            if getattr(r, "is_control"):
                continue
            seen.setdefault((r.compound_name, r.compound_code),
                            (r.compound_name, r.smiles))
        out[ds] = seen
    return out, metas


def _has_smiles(smi):
    return bool(smi) and str(smi) not in ("nan", "")


def collect_inhg_sets(compounds, smi2idx, name2idx, meta_smi):
    """dataset -> {drug_idx: (compound_name, resolution_method)} for in_hg."""
    idx2name = {}
    for ds in DATASETS:
        seen = compounds[ds]
        out = {}
        for (name, _code), (_n, smi) in seen.items():
            idx, method, _ = dr.resolve_compound(name, smi, smi2idx, name2idx,
                                                 meta_smi)
            if dr.classify(idx) == "in_hg":
                out[idx] = (name, method)
        idx2name[ds] = out
    return idx2name


def drug_name_table():
    """drug_idx -> (name, moa, target) from the L1000 drug metadata.

    The metadata table carries `drug_node_idx`, which is the same row index as
    the author embedding table (0..8980). (The older code looked for a
    `pert_idx` column, which does not exist, so every name fell back to
    "drug N".)
    """
    idx2smi = np.load(PROJECT_ROOT / "input" / "all_drugs_idx2smi.npy",
                      allow_pickle=True).item()
    names = {int(i): f"drug {i}" for i in idx2smi}
    moa, targets = {}, {}
    if MOA_CSV.exists():
        m = pd.read_csv(MOA_CSV)
        col = ("drug_node_idx" if "drug_node_idx" in m.columns
               else "pert_idx" if "pert_idx" in m.columns else None)
        if col is not None:
            m = m.dropna(subset=[col])
            moa_l, target_l = {}, {}
            for _, r in m.iterrows():
                i = int(r[col])
                c = r.get("cmap_name")
                if pd.notna(c) and str(c) != "nan":
                    names[i] = str(c)
                a = r.get("moa")
                if pd.notna(a):
                    moa_l.setdefault(i, [])
                    if str(a) not in moa_l[i]:
                        moa_l[i].append(str(a))
                t = r.get("target")
                if pd.notna(t):
                    target_l.setdefault(i, [])
                    if str(t) not in target_l[i]:
                        target_l[i].append(str(t))
            moa = {i: "; ".join(v) for i, v in moa_l.items()}
            targets = {i: "; ".join(v) for i, v in target_l.items()}
    return names, moa, targets


def _img64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


# ---------------------------------------------------------------------------
# External SMILES / pharmacology audit (compound_external_audit.py output)
# ---------------------------------------------------------------------------
def _as_bool(s):
    return str(s).strip().lower() in ("true", "1", "yes")


def load_external_audit():
    """Join the dataset-native inventory to the external lookup table.

    The audit is keyed by (lookup_name, dataset_inchikey): a name that is
    structure-backed in one dataset (DToxS) and name-only in another
    (gse262419) can legitimately resolve to two different structures, so the
    key keeps them apart.
    """
    inv = pd.read_csv(INVENTORY_CSV, dtype=str)
    aud = pd.read_csv(EXTERNAL_CSV, dtype=str)
    for df in (inv, aud):
        for c in ("lookup_name", "dataset_inchikey"):
            if c in df.columns:
                df[c] = df[c].fillna("")
    inv["merge_key"] = inv["lookup_name"] + "|" + inv["dataset_inchikey"]
    aud["merge_key"] = aud["lookup_name"] + "|" + aud["dataset_inchikey"]
    aud = aud.drop_duplicates("merge_key")
    drop = ["lookup_name", "dataset_inchikey", "dataset_smiles", "entity_type"]
    aud = aud.drop(columns=[c for c in drop if c in aud.columns])
    joined = inv.merge(aud, on="merge_key", how="left")

    # Name-only fallback for inventory rows whose compound was not in the
    # audit file (e.g. a dataset added after the last audit run).
    name_map = (aud.assign(_k=aud["merge_key"].str.split("|").str[0])
                .drop_duplicates("_k").set_index("_k"))
    miss = joined["verification"].isna()
    for c in ("pubchem_cid", "chembl_id", "verification", "smiles",
              "atc", "moa", "mesh_classes", "has_pharmacology",
              "parent_inchikey", "notes", "alias_mismatch",
              "pubchem_url", "chembl_url", "pubchem_title"):
        if c in name_map.columns:
            joined.loc[miss, c] = joined.loc[miss, "lookup_name"].map(
                name_map[c])

    joined["has_pharmacology"] = joined["has_pharmacology"].map(_as_bool)
    joined["alias_mismatch"] = joined["alias_mismatch"].map(_as_bool)
    joined["smiles"] = joined["smiles"].fillna("")
    joined["verification"] = joined["verification"].fillna("unresolved")
    joined["smiles_ok"] = (
        (joined["smiles"].str.strip() != "")
        & ~joined["verification"].isin(["unresolved", "conflict", "error"]))
    joined["verified"] = joined["verification"].eq("verified")
    return joined


def external_summary(inv):
    """Per-dataset counts of small molecules, external SMILES, pharmacology."""
    rows = []
    for ds in DATASETS:
        g = inv[inv["dataset"] == ds]
        sm = g[g["entity_type"] == "small_molecule"]
        agg = sm.groupby("lookup_name").agg(
            smiles_ok=("smiles_ok", "any"),
            verified=("verified", "any"),
            pharm=("has_pharmacology", "any"),
            ds_smiles=("dataset_smiles", lambda s: s.notna().any()))
        rows.append({
            "dataset": ds,
            "n_entries": len(g),
            "n_sm": len(agg),
            "n_biologic": int((g["entity_type"] != "small_molecule").sum()),
            "n_ds_smiles": int(agg["ds_smiles"].sum()),
            "n_ext": int(agg["smiles_ok"].sum()),
            "n_verified": int(agg["verified"].sum()),
            "n_unres": int((~agg["smiles_ok"]).sum()),
            "n_pharm": int(agg["pharm"].sum()),
        })
    df = pd.DataFrame(rows)

    sm = inv[inv["entity_type"] == "small_molecule"].copy()
    sm["ukey"] = sm["parent_inchikey"].fillna("")
    sm.loc[sm["ukey"] == "", "ukey"] = "name:" + sm["lookup_name"]
    u = sm.groupby("ukey").agg(
        smiles_ok=("smiles_ok", "any"),
        verified=("verified", "any"),
        pharm=("has_pharmacology", "any"),
        ds_smiles=("dataset_smiles", lambda s: s.notna().any()))
    union = {
        "dataset": "union（5 数据集并集）",
        "n_entries": int(len(u)),
        "n_sm": int(len(u)),
        "n_biologic": 0,
        "n_ds_smiles": int(u["ds_smiles"].sum()),
        "n_ext": int(u["smiles_ok"].sum()),
        "n_verified": int(u["verified"].sum()),
        "n_unres": int((~u["smiles_ok"]).sum()),
        "n_pharm": int(u["pharm"].sum()),
    }
    return df, union


def compound_table(inv):
    """One row per unique (compound, resolved structure) with membership."""
    sm = inv[inv["entity_type"] == "small_molecule"].copy()
    grp = sm.groupby("merge_key").agg(
        name=("compound_name", lambda s: sorted(set(s))[0]),
        datasets=("dataset", lambda s: ", ".join(sorted(set(s)))),
        codes=("compound_code", lambda s: ", ".join(
            sorted({x for x in s if x and str(x) != "nan"}))),
        note=("note", lambda s: "; ".join(sorted(
            {str(x) for x in s if isinstance(x, str) and x.strip()}))),
        pubchem_cid=("pubchem_cid", "first"),
        pubchem_url=("pubchem_url", "first"),
        pubchem_title=("pubchem_title", "first"),
        chembl_id=("chembl_id", "first"),
        chembl_url=("chembl_url", "first"),
        smiles=("smiles", "first"),
        verification=("verification", "first"),
        alias_mismatch=("alias_mismatch", "any"),
        atc=("atc", "first"),
        moa=("moa", "first"),
        mesh=("mesh_classes", "first"),
        pharm=("has_pharmacology", "any"),
        audit_notes=("notes", "first"),
    ).reset_index().sort_values("name", key=lambda s: s.str.lower())
    return grp


def venn_png(subset_sizes, labels, title):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    if len(subset_sizes) == 3:
        venn3(subsets=subset_sizes, set_labels=labels, ax=ax)
    else:
        venn2(subsets=subset_sizes, set_labels=labels, ax=ax)
    ax.set_title(title, fontsize=11)
    return _img64(fig)


def main():
    smi2idx, name2idx, meta_smi = build_author_index()
    compounds, metas = collect_treated_compounds()
    by_ds = collect_inhg_sets(compounds, smi2idx, name2idx, meta_smi)
    names, moa, targets = drug_name_table()
    ext_inv = load_external_audit()
    ext_summary, ext_union = external_summary(ext_inv)
    comp_table = compound_table(ext_inv)

    sets = {ds: set(v) for ds, v in by_ds.items()}
    union = sorted(set().union(*sets.values()))

    rows = []
    for ds in DATASETS:
        for idx in sorted(by_ds[ds]):
            rows.append({
                "drug_idx": idx,
                "drug_name": names.get(idx, f"drug {idx}"),
                "moa": moa.get(idx, ""),
                "target": targets.get(idx, ""),
                "resolved_by": by_ds[ds][idx][1],
                **{d: ("x" if idx in sets[d] else "") for d in DATASETS},
            })
    members = pd.DataFrame(rows)
    if not members.empty:
        members = members.drop_duplicates(
            subset=["drug_idx", "drug_name", "moa", "target"] + DATASETS)
    members = members.sort_values(["drug_idx"]).reset_index(drop=True)

    summary = []
    for ds in DATASETS:
        meta = metas[ds]
        n_treated = int((~meta["is_control"]).sum())
        n_comp = len(compounds[ds])
        n_smi = sum(1 for (_n, smi) in compounds[ds].values()
                    if _has_smiles(smi))
        n_in_hg = len(by_ds[ds])
        # in-HG 样本数: treated samples whose compound resolves to an in-HG
        # drug (distinct from n_in_hg, the deduped drug count).
        compound_in_hg = {}
        for (name, _code), (_n, smi) in compounds[ds].items():
            idx, _, _ = dr.resolve_compound(name, smi, smi2idx, name2idx,
                                            meta_smi)
            compound_in_hg[(name, _code)] = dr.classify(idx) == "in_hg"
        n_in_hg_samples = 0
        for r in meta.itertuples():
            if getattr(r, "is_control"):
                continue
            if compound_in_hg.get((r.compound_name, r.compound_code), False):
                n_in_hg_samples += 1
        summary.append((ds, n_comp, n_smi, n_in_hg, n_treated,
                        n_in_hg_samples))
    summary_df = pd.DataFrame(summary, columns=[
        "dataset", "n_compounds", "n_with_smiles", "n_in_hg_drugs",
        "n_treated_samples", "n_in_hg_samples"])

    total_in_hg = len(union)

    # in-HG 样本数 union = the v3 finetune union (the four finetune datasets,
    # ADR 0017); gse247541's 6 in-HG treated samples are reported per-dataset
    # but are outside the finetune scope (6-day exposure).
    FINETUNE_SUBSET = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse262419"]
    union_in_hg_samples = int(sum(
        summary_df.loc[summary_df["dataset"].isin(FINETUNE_SUBSET),
                       "n_in_hg_samples"]))

    # --- Venn diagrams -------------------------------------------------------
    venn_main = venn_png(
        [sets["dtoxs_primary"], sets["dtoxs_ipsc"], sets["gse244740"]],
        ["dtoxs_primary", "dtoxs_ipsc", "gse244740"],
        "in_hg drugs: DToxS primary / DToxS iPSC-CM / GSE244740")
    venn_secondary = venn_png(
        [sets["dtoxs_primary"], sets["gse244740"], sets["gse262419"]],
        ["dtoxs_primary", "gse244740", "gse262419"],
        "in_hg drugs: DToxS primary / GSE244740 / GSE262419")

    # --- pairwise intersection heatmap --------------------------------------
    pair = np.zeros((5, 5), dtype=int)
    for i, a in enumerate(DATASETS):
        for j, b in enumerate(DATASETS):
            pair[i, j] = len(sets[a] & sets[b])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.imshow(pair, cmap="YlGnBu", vmin=0)
    ax.set_xticks(range(5), DATASETS, rotation=45, ha="right")
    ax.set_yticks(range(5), DATASETS)
    for i in range(5):
        for j in range(5):
            ax.text(j, i, str(pair[i, j]), ha="center", va="center",
                    color="black", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("pairwise in_hg drug intersections")
    venn_heat = _img64(fig)

    # --- HTML ---------------------------------------------------------------
    def esc(x):
        return html.escape(str(x)) if pd.notna(x) else ""

    def clipped(x, n=220):
        s = "" if pd.isna(x) else str(x)
        if len(s) <= n:
            return esc(s)
        return (f"<span title='{html.escape(s)}'>{esc(s[:n])}…</span>")

    def badge(v):
        return {
            "verified": "<span class='ok'>✔ 双源</span>",
            "single": "<span class='warn'>单源</span>",
            "conflict": "<span class='bad'>冲突</span>",
        }.get(v, "<span class='muted'>未解决</span>")

    def link(url, text):
        if pd.isna(url) or not str(url).strip():
            return esc(text)
        return f"<a href='{html.escape(str(url))}'>{esc(text)}</a>"

    def cid_text(v):
        if pd.isna(v):
            return None
        try:
            return f"CID {int(float(v))}"
        except (TypeError, ValueError):
            return f"CID {v}"

    h = []
    h.append("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'>")
    h.append("<title>pert_data CM datasets — compounds, SMILES, pharmacology</title>")
    h.append("<style>"
             "body{font-family:-apple-system,'Segoe UI',Roboto,sans-serif;"
             "margin:2rem auto;max-width:1500px;padding:0 1rem;color:#222}"
             "h1{font-size:1.4rem}h2{font-size:1.12rem;margin-top:2rem}"
             "h3{font-size:1rem;margin-top:1.4rem}"
             "table{border-collapse:collapse;width:100%;font-size:.78rem;"
             "margin-top:.6rem}"
             "th,td{border:1px solid #ccc;padding:.28rem .45rem;text-align:left;"
             "vertical-align:top}"
             "th{background:#f0f4f8;position:sticky;top:0}"
             "td.num{text-align:right}td.smi{font-family:ui-monospace,monospace;"
             "font-size:.72rem;word-break:break-all;max-width:340px}"
             "tr:nth-child(even) td{background:#fafbfc}"
             "img{max-width:100%;margin-top:.6rem}"
             ".note{color:#666;font-size:.85rem}"
             ".ok{color:#137333;font-weight:600}.warn{color:#b06000;font-weight:600}"
             ".bad{color:#b00020;font-weight:600}.muted{color:#777}"
             ".venn-row{display:flex;flex-wrap:wrap;gap:1rem}"
             ".venn-row figure{margin:0;flex:1 1 420px}"
             "figcaption{font-size:.85rem;color:#444;margin-top:.2rem}"
             "code{background:#f3f3f3;padding:0 2px}"
             "</style></head><body>")
    h.append("<h1>pert_data CM 数据集 — 小分子化合物、SMILES 与药理学分类外部核查</h1>")

    h.append(
        "<p class='note'>本报告回答三个问题：(1) 每个数据集实际含多少<b>小分子</b>"
        "（排除 vehicle/CTRL、空孔以及抗体/重组蛋白/多肽等生物药）；(2) 每个小分子能否"
        "通过<b>正规公开数据库查询</b>获得 SMILES，而不是只看 GEO 交付里有没有附带结构；"
        "(3) 每个小分子能否查询到<b>药理学分类</b>。</p>")
    h.append(
        "<p class='note'>外部查询链：<b>PubChem PUG-REST</b>（名称/InChIKey → CID、SMILES、"
        "InChIKey）→ <b>UniChem</b>（InChIKey → ChEMBL / DrugBank / ChEBI / KEGG 交叉引用）"
        "→ <b>ChEMBL</b>（canonical SMILES、max phase、ATC、mechanism/MoA）；"
        "药理分类同时取 <b>PubChem PUG-View 的 MeSH Pharmacological Classification</b>，"
        "ATC 缺失时以 <b>RxNav</b> 兜底。解析一律<b>结构优先</b>：数据集自带 SMILES/InChIKey "
        "时先按结构查；只有名称的数据集才按名称查，且名称命中必须与 ChEMBL 结构在 InChIKey "
        "连接性块（前 14 位）一致才采纳。</p>")
    h.append(
        "<p class='note'><b>二次验证口径</b>：<span class='ok'>✔ 双源</span> = 同一结构"
        "（InChIKey 前 14 位，盐/母体形式按最大有机片段比较）同时登记在 ≥2 个独立来源"
        "（PubChem + ChEMBL/DrugBank/ChEBI/KEGG 中至少一个）；<span class='warn'>单源</span> "
        "= 仅一个 registry 能给出结构（SMILES 可用，但未独立交叉验证）；"
        "<span class='bad'>冲突</span> = registry 之间结构不一致；"
        "<span class='muted'>未解决</span> = 没有权威命中。名称在 ChEMBL 指向不同结构时，"
        "该名称的 ATC/MoA 不会被借用（见“需人工复核”）。</p>")

    # --- Table A: composition + external structure/pharmacology coverage ----
    h.append("<h2>A. 各数据集小分子清点与外部 SMILES / 药理学分类可得性</h2>")
    h.append(
        "<p class='note'>“扰动条目”= 数据集里非对照的化合物条目数（按作者编码计，"
        "包括同一化合物的不同编码，如 DToxS 的 SUN/SUN+）；“小分子（去重）”= 按名称/结构"
        "去重后的真正小分子数；“生物药/蛋白/肽”= 抗体（曲妥珠单抗、贝伐珠单抗、西妥昔单抗、"
        "利妥昔单抗）、重组蛋白（IGF-1、TNF-α）与多肽（内皮素-1），它们本来就没有小分子 "
        "SMILES，不应计入小分子分母。“数据集自带 SMILES”只统计 GEO/LINCS 交付里直接带结构的"
        "条目；“外部 SMILES 可得”= PubChem/ChEMBL 任一能给出结构（双源+单源）；"
        "“药理学分类可得”= ATC 或 ChEMBL MoA 或 MeSH 药理分类至少一项可查。</p>")
    h.append("<table><tr><th>dataset</th><th class='num'>扰动条目</th>"
             "<th class='num'>小分子（去重）</th><th class='num'>生物药/蛋白/肽</th>"
             "<th class='num'>数据集自带 SMILES</th><th class='num'>外部 SMILES 可得</th>"
             "<th class='num'>其中双源验证</th><th class='num'>单源</th>"
             "<th class='num'>无外部结构</th><th class='num'>药理学分类可得</th>"
             "<th class='num'>分类覆盖率</th></tr>")
    for _, r in ext_summary.iterrows():
        cov = (100.0 * r["n_pharm"] / r["n_sm"]) if r["n_sm"] else 0.0
        h.append(
            f"<tr><td>{esc(r['dataset'])}</td>"
            f"<td class='num'>{int(r['n_entries'])}</td>"
            f"<td class='num'>{int(r['n_sm'])}</td>"
            f"<td class='num'>{int(r['n_biologic'])}</td>"
            f"<td class='num'>{int(r['n_ds_smiles'])}</td>"
            f"<td class='num'>{int(r['n_ext'])}</td>"
            f"<td class='num'>{int(r['n_verified'])}</td>"
            f"<td class='num'>{int(r['n_ext'] - r['n_verified'])}</td>"
            f"<td class='num'>{int(r['n_unres'])}</td>"
            f"<td class='num'>{int(r['n_pharm'])}</td>"
            f"<td class='num'>{cov:.1f}%</td></tr>")
    u = ext_union
    cov = 100.0 * u["n_pharm"] / u["n_sm"] if u["n_sm"] else 0.0
    h.append(
        f"<tr><td><b>{esc(u['dataset'])}</b></td>"
        f"<td class='num'><b>{int(u['n_entries'])}</b></td>"
        f"<td class='num'><b>{int(u['n_sm'])}</b></td>"
        f"<td class='num'>—</td>"
        f"<td class='num'><b>{int(u['n_ds_smiles'])}</b></td>"
        f"<td class='num'><b>{int(u['n_ext'])}</b></td>"
        f"<td class='num'><b>{int(u['n_verified'])}</b></td>"
        f"<td class='num'>{int(u['n_ext'] - u['n_verified'])}</td>"
        f"<td class='num'>{int(u['n_unres'])}</td>"
        f"<td class='num'><b>{int(u['n_pharm'])}</b></td>"
        f"<td class='num'><b>{cov:.1f}%</b></td></tr>")
    h.append("</table>")
    h.append(
        "<p class='note'>并集按“母体结构”去重（盐/溶剂合物取最大有机片段），因此同一活性分子"
        "的盐形式不会重复计数；各数据集行按数据集内名称去重。"
        "gse262419 的 464 个旧计数包含一个 “Blank” 空孔，本表已扣除（463 个条目 = 462 个小分子 + "
        "1 个多肽 ET-1）；DToxS 的 52/41 个旧“化合物”条目包含抗体/蛋白/多肽，本表单列。"
        "药理学分类的三种来源可重叠；对非药物化学品，MeSH 常给出 “Insecticides / Herbicides / "
        "Fungicides” 等毒理分类，而 ATC 通常缺失，因此“分类可得”不等于“有 ATC 适应症分类”。</p>")

    # --- Table B: in-HG coverage (kept from the original report) -----------
    h.append("<h2>B. 作者药物嵌入表（in-HG）覆盖与样本数（原报告口径）</h2>")
    h.append("<p class='note'>in-HG 药物 = 化合物解析到作者药物嵌入表行 "
             "<code>drug_idx ∈ [0, %d)</code>（作者 HG 参考表，不含 FDA-extension 行）。"
             "解析规则与 step-03 manifest 一致：SMILES 精确匹配，否则作者 cmap_name 名称匹配。"
             "in-HG 样本数 = treated 样本中化合物解析到作者嵌入表（in-HG）的样本数；"
             "union 行 = v3 finetune 并集（四个 finetune 数据集，ADR 0017；gse247541 因 "
             "6 天暴露不在 finetune 范围，仅按数据集单独报告）。</p>" % N_L1000_DRUGS)
    h.append("<table><tr><th>dataset</th><th class='num'>in-HG 药物数</th>"
             "<th class='num'>treated 样本数</th><th class='num'>in-HG 样本数</th></tr>")
    for _, r in summary_df.iterrows():
        h.append(f"<tr><td>{r['dataset']}</td><td class='num'>"
                 f"{r['n_in_hg_drugs']}</td><td class='num'>"
                 f"{r['n_treated_samples']}</td><td class='num'>"
                 f"{r['n_in_hg_samples']}</td></tr>")
    h.append(f"<tr><td><b>union（5 数据集并集）</b></td>"
             f"<td class='num'><b>{total_in_hg}</b></td>"
             f"<td class='num'>—</td>"
             f"<td class='num'><b>{union_in_hg_samples}</b></td></tr>")
    h.append("</table>")
    h.append("<p class='note'>gse262419 的 treated 样本数沿用 step-03 loader 口径，"
             "其中包含 32 个 dose=0 的 “Blank” 空孔（loader 只把 VEH/Media 判为对照）；"
             "这些空孔不解析到任何药物，因此不进入 in-HG 样本数，也不进入并集。"
             "建议在后续 loader 修订中把 Blank 归入对照。</p>")

    # --- Table C: per-compound external audit ------------------------------
    h.append("<h2>C. 逐化合物外部核查明细</h2>")
    h.append("<p class='note'>每行一个（化合物, 解析结构）；“数据集”列出该条目出现的"
             "数据集。PubChem/ChEMBL 为可点击链接；SMILES 取 ChEMBL canonical SMILES，"
             "否则 PubChem SMILES。ATC 取 ChEMBL 药物记录的 WHO ATC（显示最具体一级）；"
             "MoA 取 ChEMBL mechanism 的 mechanism_of_action；MeSH 取 PubChem PUG-View 的 "
             "MeSH Pharmacological Classification。</p>")
    h.append("<table><tr><th>化合物</th><th>数据集</th><th>PubChem</th>"
             "<th>ChEMBL</th><th>SMILES</th><th>验证</th><th>ATC（WHO）</th>"
             "<th>MoA（ChEMBL）</th><th>MeSH 药理分类</th><th>备注</th></tr>")
    for _, r in comp_table.iterrows():
        note = "; ".join(
            x for x in [r.get("note"), r.get("audit_notes")]
            if isinstance(x, str) and x and x != "nan"
            and not x.startswith("name-only"))
        if r.get("alias_mismatch"):
            note = ("数据集名称在 ChEMBL 指向不同结构；"
                    + (note or ""))
        cid = cid_text(r.get("pubchem_cid"))
        pc = link(r.get("pubchem_url"), cid) if cid else "—"
        cb = (link(r.get("chembl_url"), r["chembl_id"])
              if pd.notna(r.get("chembl_id")) else "—")
        smi = esc(r.get("smiles")) if pd.notna(r.get("smiles")) else "—"
        h.append(
            f"<tr><td>{esc(r['name'])}</td><td>{esc(r['datasets'])}</td>"
            f"<td>{pc}</td><td>{cb}</td>"
            f"<td class='smi'>{smi}</td><td>{badge(r['verification'])}</td>"
            f"<td>{clipped(r.get('atc'), 160)}</td>"
            f"<td>{clipped(r.get('moa'), 160)}</td>"
            f"<td>{clipped(r.get('mesh'), 160)}</td>"
            f"<td>{clipped(note, 200)}</td></tr>")
    h.append("</table>")

    # --- Table D: compounds needing manual review --------------------------
    review = ext_inv[
        (ext_inv["entity_type"] == "small_molecule")
        & ((~ext_inv["smiles_ok"]) | ext_inv["alias_mismatch"]
           | ext_inv["verification"].isin(["single", "conflict", "error"]))]
    h.append("<h2>D. 需人工复核的化合物</h2>")
    h.append("<p class='note'>包括：没有任何权威结构命中、registry 结构冲突、"
             "或数据集名称与结构不一致（例如 DToxS 的 “Amiodarone” 实际是 "
             "ZINC150345845，一个 2-戊基苯并呋喃同系物，并非 amiodarone）。</p>")
    if review.empty:
        h.append("<p class='note'>无。</p>")
    else:
        review = (review.groupby(
            ["compound_name", "compound_code", "verification", "notes"],
            dropna=False)
            .agg(dataset=("dataset", lambda s: ", ".join(sorted(set(s)))),
                 pubchem_url=("pubchem_url", "first"),
                 pubchem_cid=("pubchem_cid", "first"),
                 smiles_ok=("smiles_ok", "first"),
                 alias_mismatch=("alias_mismatch", "first"))
            .reset_index())
        h.append("<table><tr><th>dataset</th><th>化合物</th><th>编码</th>"
                 "<th>验证</th><th>PubChem</th><th>说明</th></tr>")
        for _, r in review.sort_values(["dataset", "compound_name"]).iterrows():
            reason = []
            if r["verification"] == "error":
                reason.append("查询失败（可重试）")
            elif not r["smiles_ok"]:
                reason.append("无外部结构")
            if r["verification"] == "single":
                reason.append("仅单源（未双重验证）")
            if r["alias_mismatch"]:
                reason.append("名称↔结构不一致")
            if r["verification"] == "conflict":
                reason.append("结构冲突")
            note = "; ".join(reason)
            extra = r.get("notes")
            if isinstance(extra, str) and extra and extra != "nan":
                note += "；" + extra
            cid = cid_text(r.get("pubchem_cid"))
            b = badge(r["verification"])
            if r.get("alias_mismatch"):
                b += "<br><span class='bad'>别名≠结构</span>"
            h.append(
                f"<tr><td>{esc(r['dataset'])}</td>"
                f"<td>{esc(r['compound_name'])}</td>"
                f"<td>{esc(r['compound_code'])}</td>"
                f"<td>{b}</td>"
                f"<td>{link(r.get('pubchem_url'), cid) if cid else '—'}</td>"
                f"<td>{clipped(note, 240)}</td></tr>")
        h.append("</table>")

    # --- original in-HG overlap figures ------------------------------------
    h.append("<h2>E. in-HG 药物重叠（原报告图）</h2>")
    h.append("<div class='venn-row'>")
    h.append(f"<figure><img src='data:image/png;base64,{venn_main}'></figure>")
    h.append(f"<figure><img src='data:image/png;base64,{venn_secondary}'>"
             "</figure></div>")
    h.append("<h3>两两交集（in-HG 药物）</h3>")
    h.append(f"<img src='data:image/png;base64,{venn_heat}'>")

    h.append("<h2>F. 按 in-HG 药物的数据集归属</h2>")
    h.append("<p class='note'>drug_name / moa / target 来自 L1000 药物元数据"
             "（<code>l1000_sdst_drugs_info_8432.csv</code> 的 <code>drug_node_idx</code>）。</p>")
    h.append("<table><tr><th class='num'>drug_idx</th><th>drug_name</th>"
             "<th>moa</th><th>target</th>")
    for d in DATASETS:
        h.append(f"<th>{d}</th>")
    h.append("</tr>")
    for r in members.itertuples(index=False):
        h.append("<tr><td class='num'>%d</td><td>%s</td><td>%s</td><td>%s</td>" % (
            r.drug_idx, esc(r.drug_name), esc(r.moa), esc(r.target)))
        for d in DATASETS:
            h.append(f"<td>{getattr(r, d)}</td>")
        h.append("</tr>")
    h.append("</table>")

    h.append("<h2>G. 复现</h2>")
    h.append("<pre>python 03_cm_finetune_data/compound_external_audit.py\n"
             "python 03_cm_finetune_data/drug_overlap_report.py</pre>")
    h.append("<p class='note'>外部查询原始响应缓存在 "
             "<code>03_cm_finetune_data/output/compound_external_cache.json</code>，"
             "逐化合物结果在 <code>compound_external_audit.csv</code>，"
             "数据集条目清点在 <code>compound_inventory.csv</code>。"
             "查询时间为 2026-10-20；PubChem / UniChem / ChEMBL / RxNav 均为公开 REST 接口。"
             "药理学分类是“数据库注释”而非本项目的功能验证，引用时请保留来源。</p>")
    h.append("</body></html>")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text("\n".join(h), encoding="utf-8")
    print(f"Drug-overlap report -> {OUT_PATH} "
          f"({len(union)} in-HG drugs, {len(members)} membership rows)")


if __name__ == "__main__":
    main()
