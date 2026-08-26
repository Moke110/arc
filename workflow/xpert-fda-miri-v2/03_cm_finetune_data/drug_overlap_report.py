#!/usr/bin/env python
"""Drug-coverage + overlap report across the pert_data CM datasets.

Counts, per dataset under pert_data (excluding the L1000 pretrain corpus),
how many unique treated compounds are present in the author-provided drug
embedding table (the author HG reference, drug_idx in [0, N_L1000_DRUGS)),
and how those drug sets overlap with each other. Emits a self-contained
HTML report (Venn diagrams + pairwise intersection heatmap + per-drug
membership table).

Writes: <pert_data>/cm_drugs.html (gitignored data-lake artifact).
"""

import base64
import io
import itertools
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


def collect_inhg_sets(smi2idx, name2idx, meta_smi):
    """dataset -> {drug_idx: (compound_name, resolution_method)} for in_hg."""
    idx2name = {}
    for ds in DATASETS:
        meta, _, _ = _datasets.extract_samples(ds)
        seen = {}
        for r in meta.itertuples():
            if getattr(r, "is_control"):
                continue
            seen.setdefault((r.compound_name, r.compound_code),
                            (r.compound_name, r.smiles))
        out = {}
        for (name, _code), (_n, smi) in seen.items():
            idx, method, _ = dr.resolve_compound(name, smi, smi2idx, name2idx,
                                                 meta_smi)
            if dr.classify(idx) == "in_hg":
                out[idx] = (name, method)
        idx2name[ds] = out
    return idx2name


def drug_name_table():
    """drug_idx -> (name, moa) from author metadata (best-effort)."""
    idx2smi = np.load(PROJECT_ROOT / "input" / "all_drugs_idx2smi.npy",
                      allow_pickle=True).item()
    names = {int(i): f"drug {i}" for i in idx2smi}
    try:
        obs = sc.read_h5ad(l1000_full_h5ad(), backed="r").obs
        cmap = obs[["cmap_name", "pert_idx"]].drop_duplicates("pert_idx")
        for c, i in zip(cmap["cmap_name"], cmap["pert_idx"]):
            if str(c) != "nan":
                names.setdefault(int(i), str(c))
    except Exception:
        pass
    moa = {}
    if MOA_CSV.exists():
        m = pd.read_csv(MOA_CSV)
        if "pert_idx" in m.columns and "moa" in m.columns:
            for i, a in zip(m["pert_idx"], m["moa"]):
                if pd.notna(a):
                    moa.setdefault(int(i), str(a))
    return names, moa


def _img64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


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
    by_ds = collect_inhg_sets(smi2idx, name2idx, meta_smi)
    names, moa = drug_name_table()

    sets = {ds: set(v) for ds, v in by_ds.items()}
    union = sorted(set().union(*sets.values()))
    total_in_hg = len(union)

    rows = []
    for ds in DATASETS:
        for idx in sorted(by_ds[ds]):
            rows.append({
                "drug_idx": idx,
                "drug_name": names.get(idx, f"drug {idx}"),
                "moa": moa.get(idx, ""),
                "resolved_by": by_ds[ds][idx][1],
                **{d: ("x" if idx in sets[d] else "") for d in DATASETS},
            })
    members = pd.DataFrame(rows).sort_values(
        ["drug_idx"]).reset_index(drop=True)

    summary = []
    for ds in DATASETS:
        meta, _, _ = _datasets.extract_samples(ds)
        n_treated = int((~meta["is_control"]).sum())
        n_unique = len(by_ds[ds])
        summary.append((ds, n_treated, n_unique))
    summary_df = pd.DataFrame(
        summary, columns=["dataset", "n_treated_samples", "n_in_hg_drugs"])

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
    h = []
    h.append("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'>")
    h.append("<title>pert_data CM datasets — in-HG drug overlap</title>")
    h.append("<style>"
             "body{font-family:-apple-system,'Segoe UI',Roboto,sans-serif;"
             "margin:2rem auto;max-width:1100px;padding:0 1rem;color:#222}"
             "h1{font-size:1.4rem}h2{font-size:1.15rem;margin-top:2rem}"
             "table{border-collapse:collapse;width:100%;font-size:.85rem;"
             "margin-top:.6rem}"
             "th,td{border:1px solid #ccc;padding:.3rem .5rem;text-align:left}"
             "th{background:#f0f4f8}td.num{text-align:right}"
             "img{max-width:100%;margin-top:.6rem}"
             ".note{color:#666;font-size:.85rem}"
             ".venn-row{display:flex;flex-wrap:wrap;gap:1rem}"
             ".venn-row figure{margin:0;flex:1 1 420px}"
             "figcaption{font-size:.85rem;color:#444;margin-top:.2rem}"
             "</style></head><body>")
    h.append("<h1>pert_data CM 数据集 — 作者药物嵌入表（in-HG）覆盖与重叠</h1>")
    h.append("<p class='note'>in-HG 药物 = 化合物解析到作者药物嵌入表行 "
             "<code>drug_idx ∈ [0, %d)</code>（作者 HG 参考表，不含 FDA-extension 行）。"
             "解析规则与 step-03 manifest 一致：SMILES 精确匹配，否则作者 cmap_name 名称匹配。"
             "排除对照（vehicle/CTRL）与 l1000 预训练语料。</p>" % N_L1000_DRUGS)

    h.append("<h2>各数据集 in-HG 药物数量</h2>")
    h.append("<table><tr><th>dataset</th><th class='num'>treated 样本数</th>"
             "<th class='num'>in-HG 药物数</th></tr>")
    for _, r in summary_df.iterrows():
        h.append(f"<tr><td>{r['dataset']}</td><td class='num'>"
                 f"{r['n_treated_samples']}</td><td class='num'>"
                 f"{r['n_in_hg_drugs']}</td></tr>")
    h.append(f"<tr><td><b>union（5 数据集并集）</b></td>"
             f"<td class='num' colspan='2'><b>{total_in_hg}</b></td></tr>")
    h.append("</table>")

    h.append("<div class='venn-row'>")
    h.append(f"<figure><img src='data:image/png;base64,{venn_main}'></figure>")
    h.append(f"<figure><img src='data:image/png;base64,{venn_secondary}'>"
             "</figure></div>")

    h.append("<h2>两两交集（in-HG 药物）</h2>")
    h.append(f"<img src='data:image/png;base64,{venn_heat}'>")

    h.append("<h2>按药物的数据集归属</h2>")
    h.append("<table><tr><th class='num'>drug_idx</th><th>drug_name</th>"
             "<th>moa</th>")
    for d in DATASETS:
        h.append(f"<th>{d}</th>")
    h.append("</tr>")
    for r in members.itertuples(index=False):
        h.append("<tr><td class='num'>%d</td><td>%s</td><td>%s</td>" % (
            r.drug_idx, r.drug_name, r.moa))
        for d in DATASETS:
            h.append(f"<td>{getattr(r, d)}</td>")
        h.append("</tr>")
    h.append("</table>")
    h.append("</body></html>")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text("\n".join(h), encoding="utf-8")
    print(f"Drug-overlap report -> {OUT_PATH} "
          f"({len(union)} in-HG drugs, {len(members)} membership rows)")


if __name__ == "__main__":
    main()
