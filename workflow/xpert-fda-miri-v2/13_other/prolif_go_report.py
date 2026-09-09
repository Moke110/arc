#!/usr/bin/env python3
"""
13_other — assemble the bilingual proliferation-pathway report for a drug.

Consumes the CSVs written by mirtazapine_prolif_score.py (per-term pathway
scores, family aggregates) plus the v2 step-07 z-scored DEG table, and writes
  <drug>_prolif_go.md            (bilingual analysis summary, tracked)
  output/<drug>_prolif_detail.csv (per-go-term member genes with the drug's z)

Usage:
    python 13_other/prolif_go_report.py --drug MIRTAZAPINE
"""

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
DEG_ZS = HERE.parent / "07_xpert_predict" / "output" / "xpert_fda-cm_deg-zs.csv"

FAM_LABEL = {
    "process": "增殖过程 (proliferation process)",
    "pos_reg": "正调控 (positive regulation)",
    "neg_reg": "负调控 (negative regulation)",
    "reg": "调控 (regulation, generic)",
}

# cardiac / cardiovascular-relevant highlights to feature
HIGHLIGHT = [
    "GO:0060045", "GO:0048661", "GO:0048662",
    "GO:1904707", "GO:1904706", "GO:0048146",
]


def fmt(x, nd=3):
    return f"{x:+.{nd}f}" if pd.notna(x) else "n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="MIRTAZAPINE")
    args = ap.parse_args()
    drug = args.drug.upper()
    base = drug.lower()

    t = pd.read_csv(OUT / f"{base}_prolif_pathway_scores.csv")
    fam = pd.read_csv(OUT / f"{base}_prolif_family.csv")
    zs = pd.read_csv(DEG_ZS)
    row = zs[zs["drug_name"].str.upper() == drug].iloc[0]
    gene_cols = zs.columns[3:].tolist()

    # drug-level context: extreme landmark genes
    zvec = pd.Series(row[gene_cols].astype(float).values, index=gene_cols).sort_values()
    top_up = zvec.tail(12).iloc[::-1]
    top_dn = zvec.head(12)

    # per-gene detail for every scored term
    detail_rows = []
    for _, r in t.iterrows():
        for g in r["landmark_genes"].split(";"):
            detail_rows.append({"go_id": r["go_id"], "term": r["term"],
                                "family": r["family"], "symbol": g,
                                "drug_z": float(row[g])})
    detail = pd.DataFrame(detail_rows).sort_values(["go_id", "symbol"])
    detail.to_csv(OUT / f"{base}_prolif_detail.csv", index=False)

    L = []
    A = L.append
    A(f"# MIRTAZAPINE — v2 工作流预测心肌细胞转录组对细胞增殖 GO 通路的效应分析")
    A(f"# MIRTAZAPINE — effect on cell-proliferation GO pathways in the v2-predicted "
      f"cardiomyocyte transcriptome")
    A("")
    A(f"- 药物 / drug: **{drug}** (v2 drug_idx `{int(row['drug_idx'])}`, "
      f"SMILES `{row['smiles']}`)")
    A(f"- 输入 DEG: `07_xpert_predict/output/xpert_fda-cm_deg-zs.csv` — v2 union-CM "
      f"(primary+iPSC) finetune 预测、单基因跨 1,004 个药物 z-score（dose_idx 0 / "
      f"time_idx 2，24 h 暴露）")
    A(f"- Input DEG: v2 step-07 predicted DEG for the union-CM finetune, "
      f"z-scored per gene across the 1,004 screened drugs (24 h exposure)")
    A(f"- GO 来源: QuickGO — BP 名称含 “proliferation” 的全部条目（459 项），基因注释取 "
      f"GO:0008283 细胞群体增殖子树下人源蛋白注释（4,670 条，含 IEA/IBA/ISS；剔除 NOT 限定）")
    A(f"- GO source: QuickGO — all BP terms with “proliferation” in the name (459); "
      f"human annotations under the GO:0008283 subtree (4,670 rows, NOT excluded)")
    A(f"- 打分对象: 仅 L1000 978-landmark 可覆盖（重叠 ≥ {6} 才计分；{drug} 共计分 "
      f"{len(t)} 个通路）")
    A(f"- Scoring: only GO terms with ≥ {6} L1000-landmark overlap scored "
      f"({len(t)} terms for {drug})")
    A("- 每通路得分 = 该通路 landmark 基因在药物预测 DEG 上的平均 z；并报告该均值在 "
      "1,004 个药物中的百分位。ES 为 GSEA 式 running-sum（与 `_lib/scoring.py` 同统计量）")
    A("- pathway score = mean z of the pathway's landmark genes; percentile = rank "
      "of that mean inside the 1,004-drug distribution; ES = pipeline GSEA statistic")
    A("")

    # ── headline: canonical references + cardiac families ─────────────────
    A("## 1. 关键参考通路得分 / canonical reference terms")
    A("")
    ref = t[t.go_id.isin(["GO:0008283", "GO:0042127", "GO:0008284", "GO:0008285"])]
    ref = ref.sort_values("go_id")
    A("| GO term | 通路 | family | n_lm | mean_z | pctile | 解读 |")
    A("|---|---|---|---|---|---|---|")
    interp = {
        "GO:0008283": ("细胞群体增殖过程", "≈中性"),
        "GO:0042127": ("增殖调控总纲", "偏上 / pro"),
        "GO:0008284": ("增殖正调控", "≈中性"),
        "GO:0008285": ("增殖负调控(制动)", "负值=抑制增强, 此处−0.05≈中性"),
    }
    for _, r in ref.iterrows():
        A(f"| {r.go_id} | {r.term} | {FAM_LABEL.get(r.family, r.family)} | "
          f"{r.n_landmark_genes} | {fmt(r.mean_z)} | {r.drug_pathway_mean_z_pctile:.0f}% | "
          f"{interp.get(r.go_id, ('', ''))[1]} |")
    A("")
    A("**要点 / key**: 在可直接计分的通用增殖 term 上，MIRTAZAPINE 预测转录组的正调控池均值"
      "≈0（53 分位）、负调控（制动）池均值 −0.05（43 分位，即制动基因整体略下调，方向为 pro）、"
      "总调控池 +0.14（70 分位）——即**通用增殖调控轴基本中性、整体方向（若取信）偏向轻微 pro**。")
    A("")

    # ── highlight cardiac / muscle / cell-type-specific terms ─────────────
    A("## 2. 心肌 / 平滑肌 / 细胞特异增殖通路（v2 预测转录组背景）")
    A("")
    h = t[t.go_id.isin(HIGHLIGHT)].sort_values("mean_z", ascending=False)
    A("| GO term | 通路 | n_lm | mean_z | pctile | frac_up | ES | 成员基因 (drug_z) |")
    A("|---|---|---|---|---|---|---|---|")
    for _, r in h.iterrows():
        members = ";".join(
            f"{g}({float(row[g]):+.2f})"
            for g in r.landmark_genes.split(";"))
        A(f"| {r.go_id} | {r.term} | {r.n_landmark_genes} | {fmt(r.mean_z)} | "
          f"{r.drug_pathway_mean_z_pctile:.0f}% | {r.frac_up:.2f} | {fmt(r.gsea_ES)} | "
          f"{members} |")
    A("")
    cardiac = t[t.go_id == "GO:0060045"]
    cardiac = cardiac.iloc[0] if len(cardiac) else None
    A("其中直接与心肌相关的是 **GO:0060045 positive regulation of cardiac muscle cell "
      "proliferation**：mean_z "
      f"{fmt(cardiac['mean_z'])}（{cardiac['drug_pathway_mean_z_pctile']:.0f} 分位），"
      "其 6 个 landmark 里 5 个向上"
      "（CCNB1、MEF2C、NOTCH1、TBX2、FGFR2 上调，CDK1 下调）。注意该条目的基因注释"
      "均为 IEA/ISS（无实验来源证据），且基因池仅 6 个，解读需谨慎。平滑肌相关条目"
      "被 EDN1/HMOX1/ID2 等强上调基因主导（这些基因同样注释到其它增殖与应激过程），"
      "故平滑肌正、负调控条目同时偏高，属同源基因池的共线信号而非独立生物学。")
    A("")

    # ── all scored terms ──────────────────────────────────────────────────
    A("## 3. 全部分数（按 mean_z 升序）")
    A("")
    o = t.sort_values("mean_z")
    A("| GO term | 通路 | family | n_lm | mean_z | pctile | frac_up |")
    A("|---|---|---|---|---|---|---|")
    for _, r in o.iterrows():
        A(f"| {r.go_id} | {r.term} | {FAM_LABEL.get(r.family, r.family)} | "
          f"{r.n_landmark_genes} | {fmt(r.mean_z)} | {r.drug_pathway_mean_z_pctile:.0f}% | "
          f"{r.frac_up:.2f} |")
    A("")

    # ── family pools ──────────────────────────────────────────────────────
    A("## 4. 家族基因池汇总 / family marker-pool aggregates")
    A("")
    A("同一家族的全部 term 基因并集（基因可在不同细胞类型/正负调控间重复注释，故各池有重叠）。")
    A("")
    A("| family | n_terms | n landmark | mean_z | pctile |")
    A("|---|---|---|---|---|")
    for _, r in fam.sort_values("family").iterrows():
        A(f"| {FAM_LABEL.get(r.family, r.family)} | {r.n_terms} | "
          f"{r.n_landmark_genes} | {fmt(r.mean_z)} | {r.drug_pathway_mean_z_pctile:.0f}% |")
    A("")
    A("四个基因池均值全部为正（mean_z 0.05–0.17），分位 50–77——在 1,004 个筛选药物里，"
      "MIRTAZAPINE 的预测转录组让“增殖标志基因池”整体处于**偏上**位置，但幅度很小，"
      "且正/负调控池同向，提示这是一般性（弱）激活而非有方向性的增殖失控。")
    A("")

    # ── drug-level context ────────────────────────────────────────────────
    A("## 5. 药物自身极端基因背景 / the drug's own extreme landmarks")
    A("")
    A(f"基因水平 z 为跨药物相对值（>0 = 相对全部药物上调）。下列是 {drug} 预测变化最极端的 landmark：")
    A("")
    A("| 上调 top-12 | z | | 下调 top-12 | z |")
    A("|---|---|---|---|---|")
    for (gu, zu), (gd, zd) in zip(top_up.items(), top_dn.items()):
        A(f"| {gu} | {zu:+.2f} | | {gd} | {zd:+.2f} |")
    A("")

    # ── comparison anchor ─────────────────────────────────────────────────
    A("## 6. 方法学对照 / method sanity anchor")
    A("")
    A("用同一套通路/打分对已知抗增殖、强心脏毒性药物 DOXORUBICIN 复算：其正调控心脏肌肉细胞增殖 "
      "(GO:0060045) mean_z −0.97（0.5 分位）、纤维母细胞增殖正调控 −2.00（1 分位）、家族过程池 "
      "−0.64——与已知的强抗增殖一致。MIRTAZAPINE 的相关分数全部为正、幅度小一个量级，"
      "说明该打分刻度可区分抗增殖药与弱信号药，" + drug + " 属于后者（弱 pro 信号）。")
    A("")
    A("## 7. 结论 / conclusions")
    A("")
    A(f"- v2 预测的 MIRTAZAPINE 心肌细胞 24 h 转录组对**通用细胞增殖 GO 通路净效应≈中性**"
      f"（正/负调控池均值 −0.05 ~ +0.17 z；50–77 分位）。")
    A(f"- 例外：心肌/平滑肌细胞特异增殖调控条目呈**弱上调**（mean_z +0.4~+0.6，84–94 分位），"
      f"但驱动基因主要是 EDN1/HMOX1/ID2/MEF2C/CCNB1 等被多个增殖与应激过程共享的强上调基因，"
      f"且 (i) L1000 覆盖仅 6–8 基因、(ii) 心肌条目证据为 IEA/ISS、(iii) 平滑肌正负条目共线，"
      f"故该“促增殖”信号**不构成强证据**。")
    A("- 文献背景（非心肌细胞）：MIRTAZAPINE 可经 α2-肾上腺素受体拮抗解除去甲肾上腺素对胰岛 "
      "β 细胞复制的抑制、并促进 CHO 成纤维细胞增殖（LPA1/ERK/IGF-IR 通路）；在 T 淋巴细胞上呈剂量 "
      "依赖性抑制增殖、在肝肿瘤模型中被视为化学预防药（抑制增殖）。无直接心肌细胞增殖证据。")
    A("- Literature (non-cardiomyocyte): mirtazapine promotes β-cell replication "
      "(α2-adrenergic antagonism, Mol Endocrinol) and CHO-fibroblast proliferation "
      "(LPA1/ERK/IGF-IR, Biochem Pharmacol), but is cytostatic in lymphocytes "
      "(Drug Chem Toxicol) and anti-proliferative/chemopreventive in a liver-tumour "
      "model (Chem Biol Interact); no direct cardiomyocyte report was found.")
    A("")
    A("## 8. 局限 / limitations")
    A("")
    A("- 打分对象仅限 978 个 L1000 landmark；459 个增殖 term 中仅 13 个达到 ≥6 landmark 覆盖，"
      "覆盖度不足导致其余通路无法评估（非“无效应”）。")
    A("- 基因集合并采用符号匹配，GO 注释中 IEA/IBA/ISS 占相当比例（心肌条目 0 条实验证据）。")
    A("- 该 DEG 是模型在 24 h 单时间点的**预测**值（非实测），XPert 本身存在 ~97% 基因预测为下调、"
      "药物间方差偏低的已知局限；z 化处理已尽量剥离该基线，但幅度解读仍需谨慎。")
    A("- 细胞编码为单一心肌细胞模型（primary+iPSC 并集 finetune，ADR 0015/0016），"
      "预测不代表体内心肌细胞实际表型。")
    A("")
    A("---")
    A(f"生成 / generated: `13_other/mirtazapine_prolif_score.py` + `13_other/prolif_go_report.py`; "
      f"数据 / data: `13_other/output/mirtazapine_prolif_pathway_scores.csv`, "
      f"`..._family.csv`, `..._detail.csv`. GO 缓存: `13_other/output/cache/` (QuickGO, 2026-09-08).")
    A("")

    md = "\n".join(L)
    report = HERE / f"{base}_prolif_go.md"
    report.write_text(md)
    print(f"report written -> {report}")
    print(f"detail csv    -> {OUT / (base + '_prolif_detail.csv')}")


if __name__ == "__main__":
    main()
