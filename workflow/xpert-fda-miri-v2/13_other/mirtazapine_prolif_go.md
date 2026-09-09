# MIRTAZAPINE — v2 工作流预测心肌细胞转录组对细胞增殖 GO 通路的效应分析
# MIRTAZAPINE — effect on cell-proliferation GO pathways in the v2-predicted cardiomyocyte transcriptome

- 药物 / drug: **MIRTAZAPINE** (v2 drug_idx `2969`, SMILES `CN1CCN2c3ncccc3Cc3ccccc3C2C1`)
- 输入 DEG: `07_xpert_predict/output/xpert_fda-cm_deg-zs.csv` — v2 union-CM (primary+iPSC) finetune 预测、单基因跨 1,004 个药物 z-score（dose_idx 0 / time_idx 2，24 h 暴露）
- Input DEG: v2 step-07 predicted DEG for the union-CM finetune, z-scored per gene across the 1,004 screened drugs (24 h exposure)
- GO 来源: QuickGO — BP 名称含 “proliferation” 的全部条目（459 项），基因注释取 GO:0008283 细胞群体增殖子树下人源蛋白注释（4,670 条，含 IEA/IBA/ISS；剔除 NOT 限定）
- GO source: QuickGO — all BP terms with “proliferation” in the name (459); human annotations under the GO:0008283 subtree (4,670 rows, NOT excluded)
- 打分对象: 仅 L1000 978-landmark 可覆盖（重叠 ≥ 6 才计分；MIRTAZAPINE 共计分 13 个通路）
- Scoring: only GO terms with ≥ 6 L1000-landmark overlap scored (13 terms for MIRTAZAPINE)
- 每通路得分 = 该通路 landmark 基因在药物预测 DEG 上的平均 z；并报告该均值在 1,004 个药物中的百分位。ES 为 GSEA 式 running-sum（与 `_lib/scoring.py` 同统计量）
- pathway score = mean z of the pathway's landmark genes; percentile = rank of that mean inside the 1,004-drug distribution; ES = pipeline GSEA statistic

## 1. 关键参考通路得分 / canonical reference terms

| GO term | 通路 | family | n_lm | mean_z | pctile | 解读 |
|---|---|---|---|---|---|---|
| GO:0008283 | cell population proliferation | 增殖过程 (proliferation process) | 9 | -0.007 | 51% | ≈中性 |
| GO:0008284 | positive regulation of cell population proliferation | 正调控 (positive regulation) | 58 | -0.007 | 53% | ≈中性 |
| GO:0008285 | negative regulation of cell population proliferation | 负调控 (negative regulation) | 38 | -0.049 | 43% | 负值=抑制增强, 此处−0.05≈中性 |
| GO:0042127 | regulation of cell population proliferation | 调控 (regulation, generic) | 19 | +0.138 | 70% | 偏上 / pro |

**要点 / key**: 在可直接计分的通用增殖 term 上，MIRTAZAPINE 预测转录组的正调控池均值≈0（53 分位）、负调控（制动）池均值 −0.05（43 分位，即制动基因整体略下调，方向为 pro）、总调控池 +0.14（70 分位）——即**通用增殖调控轴基本中性、整体方向（若取信）偏向轻微 pro**。

## 2. 心肌 / 平滑肌 / 细胞特异增殖通路（v2 预测转录组背景）

| GO term | 通路 | n_lm | mean_z | pctile | frac_up | ES | 成员基因 (drug_z) |
|---|---|---|---|---|---|---|---|
| GO:0048662 | negative regulation of smooth muscle cell proliferation | 7 | +0.485 | 92% | 0.71 | +0.294 | APOE(+0.47);BMP4(-0.50);HMOX1(+1.55);IGFBP3(-0.12);PPARD(+0.37);PPARG(+0.32);TRIB1(+1.30) |
| GO:0060045 | positive regulation of cardiac muscle cell proliferation | 6 | +0.432 | 94% | 0.83 | -0.384 | CCNB1(+0.80);CDK1(-0.40);FGFR2(+0.35);MEF2C(+1.00);NOTCH1(+0.51);TBX2(+0.34) |
| GO:0048661 | positive regulation of smooth muscle cell proliferation | 7 | +0.400 | 85% | 0.57 | -0.225 | AKT1(-0.67);BMP4(-0.50);EDN1(+1.52);HMOX1(+1.55);ID2(+1.09);MNAT1(-0.32);PIK3CA(+0.14) |
| GO:1904707 | positive regulation of vascular associated smooth muscle cell proliferation | 8 | +0.190 | 64% | 0.75 | +0.332 | DNMT1(-0.67);EDN1(+1.52);FGFR2(+0.35);GNAI2(+0.22);JUN(-1.32);MMP2(+0.35);PAK1(+0.08);TERT(+1.00) |
| GO:0048146 | positive regulation of fibroblast proliferation | 9 | +0.115 | 58% | 0.67 | +0.352 | ABL1(+0.45);CCNA2(+0.96);CDK4(+0.46);CDK6(+0.73);CDKN1A(-1.42);EGFR(-0.32);MIF(+0.39);MYC(-0.67);PDGFA(+0.46) |
| GO:1904706 | negative regulation of vascular associated smooth muscle cell proliferation | 8 | +0.056 | 56% | 0.62 | -0.351 | BMP4(-0.50);CDKN1A(-1.42);CDKN1B(+0.50);GPER1(+0.54);MEF2C(+1.00);PPARG(+0.32);TGFB3(-0.63);TPM1(+0.62) |

其中直接与心肌相关的是 **GO:0060045 positive regulation of cardiac muscle cell proliferation**：mean_z +0.432（94 分位），其 6 个 landmark 里 5 个向上（CCNB1、MEF2C、NOTCH1、TBX2、FGFR2 上调，CDK1 下调）。注意该条目的基因注释均为 IEA/ISS（无实验来源证据），且基因池仅 6 个，解读需谨慎。平滑肌相关条目被 EDN1/HMOX1/ID2 等强上调基因主导（这些基因同样注释到其它增殖与应激过程），故平滑肌正、负调控条目同时偏高，属同源基因池的共线信号而非独立生物学。

## 3. 全部分数（按 mean_z 升序）

| GO term | 通路 | family | n_lm | mean_z | pctile | frac_up |
|---|---|---|---|---|---|---|
| GO:0008285 | negative regulation of cell population proliferation | 负调控 (negative regulation) | 38 | -0.049 | 43% | 0.47 |
| GO:0008284 | positive regulation of cell population proliferation | 正调控 (positive regulation) | 58 | -0.007 | 53% | 0.52 |
| GO:0008283 | cell population proliferation | 增殖过程 (proliferation process) | 9 | -0.007 | 51% | 0.56 |
| GO:1904706 | negative regulation of vascular associated smooth muscle cell proliferation | 负调控 (negative regulation) | 8 | +0.056 | 56% | 0.62 |
| GO:0050679 | positive regulation of epithelial cell proliferation | 正调控 (positive regulation) | 9 | +0.072 | 62% | 0.56 |
| GO:0048146 | positive regulation of fibroblast proliferation | 正调控 (positive regulation) | 9 | +0.115 | 58% | 0.67 |
| GO:0042127 | regulation of cell population proliferation | 调控 (regulation, generic) | 19 | +0.138 | 70% | 0.58 |
| GO:1904707 | positive regulation of vascular associated smooth muscle cell proliferation | 正调控 (positive regulation) | 8 | +0.190 | 64% | 0.75 |
| GO:0048661 | positive regulation of smooth muscle cell proliferation | 正调控 (positive regulation) | 7 | +0.400 | 85% | 0.57 |
| GO:0060045 | positive regulation of cardiac muscle cell proliferation | 正调控 (positive regulation) | 6 | +0.432 | 94% | 0.83 |
| GO:0042100 | B cell proliferation | 增殖过程 (proliferation process) | 6 | +0.484 | 89% | 0.67 |
| GO:0048662 | negative regulation of smooth muscle cell proliferation | 负调控 (negative regulation) | 7 | +0.485 | 92% | 0.71 |
| GO:0030890 | positive regulation of B cell proliferation | 正调控 (positive regulation) | 8 | +0.648 | 85% | 1.00 |

## 4. 家族基因池汇总 / family marker-pool aggregates

同一家族的全部 term 基因并集（基因可在不同细胞类型/正负调控间重复注释，故各池有重叠）。

| family | n_terms | n landmark | mean_z | pctile |
|---|---|---|---|---|
| 负调控 (negative regulation) | 83 | 56 | +0.096 | 70% |
| 正调控 (positive regulation) | 86 | 97 | +0.100 | 69% |
| 增殖过程 (proliferation process) | 160 | 42 | +0.168 | 77% |
| 调控 (regulation, generic) | 130 | 29 | +0.048 | 50% |

四个基因池均值全部为正（mean_z 0.05–0.17），分位 50–77——在 1,004 个筛选药物里，MIRTAZAPINE 的预测转录组让“增殖标志基因池”整体处于**偏上**位置，但幅度很小，且正/负调控池同向，提示这是一般性（弱）激活而非有方向性的增殖失控。

## 5. 药物自身极端基因背景 / the drug's own extreme landmarks

基因水平 z 为跨药物相对值（>0 = 相对全部药物上调）。下列是 MIRTAZAPINE 预测变化最极端的 landmark：

| 上调 top-12 | z | | 下调 top-12 | z |
|---|---|---|---|---|
| ARHGEF12 | +3.06 | | RELB | -3.68 |
| HSPA1A | +2.94 | | DYRK3 | -2.87 |
| TNFRSF21 | +2.42 | | NIT1 | -2.65 |
| SOCS2 | +2.20 | | CD44 | -2.52 |
| INSIG1 | +2.09 | | PMM2 | -2.29 |
| TCTA | +2.08 | | FOS | -2.20 |
| NR2F6 | +1.90 | | TES | -2.07 |
| PRKX | +1.90 | | UBE2L6 | -2.05 |
| TARBP1 | +1.84 | | KLHDC2 | -2.05 |
| KIF20A | +1.68 | | ZNF395 | -1.95 |
| SLC27A3 | +1.67 | | APP | -1.94 |
| ARNT2 | +1.61 | | PTK2 | -1.83 |

## 6. 方法学对照 / method sanity anchor

用同一套通路/打分对已知抗增殖、强心脏毒性药物 DOXORUBICIN 复算：其正调控心脏肌肉细胞增殖 (GO:0060045) mean_z −0.97（0.5 分位）、纤维母细胞增殖正调控 −2.00（1 分位）、家族过程池 −0.64——与已知的强抗增殖一致。MIRTAZAPINE 的相关分数全部为正、幅度小一个量级，说明该打分刻度可区分抗增殖药与弱信号药，MIRTAZAPINE 属于后者（弱 pro 信号）。

## 7. 结论 / conclusions

- v2 预测的 MIRTAZAPINE 心肌细胞 24 h 转录组对**通用细胞增殖 GO 通路净效应≈中性**（正/负调控池均值 −0.05 ~ +0.17 z；50–77 分位）。
- 例外：心肌/平滑肌细胞特异增殖调控条目呈**弱上调**（mean_z +0.4~+0.6，84–94 分位），但驱动基因主要是 EDN1/HMOX1/ID2/MEF2C/CCNB1 等被多个增殖与应激过程共享的强上调基因，且 (i) L1000 覆盖仅 6–8 基因、(ii) 心肌条目证据为 IEA/ISS、(iii) 平滑肌正负条目共线，故该“促增殖”信号**不构成强证据**。
- 文献背景（非心肌细胞）：MIRTAZAPINE 可经 α2-肾上腺素受体拮抗解除去甲肾上腺素对胰岛 β 细胞复制的抑制、并促进 CHO 成纤维细胞增殖（LPA1/ERK/IGF-IR 通路）；在 T 淋巴细胞上呈剂量 依赖性抑制增殖、在肝肿瘤模型中被视为化学预防药（抑制增殖）。无直接心肌细胞增殖证据。
- Literature (non-cardiomyocyte): mirtazapine promotes β-cell replication (α2-adrenergic antagonism, Mol Endocrinol) and CHO-fibroblast proliferation (LPA1/ERK/IGF-IR, Biochem Pharmacol), but is cytostatic in lymphocytes (Drug Chem Toxicol) and anti-proliferative/chemopreventive in a liver-tumour model (Chem Biol Interact); no direct cardiomyocyte report was found.

## 8. 局限 / limitations

- 打分对象仅限 978 个 L1000 landmark；459 个增殖 term 中仅 13 个达到 ≥6 landmark 覆盖，覆盖度不足导致其余通路无法评估（非“无效应”）。
- 基因集合并采用符号匹配，GO 注释中 IEA/IBA/ISS 占相当比例（心肌条目 0 条实验证据）。
- 该 DEG 是模型在 24 h 单时间点的**预测**值（非实测），XPert 本身存在 ~97% 基因预测为下调、药物间方差偏低的已知局限；z 化处理已尽量剥离该基线，但幅度解读仍需谨慎。
- 细胞编码为单一心肌细胞模型（primary+iPSC 并集 finetune，ADR 0015/0016），预测不代表体内心肌细胞实际表型。

---
生成 / generated: `13_other/mirtazapine_prolif_score.py` + `13_other/prolif_go_report.py`; 数据 / data: `13_other/output/mirtazapine_prolif_pathway_scores.csv`, `..._family.csv`, `..._detail.csv`. GO 缓存: `13_other/output/cache/` (QuickGO, 2026-09-08).
