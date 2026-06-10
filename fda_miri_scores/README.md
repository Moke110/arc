# MIRI 药物疗效评分：方法框架与实施

> **最后更新**: 2026-06-08 &nbsp;|&nbsp; **状态**: 3 条评分管线已部署并执行
> **依赖**: [XPert](../xpert/) 预测的 FDA 药物心肌细胞转录组 DEG &nbsp;|&nbsp; GEO MIRI 转录组实验数据

---

## 目录

1. [核心原理：转录组连接性](#1-核心原理转录组连接性)
2. [三条评分管线](#2-三条评分管线)
3. [使用方法](#3-使用方法)
4. [评分结果](#4-评分结果)
5. [参考文献](#5-参考文献)

---

## 1. 核心原理：转录组连接性

### 1.1 基本假设

药物疗效通过**转录组连接性**（Transcriptional Connectivity）判断：

- MIRI 导致某些基因**上调**（炎症、凋亡、纤维化相关基因）→ 有效药物应**下调**这些基因
- MIRI 导致某些基因**下调**（线粒体功能、能量代谢、离子通道基因）→ 有效药物应**上调**这些基因
- **药物转录扰动方向与疾病转录失调方向相反 → 药物具有潜在疗效**

### 1.2 输入与输出

```
输入：
  ├── 药物转录组效应 P_drug ∈ R^978    ← XPert 预测的 978 个 landmark 基因 logFC
  └── MIRI 疾病特征 S_disease          ← GEO 实验数据 DEG → human L1000 映射
       ├── UP 基因集 (在 MIRI 中上调)
       ├── DOWN 基因集 (在 MIRI 中下调)
       └── 通路基因集 + 方向 (16 条 MIRI 相关 KEGG 通路)

输出：
  └── 药效评分 s ∈ R    ← 4 个标量评分 + 集成排名，量化药物逆转疾病转录状态的潜力
```

### 1.3 参照数据类型

本项目使用三种类型的疾病参照数据：

| 类型 | 描述 | 格式 | 适用范围 |
|------|------|------|----------|
| **Type A** | 疾病上调/下调基因集 | UP = {gene_1, ...}, DOWN = {gene_1, ...} | WCS, RGES |
| **Type B** | 疾病全基因排序列表 | 按 |disease logFC| 排序 + 带符号秩 | Zhang Score |
| **Type C** | 通路基因集 + 方向 | Pathway_i + direction ∈ {UP, DOWN} | Pathway RGES |

**本项目的参照数据**：

| 来源 | 基因集 | 基因数 | 备注 |
|------|--------|--------|------|
| GSE160516 (mouse) | UP + DOWN | 75 + 32 L1000 | 主要参照 |
| GSE108940 (rat) | UP + DOWN | 4 + 6 L1000 | 阵列注释限制映射 |
| GSE137293 (pig) | UP + DOWN | 2 + 1 L1000 | 无 sham 对照 |
| MIRI_up-down.json | Union | 78 UP + 36 DOWN | 三物种多数投票 |
| MIRI_rank.json | Ranked | 107 genes | Zhang Score Type B 参照 |
| MIRI_pathways.json | Pathways | 16 通路 (14 UP + 2 DOWN) | KEGG + 文献 curated |

---

## 2. 三条评分管线

所有管线均已通过 `scripts/` 中的独立脚本部署，并通过 `scripts/run_all.py` 编排。

### 管线 1: WCS + RGES（基因集法）

- **参考文献**: Subramanian et al. 2017 *Cell*; Chen et al. 2017 *Nat. Commun.*
- **参照类型**: Type A (MIRI_up-down.json)
- **核心算法**: 加权双样本 KS 富集统计量
  - **WCS** (Weighted Connectivity Score): (ES_up − ES_down)/2，仅当 ES_up 与 ES_down 异号时非零
  - **RGES** (Reverse Gene Expression Score): |ES_up| − |ES_down|，始终定义，保留更多信息
- **关键阈值**: WCS < −0.2 且 RGES < −0.2 → 强候选
- **脚本**: `scripts/score_wcs_rges.py`

### 管线 2: Zhang Score / ssCMap（秩乘积法）

- **参考文献**: Zhang & Gant 2008 *BMC Bioinformatics*
- **参照类型**: Type B (MIRI_rank.json, 推荐) 或 Type A
- **核心算法**: 带符号秩乘积
  - 药物基因按 |logFC| 降序赋带符号秩 R(g)
  - 疾病基因预计算带符号秩 s(g)（保留真实疾病效应幅度）
  - c = Σ R(g)·s(g) / C_max ∈ [−1, 1]
- **关键阈值**: c < −0.3 → 强候选
- **脚本**: `scripts/score_zhang.py`

### 管线 3: Pathway-Level RGES（通路法）

- **参考文献**: Chen et al. 2017 *Nat. Commun.* (扩展)
- **参照类型**: Type C (MIRI_pathways.json, 16 通路)
- **核心算法**: 每条通路独立计算 RGES
  - 疾病-UP 通路 → 药物应下调 → efficacy = −ES
  - 疾病-DOWN 通路 → 药物应上调 → efficacy = ES
  - 汇总: mean_efficacy = Σ efficacy / 16
- **关键指标**: mean_efficacy 越高越好。n_effective = efficacy > 0.1 的通路数
- **脚本**: `scripts/score_pathway_rges.py`

### 管线 4: 集成排名

- **方法**: 各管线独立排名求平均
  - rank_WCS, rank_RGES, rank_Zhang: 越低越好（负值 = 逆转）
  - rank_pathway: 升序排名（mean_efficacy 越高越好）
  - rank_integrated = mean(rank_WCS, rank_RGES, rank_Zhang, rank_pathway)
- **脚本**: `scripts/run_all.py`

### 方法对比

| 方法 | 区分度 | 可解释性 | 参考DB需求 | 2023-2025 应用广泛度 |
|------|--------|---------|-----------|---------------------|
| WCS + RGES | ★★★☆☆ | ★★★☆☆ | 无 | ★★★★★ (最广) |
| Zhang Score | ★★★★★ | ★★★☆☆ | 无 | ★★★☆☆ |
| Pathway RGES | ★★★☆☆ | ★★★★★ | 无 | ★★★☆☆ |

---

## 3. 使用方法

```bash
# 完整管线（推荐：使用 z-score 归一化）
cd fda_miri_scores/
python scripts/run_all.py --zscore

# 单独评分
python scripts/score_wcs_rges.py --deg-csv xpert_fda-cm_deg.csv \
    --gene-set reference-genesets/MIRI_up-down.json --zscore
python scripts/score_zhang.py --deg-csv xpert_fda-cm_deg.csv \
    --rank-ref reference-genesets/MIRI_rank.json --zscore
python scripts/score_pathway_rges.py --deg-csv xpert_fda-cm_deg.csv \
    --pathways reference-genesets/MIRI_pathways.json --zscore

# 交互式分析
jupyter notebook scoring.ipynb

# 重新生成参照基因集
python scripts/generate_reference_genesets.py
python scripts/build_union_and_rank.py
```

**评分 Notebook** (`scoring.ipynb`) 包含：
- 完整管线执行（3 种方法）
- 9 张可视化图表（分布、排名、一致性、通路热图等）
- 结果汇总导出

---

## 4. 评分结果

### 4.1 评分区间（z-score 归一化后）

| 方法 | 最小值 | 最大值 | 均值 | 区分度 |
|------|--------|--------|------|--------|
| WCS | −0.408 | +0.396 | −0.027 | 315/3127 非零 |
| RGES | −0.292 | +0.218 | −0.040 | 全区分 |
| Zhang c | −0.604 | +0.443 | −0.076 | 696 个 c < −0.3 |
| Pathway Efficacy | −0.127 | +0.513 | +0.077 | 64 个 efficacy > 0.3 |

### 4.2 跨方法一致性（Spearman ρ）

方法间一致性总体较低（ρ: 0.01–0.27），各方法捕获药物反应的不同方面：
- WCS–Zhang: ρ = 0.046
- WCS–RGES: ρ = 0.230

### 4.3 Top 候选药物

| 排名 | 药物 | 类型 | 通路翻转 (n/16) |
|------|------|------|----------------|
| 1 | Cabozantinib | 多靶点 TKI | 16/16 |
| 2 | Ponatinib | 多靶点 TKI | 14/16 |
| 3 | Ripretinib | 多靶点 TKI | 15/16 |
| 4 | Fostamatinib | SYK 抑制剂 | 15/16 |
| 5 | Cobimetinib | MEK 抑制剂 | 14/16 |

多靶点 TKI 类药物的富集与 MIRI 涉及多通路失调的病理机制一致。

### 4.4 已知局限

- XPert 预测 97.2% 的 DEG 为下调方向，药物间方差极低（per-gene std ≈ 0.01）
- z-score 归一化**必须**使用以求药物区分
- 3 种方法无药物进入共识 top-50（方法互补性高）
- 上游改进方向：提高 XPert 模型的药物间 DEG 方差

---

## 5. 参考文献

1. **Subramanian A**, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*, 2017. DOI: [10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)
2. **Chen B**, et al. Reversal of cancer gene expression correlates with drug efficacy and reveals therapeutic targets. *Nat. Commun.*, 2017. DOI: [10.1038/ncomms16022](https://doi.org/10.1038/ncomms16022)
3. **Zhang SD**, Gant TW. A simple and robust method for connecting small-molecule drugs using gene-expression signatures. *BMC Bioinformatics*, 2008. DOI: [10.1186/1471-2105-9-258](https://doi.org/10.1186/1471-2105-9-258)
4. **Lamb J**, et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science*, 2006. DOI: [10.1126/science.1132939](https://doi.org/10.1126/science.1132939)
5. **Subramanian A**, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 2005. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
6. **Tham E**, Langley SR. Evaluating the robustness of connectivity methods to noise for in silico drug repurposing studies. *Front. Syst. Biol.*, 2022. DOI: [10.3389/fsysb.2022.1050730](https://doi.org/10.3389/fsysb.2022.1050730)
7. **Samart K**, et al. Reconciling multiple connectivity scores for drug repurposing. *Brief. Bioinform.*, 2021. DOI: [10.1093/bib/bbab161](https://doi.org/10.1093/bib/bbab161)

---

> **文件说明**: 本文件为人类可读中文文档。LLM 上下文索引见 [CLAUDE.md](CLAUDE.md)。

