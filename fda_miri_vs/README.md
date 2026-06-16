# MIRI 药物疗效评分：方法框架与实施

> **最后更新**: 2026-06-12 &nbsp;|&nbsp; **状态**: 药物筛选管线已部署（评分→过滤→排名）
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
       └── 通路基因集 + 方向 (60 条 MIRI KEGG/GO 通路，30 上调 + 30 下调)

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
| MIRI_pathways.json | Pathways | 40 通路 (20 UP + 20 DOWN) | KEGG/GO + GEO GSEA + 文献，程序化构建 |

### 1.4 筛选管线

**评分 Notebook** (`fda_screen.ipynb`) 包含三阶段筛选管线：

1. **Scoring** — 3 种方法 × 4 个 Z-score 分量 + 6 条关键通路 efficacy Z-score
2. **Filtering** — WCS < 0 且≥4/6 关键通路 efficacy > 0 → `scores/fda_filtered_scores.csv`
3. **Ranking** — integrated_score = −0.3×RGES_zs − 0.1×Zhang_zs + 0.3×pathway_reversal_zs + 0.05×Σcritical_efficacy_zs → `fda_screened.csv`

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
- **参照类型**: Type C (MIRI_pathways.json, 40 通路: 20 UP + 20 DOWN)
- **核心算法**: 每条通路独立计算 RGES
  - 疾病-UP 通路 → 药物应下调 → efficacy = −ES
  - 疾病-DOWN 通路 → 药物应上调 → efficacy = ES
  - 汇总: mean_efficacy = Σ efficacy / 40, pathway_reversal = n_effective / 40
- **关键指标**: mean_efficacy 越高越好。pathway_reversal = efficacy > 0.1 的通路占比
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

# 交互式分析 / 药物筛选
jupyter notebook fda_screen.ipynb

# 重新生成参照基因集
python scripts/generate_reference_genesets.py
python scripts/build_union_and_rank.py
```

**评分 Notebook** (`fda_screen.ipynb`) 包含：
- 完整管线执行（3 种方法）
- 4 个 Z-score 分量可视化（RainCloud, Z-score distributions）
- 关键通路逆转分布（6 条 critical 通路）
- 过滤：WCS < 0 且≥4/6 关键通路 efficacy > 0
- 排名：integrated_score = −0.3×RGES_zs − 0.1×Zhang_zs + 0.3×pathway_reversal_zs + 0.05×Σcritical_efficacy_zs
- 结果汇总导出（fda_filtered_scores.csv, fda_screened.csv）

---

## 4. 评分结果

### 4.1 评分区间（z-score 归一化后）

| 方法 | 最小值 | 最大值 | 均值 | 区分度 |
|------|--------|--------|------|--------|
| WCS | −0.372 | +0.337 | −0.027 | 315/3127 非零 |
| RGES | −0.393 | +0.258 | −0.040 | 全区分 |
| Zhang c | −0.462 | +0.506 | −0.076 | 696 个 c < −0.3 |
| pathway_reversal | 0.175 | 0.825 | — | 40 通路上计算 |

### 4.2 跨方法一致性（Spearman ρ）

方法间一致性总体较低（ρ: 0.01–0.27），各方法捕获药物反应的不同方面：
- WCS–Zhang: ρ = 0.046
- WCS–RGES: ρ = 0.230

### 4.3 筛选管线

**过滤条件**: WCS < 0 AND ≥4/6 关键通路 efficacy > 0

**通过药物**: 199/3127

6 条关键通路（MIRI_pathways.json, `"critical": true`）:
- KEGG-Apoptosis (UP, 36 genes)
- GO-Cellular Response To Oxidative Stress (UP, 23 genes)
- KEGG-Necroptosis (UP, 24 genes)
- GO-Regulation Of Cell Death (UP, 6 genes)
- KEGG-Calcium signaling pathway (DOWN, 25 genes)
- GO-Mitochondrial Gene Expression (DOWN, 8 genes)

**排名公式**: integrated_score = −0.3×RGES_zs − 0.1×Zhang_zs + 0.3×pathway_reversal_zs + 0.05×Σcritical_efficacy_zs

符号约定: RGES_zs, Zhang_zs 越负越好（逆转疾病）→ 取反（−0.3, −0.1）。pathway_reversal_zs 越正越好（更多通路逆转）→ +0.3。Σcritical_efficacy_zs 越正越好 → +0.05 bonus。

**Top 4 候选药物**:
| 排名 | 药物 | integrated_score | −RGES_zs | −Zhang_zs | PR_zs |
|------|------|-----------------|----------|-----------|-------|
| 1 | Terazosin | +1.429 | +1.507 | +0.973 | +2.603 |
| 2 | Zolpidem Tartrate | +1.316 | +2.350 | +1.376 | +1.162 |
| 3 | Cefiderocol Sulfate Tosylate | +1.312 | +2.085 | +0.672 | +1.402 |
| 4 | Loracarbef | +1.303 | +2.154 | +0.196 | +1.402 |
| 5 | Erythromycin Glucepate | +1.294 | +3.611 | +0.691 | −0.039 |

**输出文件**:
- `scores/fda_filtered_scores.csv` — 通过过滤的药物 + Z-score 各分量
- `fda_screened.csv` — 最终排名候选列表

### 4.4 已知局限

- XPert 预测 97.2% 的 DEG 为下调方向，药物间方差极低（per-gene std ≈ 0.01）
- z-score 归一化**必须**使用以求药物区分
- 3 种方法无药物进入共识 top-50（方法互补性高）
- 上游改进方向：提高 XPert 模型的药物间 DEG 方差

---

## 5. MIRI 通路重建（2026-06-11）

### 5.1 背景与目的

原 MIRI_pathways.json 含 16 条手动 curated 通路（14 上调 + 2 下调）。为提供更系统、可复现的通路参照，基于 GEO 数据与文献调研程序化重建通路基因集。

### 5.2 数据来源

| 来源 | 描述 | 通路数 |
|------|------|--------|
| **KEGG_2021_Human** (Enrichr) | 人类 KEGG 通路基因集（HUGO 符号） | 320 |
| **GO_Biological_Process_2023** (Enrichr) | 人类 GO BP 术语基因集（HUGO 符号） | 5,407 |
| **GSE160516** (mouse MIRI) | 主要 GEO 数据集，20,710 DEGs，mygene+homologene → human | 权重 1.0 |
| **GSE108940** (rat MIRI) | 辅助验证，28,036 DEGs | 权重 0.5 |
| **GSE137293** (pig scar vs border) | 辅助验证，13,117 DEGs | 权重 0.3 |
| **gpt-researcher** | 文献调研，深度搜索 MIRI 通路文献 | 补充验证 |

### 5.3 构建流程

```
1. 从 Enrichr 获取 KEGG/GO 基因集 → 缓存至 .claude/
2. 加载 GEO DEG → mygene+homologene 映射至 human → 按 logFC 排序
3. Preranked GSEA (enrichment_score) → 每条通路计算 ES + 预测方向
4. 与 L1000 978 landmark 基因取交集
5. 双重筛选：
   a. n_genes > 5（L1000 交集后）
   b. n_l1000/n_total ≥ 978/20000 = 0.0489（覆盖率不低于全基因组均值）
6. 合并基因集完全相同的通路（KEGG+GO → KEGG&GO-）
7. 综合评分排序（GEO ES × 数据集权重 + 文献支持 + 基因数奖励）
8. 两阶段选取前 ~30 个上调 + ~30 个下调通路（文献支持通路优先，GEO 填充补齐）
9. 输出 MIRI_pathways.json
```

### 5.4 筛选参数

| 参数 | 值 | 说明 |
|------|-----|------|
| total_human_genome_genes | 20,000 | 人类蛋白编码基因估计数 |
| n_landmark | 978 | L1000 landmark 基因数 |
| coverage_factor | 1.0 | 不放宽覆盖率要求 |
| coverage_threshold | 0.0489 | 978/20000 |
| min_l1000_genes | 6 (>5) | 与 utils.py 阈值一致 |

### 5.5 重建结果

- **候选通路总数**: 8,619（来自 GEO KEGG + GO GSEA + 文献）
- **通过筛选**: 1,660 条
- **合并**: 53 条通路合并至已有通路
- **最终选取**: 40 条（20 UP + 20 DOWN），其中 20 条有文献支持（13 UP + 7 DOWN）
- **选取策略**: 两阶段优先制 — 文献支持通路优先入选，GEO 评分通路填充剩余名额
- **上调通路主题**: 炎症/免疫信号（TNF、NF-κB、MAPK、趋化因子）、凋亡/坏死性凋亡、氧化应激、伤口愈合
- **下调通路主题**: 心肌收缩、钙信号、肾上腺素信号、粘着斑/ECM、肌动蛋白骨架、扩张型心肌病

#### 5.5.1 上调通路（20 条，MIRI 中激活，药物应下调）

★ = 有文献支持

| ★ | 通路名称 | 基因数 | 简介 |
|----|----------|--------|------|
| ★ | KEGG-Chemokine signaling pathway | 35 | 趋化因子信号通路，介导中性粒细胞、单核细胞向缺血心肌趋化迁移，MIRI 无菌性炎症的核心驱动通路 |
| ★ | GO-Wound Healing | 15 | 伤口愈合过程，包含凝血、纤维蛋白溶解、ECM 重塑等 MIRI 后期组织修复机制 |
| ★ | GO-Inflammatory Response | 20 | 炎症反应，MIRI 中最显著上调的生物学过程，涵盖 TNF、IL-1、IL-6 等促炎因子的合成与分泌 |
| ★ | GO-Regulation Of Inflammatory Response | 22 | 炎症反应调控，包含 TRAF1、ZNF692 等 MIRI 炎症调节因子，地塞米松的潜在作用靶点 |
| ★ | GO-Cellular Response To Tumor Necrosis Factor | 21 | 细胞对 TNF 的应答，NF-κB 通路核心下游，Ting et al. (2026) 实验验证的 MIRI 生物标志物通路 |
| ★ | GO-Cytokine-Mediated Signaling Pathway | 27 | 细胞因子介导的信号通路，涵盖 IL-6、IL-8、IL-1β 等 MIRI 关键促炎因子的信号传导 |
| ★ | KEGG-TNF signaling pathway | 23 | TNF 信号通路，MIRI 炎症核心轴，Nfkbia 和 Icam1 上调经 RT-qPCR 验证（Ting et al., 2026） |
| ★ | KEGG-Apoptosis | 36 | 凋亡通路，MIRI 中促凋亡调节因子和执行 Caspase 上调，导致心肌细胞丢失 |
| ★ | GO-Cellular Response To Oxidative Stress | 23 | 细胞氧化应激应答，反映 MIRI 再灌注时 ROS 爆发引起的抗氧化和应激基因表达变化 |
| ★ | KEGG-Necroptosis | 24 | 坏死性凋亡通路，MIRI 中程序性坏死相关基因富集，较凋亡更具促炎特性 |
| ★ | GO-Positive Regulation Of NF-kappaB Transcription Factor Activity | 17 | NF-κB 转录因子活性正调控，RAGE/NF-κB/NLRP3 炎症轴的关键环节（Liu Y. et al., 2025） |
| ★ | KEGG-MAPK signaling pathway | 54 | MAPK 信号通路，与 TNF/NF-κB 协同激活，S100A8/A9 通过 p38/JNK 促进心肌肥大 |
| ★ | GO-Regulation Of Cell Death | 6 | 细胞死亡调控，涵盖凋亡和坏死性凋亡的交叉调节机制 |
| | GO-Cellular Response To Lectin | 6 | 细胞对凝集素的响应，C 型凝集素受体（CLR）信号介导 MIRI 固有免疫激活 |
| | GO-Fc Receptor Mediated Stimulatory Signaling Pathway | 8 | Fc 受体介导的刺激性信号通路，参与免疫复合物介导的心肌炎症 |
| | GO-Positive Regulation Of Nitric Oxide Metabolic Process | 6 | 一氧化氮代谢正调控，iNOS 催化产生过量 NO 导致 MIRI 氧化/硝化损伤 |
| | GO-Negative Regulation Of Insulin Receptor Signaling Pathway | 8 | 胰岛素受体信号负调控，MIRI 中胰岛素抵抗加重心肌能量代谢障碍 |
| | GO-Positive Regulation Of Interleukin-1 Beta Production | 6 | IL-1β 生成正调控，NLRP3 炎性小体激活后 IL-1β 分泌增多，放大炎症级联 |
| | GO-Positive Regulation Of Tumor Necrosis Factor Production | 8 | TNF 生成正调控，DAMP 分子 S100A8/A9 通过 TLR4/RAGE 刺激 TNF 分泌 |
| | GO-Extracellular Structure Organization | 6 | 细胞外结构组织，MIRI 后 ECM 重塑的早期阶段，涉及胶原和纤维连接蛋白 |

#### 5.5.2 下调通路（20 条，MIRI 中抑制，药物应上调）

| ★ | 通路名称 | 基因数 | 简介 |
|----|----------|--------|------|
| ★ | KEGG-Adrenergic signaling in cardiomyocytes | 16 | 心肌细胞肾上腺素能信号，β-AR 及其下游效应分子表达下降，导致心肌收缩功能障碍 |
| ★ | GO-Regulation Of Protein Kinase Activity | 20 | 蛋白激酶活性调控，多种激酶级联负调控因子下调，影响心肌保护信号传导 |
| ★ | KEGG-Focal adhesion | 37 | 粘着斑，心肌细胞-ECM 机械力传导的关键结构，下调导致心肌结构完整性受损 |
| ★ | KEGG-Regulation of actin cytoskeleton | 30 | 肌动蛋白骨架调控，心肌细胞骨架紊乱的转录组反映，影响收缩单元排列 |
| ★ | KEGG-Calcium signaling pathway | 25 | 钙信号通路，钙离子处理蛋白表达下降，经典"钙超载"范式的转录组对应 |
| ★ | GO-Extracellular matrix organization | 10 | ECM 组织，胶原和纤维连接蛋白基因网络被抑制，基质金属蛋白酶参与心室重构 |
| ★ | KEGG-Dilated cardiomyopathy | 7 | 扩张型心肌病，心肌病相关基因表达降低，反映 MIRI 后收缩功能损伤 |
| | KEGG-Valine, leucine and isoleucine degradation | 9 | 支链氨基酸（BCAA）降解，MIRI 中代谢从脂肪酸氧化转向糖酵解，线粒体 BCAA 分解减弱 |
| | GO-Peptidyl-Threonine Phosphorylation | 9 | 肽基苏氨酸磷酸化，涉及 mTOR/p70S6K 等代谢调节激酶的苏氨酸修饰，MIRI 中信号减弱 |
| | GO-Cellular Response To Retinoic Acid | 6 | 细胞对视黄酸的应答，视黄酸信号参与心肌分化与代谢调节，MIRI 中受抑 |
| | GO-Fatty Acid Oxidation | 7 | 脂肪酸氧化，心肌主要能量来源，MIRI 中脂质氧化酶基因下调，代谢底物利用障碍 |
| | GO-Mitochondrial Gene Expression | 8 | 线粒体基因表达，氧化磷酸化复合体亚基转录下降，线粒体功能障碍的核心表现 |
| | GO-Fatty Acid Catabolic Process | 8 | 脂肪酸分解代谢，涵盖 FAO 和脂解的多个步骤，MIRI 中脂质代谢全面抑制 |
| | GO-Positive Regulation Of Cyclin-Dependent Protein Kinase Activity | 6 | CDK 活性正调控，心肌细胞周期停滞，终末分化心肌细胞丧失增殖能力 |
| | GO-Regulation Of G0 To G1 Transition | 6 | G0/G1 期转换调控，DNA 复制和细胞周期基因下调，MIRI 后心肌再生受限 |
| | KEGG-Pyruvate metabolism | 6 | 丙酮酸代谢，糖酵解产物进入 TCA 循环的节点，MIRI 中丙酮酸脱氢酶活性受抑 |
| | KEGG-Type II diabetes mellitus | 7 | II 型糖尿病通路，胰岛素信号和葡萄糖代谢基因下调，代谢综合征加重 MIRI 易感性 |
| | GO-Fatty Acid Beta-Oxidation | 7 | 脂肪酸 β-氧化，线粒体基质中的脂肪酸降解步骤，MIRI 中最显著下调的代谢通路之一 |
| | GO-Positive Regulation Of Stem Cell Population Maintenance | 6 | 干细胞群体维持正调控，心肌祖细胞维持信号下调，影响内源性心肌修复 |
| | GO-Positive Regulation Of Lipid Biosynthetic Process | 6 | 脂质生物合成正调控，胆固醇和磷脂合成基因下调，膜修复能力减弱 |

### 5.6 使用方法

```bash
# 重建通路（含文献调研）
python scripts/rebuild_pathways.py --target-up 30 --target-down 30

# 仅 GEO 富集（跳过文献调研）
python scripts/rebuild_pathways.py --skip-literature

# 强制刷新 Enrichr 缓存
python scripts/rebuild_pathways.py --refresh-cache

# 自定义人类基因组大小
python scripts/rebuild_pathways.py --total-human-genes 25000
```

中间产物保存至 `fda_miri_scores/.claude/`。

---

## 6. 参考文献

1. **Subramanian A**, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*, 2017. DOI: [10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)
2. **Chen B**, et al. Reversal of cancer gene expression correlates with drug efficacy and reveals therapeutic targets. *Nat. Commun.*, 2017. DOI: [10.1038/ncomms16022](https://doi.org/10.1038/ncomms16022)
3. **Zhang SD**, Gant TW. A simple and robust method for connecting small-molecule drugs using gene-expression signatures. *BMC Bioinformatics*, 2008. DOI: [10.1186/1471-2105-9-258](https://doi.org/10.1186/1471-2105-9-258)
4. **Lamb J**, et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science*, 2006. DOI: [10.1126/science.1132939](https://doi.org/10.1126/science.1132939)
5. **Subramanian A**, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 2005. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
6. **Tham E**, Langley SR. Evaluating the robustness of connectivity methods to noise for in silico drug repurposing studies. *Front. Syst. Biol.*, 2022. DOI: [10.3389/fsysb.2022.1050730](https://doi.org/10.3389/fsysb.2022.1050730)
7. **Samart K**, et al. Reconciling multiple connectivity scores for drug repurposing. *Brief. Bioinform.*, 2021. DOI: [10.1093/bib/bbab161](https://doi.org/10.1093/bib/bbab161)

---

> **文件说明**: 本文件为人类可读中文文档。LLM 上下文索引见 [CLAUDE.md](CLAUDE.md)。

