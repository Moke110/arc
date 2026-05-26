# DToxS Level3 本地数据：孔级样本字段覆盖小结

本文档结合 `DToxS_article.md`（*Nature Communications* 2020；人原代 PromoCell 心肌样细胞，3′-DGE RNA-seq，48 h 处理，对照与多种扰动）与当前目录下 **Level3** 配对文件（各子目录中的 `PromoCell-DGE-Experiment-Design*.tsv` + `PromoCell-DGE-Read-Counts*.tsv`）及 `Metadata/` 元数据，对「**单个培养孔**」一条实验设计行（`Species`–`Cell`–`State`–`Plate`–`Well` 等）进行统计。

**样本量**：在 `Data/Level3` 全部 8 份实验设计表中共有 **956** 条孔级记录（与包内 `ReadMe.txt` 所述 914 个 FASTQ 样本数不必一一相等，可能反映不同交付切片或多轮子实验汇总方式）。

**细胞系编码**：部分设计表使用 `PC-A`/`PC-B` 等形式，与 `Primary_Cell_Metadata.txt` 中的 `pc_center_canonical_id`（`A`/`B`/`D`/`E`）对齐时按 **`PC-` 前缀剥离** 后关联。

**对照判定**：`State` 为 `CTRL`、`CTRL-EtOH`、`CTRL-DMSO` 视为**对照/溶媒类**（共 **206** 孔）；其余 **750** 孔视为**给药或其它扰动**。

**扰动元数据匹配**：给药孔的 `State` 三字母（或表内编码）优先在 `Small_Molecule_Metadata.txt` 的 `sm_center_batch_id` 中查找，否则在 `Antibody_Metadata.txt`（`ar_center_batch_id`，如 `TRS`）与 `Protein_Metadata.txt`（`pr_center_canonical_id` / `pr_center_batch_id`，如 `IGF`、`TNF`）中查找。

---

## 1. 转录组（Level3 基因表达矩阵列）

论文方法学：TRIzol 抽提 RNA，3′ 数字基因表达（DGE），比对 hg38，经质控与归一化后得到 fold-change / 读数矩阵；Level3 在本仓库中为 **基因 × 样本孔** 的 `PromoCell-DGE-Read-Counts*.tsv`，列名为孔位 `Well`。

| 项目 | 含 | 不含 |
|------|----|------|
| 转录组（配对读数表中存在该孔对应列） | **956** | **0** |

---

## 2. 培养板号 / 批次号

| 项目 | 含 | 不含 | 说明 |
|------|----|------|------|
| **培养板号**（实验设计 `Plate` 非空） | **956** | **0** | 与 `Well` 共同定位孔位。 |
| **细胞系供应商 Lot / 批次**（`Primary_Cell_Metadata.txt` 中 `pc_provider_batch_id`，与论文 Methods 中细胞 Lot# 一致） | **956** | **0** | 由 `Cell`（规范化后 A/B/D/E）映射。 |
| **测序/处理子实验批次**（数据路径如 `Sequencing Experiment Number N` 及子文件夹名） | **956** | **0** | 每条记录均可归属到具体 Level3 子目录。 |
| **给药孔的试剂供应商批号**（小分子 `sm_provider_batch_id` / 抗体 `ar_provider_batch_id` / 蛋白 `pr_provider_batch_id` 任一类非空） | **721** | **29** | 仅统计 **750** 个给药孔；**29** 孔为 `END`、`SUN+`、`DIC` 等在本包元数据表中**无精确编码匹配**（见下节备注）。 |

若将「培养板号/批次号」理解为 **同时具备板号 +（细胞 Lot 或 测序子实验归属）**，则 **956 / 0**。

---

## 3. 处理条件：对照 / 给药；药物（扰动）信息与 canonical SMILES

### 3.1 对照 vs 给药（全 956 孔）

| 项目 | 孔数 |
|------|------|
| 对照 / 溶媒类（`CTRL`、`CTRL-EtOH`、`CTRL-DMSO`） | **206** |
| 给药或其它扰动（其余 `State`） | **750** |

### 3.2 给药孔（750）：是否能在元数据中定位扰动物

| 项目 | 含 | 不含 |
|------|----|------|
| 扰动信息（小分子 / 抗体 / 蛋白元数据任一命中） | **721** | **29** |

**不含的 29 孔** 的 `State` 均为：`END`（14）、`SUN+`（12）、`DIC`（3）。其中 `END` 与内皮素条目在 `Small_Molecule_Metadata.txt` 中中心批号为 `EDN-1`（canonical `EDN`），**表内编码不一致**；`SUN+` 与 `SUN`（Sunitinib）同理未建别名；`DIC` 在本包小分子表中无同码条目——故严格按字符串匹配计为「不含」。

### 3.3 给药孔（750）：canonical SMILES（LINCS 小分子字段）

本仓库中 **canonical / 批次 SMILES** 仅出现在 `Small_Molecule_Metadata.txt` 的 `sm_smiles_parent` 或 `sm_smiles_batch`（非空且非 `null`）。**抗体/蛋白记录不提供 SMILES**（与论文 Methods 中「小分子指纹用 PubChem SMILES」的适用范围一致）。

| 项目 | 含 | 不含 |
|------|----|------|
| **给药孔**中具有 **非空 SMILES** 的小分子扰动（即命中 `sm_*` 且 SMILES 字段可用） | **672** | **78** |

**不含 78** 的构成：**29** 孔无扰动元数据匹配；**49** 孔为抗体 **TRS**（19）或重组蛋白 **IGF**（16）、**TNF**（14）——有蛋白/抗体元数据与序列类信息，但**无 SMILES 字段**。

---

## 4. 细胞种类

| 项目 | 含 | 不含 |
|------|----|------|
| 细胞种类（`Primary_Cell_Metadata.txt` 中 `pc_cell_type` 非空；论文与元数据均为 **Cardiomyocyte** / 人心肌相关描述） | **956** | **0** |

---

## 5. 与论文表述的对照（简要）

- **转录组**：论文核心为孔级 DGE 谱；本 Level3 读数表与实验设计一一对应列，**全部孔**均有矩阵列，与论文数据类型一致。
- **板号与批次**：论文与 Methods 给出细胞 **Lot#**；本包以 `Primary_Cell_Metadata` 的 `pc_provider_batch_id` 及实验设计 `Plate`/`Well`、目录结构承载技术批次。
- **药物与 SMILES**：论文主文侧重 **激酶抑制剂** 与 FAERS 风险；本 LINCS/DToxS 包另含多种非 KI 小分子、抗体（如 Trastuzumab）及蛋白因子（IGF-1、TNF-α）等，**SMILES 仅覆盖小分子元数据**。

---

*统计脚本逻辑：遍历 `Data/Level3/**/PromoCell-DGE-Experiment-Design*.tsv`，与同目录同名后缀的 `PromoCell-DGE-Read-Counts*.tsv` 首行列名比对 `Well`；元数据表路径为 `Metadata/Small_Molecule_Metadata.txt`、`Antibody_Metadata.txt`、`Protein_Metadata.txt`、`Primary_Cell_Metadata.txt`。*
