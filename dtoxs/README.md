# DToxS 数据集文档

## 数据集概述

DToxS（Drug Toxicity Signature，LINCS 数据集编号 LDS-1593）是一个人源原代心肌细胞药物扰动转录组数据集。该数据集由 Nature Communications 于 2020 年发表，旨在通过高通量转录组学方法评估 FDA 批准药物（尤其是激酶抑制剂）对心肌细胞的潜在毒性。

本数据集在当前项目中的角色：作为 XPert 药物转录组扰动预测模型的**微调目标数据集**，使模型从通用 L1000 转录组空间适配到心肌细胞特异性药物响应预测空间。

## 数据来源与实验设计

### 基本信息

| 项目 | 内容 |
|------|------|
| LINCS 数据集 ID | LDS-1593 |
| GEO 登录号 | GSE146096, GSE146097 |
| 测序平台 | 3'-DGE RNA-seq（Illumina HiSeq 2500） |
| 读长 | 单端 46 bp |
| 文库构建 | TruSeq Stranded mRNA Library Prep Kit |
| 数据获取日期 | 2020-10-17 |
| 参考文献 | Nature Communications (2020)，详见 `DToxS_article.md` |

### 细胞信息

| 项目 | 内容 |
|------|------|
| 细胞类型 | 人源原代心肌细胞（Human Primary Cardiomyocyte） |
| 供应商 | PromoCell（货号 C-12810） |
| 供体系数 | 4 条（PC-A, PC-B, PC-D, PC-E） |
| 供体性别 | 2 男 (A, E)，2 女 (B, D) |
| 供体年龄 | 54–62 岁 |
| 供体种族 | 高加索人（Caucasian） |
| 健康状态 | 健康（healthy），无已知疾病 |

### 实验处理

| 项目 | 内容 |
|------|------|
| 处理方式 | 药物暴露 48 小时 |
| FASTQ 样本数 | 914 个 |
| 孔级记录数 | 956 条（实验设计表中记录） |
| 对照孔 | 206 个（CTRL, CTRL-DMSO, CTRL-EtOH） |
| 给药孔 | 750 个（小分子 + 抗体 + 蛋白） |
| 实验批次 | 5 个测序实验（Exp 1, 2, 4, 5；其中 Exp 5 含 Cor / Spec-Fixed / Spec-Var-1/2/3 子实验） |
| 归一化方法 | RMA 归一化（log2 表达值） |

### 扰动类型分布

| 扰动类型 | 孔数 | 说明 |
|----------|------|------|
| 小分子药物 | 672 | 具有有效 canonical SMILES，36 种小分子 |
| 抗体 | 19 | TRS（Trastuzumab） |
| 重组蛋白 | 16 + 14 | IGF-1（16 孔），TNF-alpha（14 孔） |
| 无元数据匹配 | 29 | END（14 孔，应为 EDN/内皮素）、SUN+（12 孔，Sunitinib 变体代码）、DIC（3 孔，Diclofenac？无匹配） |
| **给药孔合计** | **750** | |

## 目录结构

```
dtoxs/
├── Data/
│   ├── Level3/                                    # 原始 read count 数据 (TSV)
│   │   ├── Sequencing Experiment Number 1/        # 每个子实验含实验设计表 + 读数表
│   │   ├── Sequencing Experiment Number 2/
│   │   ├── Sequencing Experiment Number 4/
│   │   └── Sequencing Experiment Number 5/
│   │       ├── Cor/                               # 核心实验
│   │       ├── Spec-Fixed/
│   │       ├── Spec-Var-1/
│   │       ├── Spec-Var-2/
│   │       └── Spec-Var-3/
│   └── PromoCell-DGE-RNAseq-Experiment-Design.DCIC.xlsx  # 总体实验设计 Excel
├── Metadata/
│   ├── Small_Molecule_Metadata.txt                # 小分子元数据（48 种，含 SMILES/InChI）
│   ├── Antibody_Metadata.txt                      # 抗体元数据（如 Trastuzumab）
│   ├── Protein_Metadata.txt                       # 蛋白元数据（如 IGF-1, TNF-alpha）
│   ├── Primary_Cell_Metadata.txt                  # 原代细胞供体元数据（4 条）
│   ├── ExperimentalMetadata.txt                   # 实验参数（仪器、日期等）
│   ├── DatasetMetadata.txt                        # 数据集总体描述
│   ├── MetadataPropertyMapping.txt                # 元数据字段说明映射
│   ├── UCSC-hg38-Gene-HGNC-Entrez-Map.tsv          # 基因符号/HGNC/Entrez ID 映射表
│   └── small_mol_category.md                      # 36 种小分子按化学结构分 6 类
├── ProcessingPipeline/
│   └── ProcessingPipelineSpecifications.txt       # 原始数据处理管线说明
├── DToxS_article.md                               # 原始论文摘要
├── DToxS_summary.md                               # 孔级元数据覆盖度分析报告
├── ReadMe.txt                                     # GEO 原始 README
├── README.md                                      # 本文件（中文使用文档）
└── CLAUDE.md                                      # LLM 上下文索引
```

## 小分子药物分类

DToxS 中 36 种小分子按化学母核/药理合并为 **6 大类**，用于 XPert 微调的 cold-cate（按类别留出）交叉验证策略。

### Category 1: 喹唑啉类激酶抑制剂 (Quinazoline-based KIs)

共享 4-苯胺基喹唑啉 (4-anilinoquinazoline) 母核，主要靶向 EGFR/HER2。

| 代码 | 名称 | 靶点 | 样本数 |
|------|------|------|--------|
| AFA | Afatinib | ErbB2, EGFR | 16 |
| ERL | Erlotinib | ErbB1 (EGFR) | 14 |
| GEF | Gefitinib | ErbB1 (EGFR) | 18 |
| LAP | Lapatinib | ErbB1 (EGFR), HER2 | 27 |
| VAN | Vandetanib | VEGFR, EGFR, RET | 16 |

**5 种药物，91 样本**

### Category 2: 嘧啶母核激酶抑制剂 (Pyrimidine-core KIs)

以嘧啶/氨基嘧啶为母核，靶向 BCR-ABL、JAK、ALK、MEK、B-Raf、VEGFR 等。

| 代码 | 名称 | 靶点 | 样本数 |
|------|------|------|--------|
| CER | Ceritinib | ALK | 16 |
| DAB | Dabrafenib | BRAF | 16 |
| DAS | Dasatinib | ABL, SRC, KIT | 16 |
| IMA | Imatinib | BCR-ABL | 20 |
| NIL | Nilotinib | BCR-ABL | 16 |
| PAZ | Pazopanib | VEGFR, PDGFR | 16 |
| PON | Ponatinib | BCR-ABL, SRC | 16 |
| RUX | Ruxolitinib | JAK | 16 |
| TOF | Tofacitinib | JAK | 18 |
| TRA | Trametinib | MEK1, MEK2 | 16 |

**10 种药物，166 样本**

### Category 3: 其他杂环激酶抑制剂 (Other heterocyclic KIs)

母核多样（二苯脲、吡啶、吲哚、吲唑、喹啉、噻唑/吡啶并嘧啶、磺酰胺等），靶向 VEGFR/Raf/c-Met/PDGFR 等多靶点。

| 代码 | 名称 | 靶点 | 样本数 |
|------|------|------|--------|
| AXI | Axitinib | VEGFR, PDGFR | 16 |
| BOS | Bosutinib | BCR-ABL, SRC | 16 |
| CAB | Cabozantinib | c-Met, VEGFR2 | 14 |
| CRI | Crizotinib | ALK, HGFR (c-Met) | 19 |
| REG | Regorafenib | VEGFR, PDGFR, RET | 20 |
| SOR | Sorafenib | BRAF, VEGFR, PDGFR | 28 |
| SUN | Sunitinib | VEGFR, PDGFR, KIT | 22 |
| VEM | Vemurafenib | BRAF V600E | 16 |

**8 种药物，151 样本**

### Category 4: 类固醇及相关 (Steroid-related)

| 代码 | 名称 | 类别 | 样本数 |
|------|------|------|--------|
| EST | Estradiol | 雌激素 (steroid hormone) | 18 |
| URS | Ursodiol | 胆汁酸 (bile acid, steroid-related) | 12 |

**2 种药物，30 样本**

### Category 5: 神经系统药物 (Neuro/psychoactive drugs)

| 代码 | 名称 | 类别 | 样本数 |
|------|------|------|--------|
| DOM | Domperidone | 多巴胺 D2 拮抗剂 (benzimidazole) | 3 |
| DTP | Diethylpropion | 中枢兴奋剂 (phenethylamine) | 3 |
| FLE | Flecainide | 钠通道阻滞剂 (benzamide) | 9 |
| LOP | Loperamide | 阿片受体激动剂 (piperidine) | 4 |
| PAR | Paroxetine | SSRI (benzodioxole) | 7 |

**5 种药物，26 样本**

### Category 6: 其他治疗药物 (Other diverse therapeutics)

| 代码 | 名称 | 类别 | 样本数 |
|------|------|------|--------|
| ALE | Alendronic acid | 双膦酸盐 (bisphosphonate) | 3 |
| AMI | Amiodarone | 抗心律失常 (benzofuran) | 8 |
| CYC | Cyclosporine | 免疫抑制 (cyclic peptide) | 17 |
| ENT | Entecavir | 抗病毒 (nucleoside analog) | 3 |
| MIL | Milrinone | PDE3 抑制剂 (pyridine-carbonitrile) | 8 |
| OLM | Olmesartan | 血管紧张素 II 受体拮抗剂 (biphenyl-tetrazole) | 18 |

**6 种药物，57 样本**

### 总计

| 分类 | 药物数 | 样本数 |
|------|--------|--------|
| Category 1: 喹唑啉类 KI | 5 | 91 |
| Category 2: 嘧啶母核 KI | 10 | 166 |
| Category 3: 其他杂环 KI | 8 | 151 |
| Category 4: 类固醇及相关 | 2 | 30 |
| Category 5: 神经系统药物 | 5 | 26 |
| Category 6: 其他治疗药物 | 6 | 57 |
| **合计** | **36** | **521** |

注：521 样本仅统计小分子药物样本，不包含抗体（19 孔 TRS）、蛋白（30 孔 IGF+TNF）及 29 孔无元数据匹配的扰动。

## 数据处理流程

DTOXS 原始 Level3 read count 数据经以下步骤处理为 XPert 模型的微调输入：

```
Raw Level3 TSV (read counts per gene per well)
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│ Step 1: Gene name mapping                               │
│  DTOXS gene symbols → L1000 978 landmark genes          │
│  Using UCSC-hg38-Gene-HGNC-Entrez-Map.tsv               │
│  27 genes need special mapping (18 renamed, 9 direct)   │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│ Step 2: Normalization                                   │
│  Compute: log2(CPM + 1) per well                        │
│  Output: 978 genes x 956 wells expression matrix        │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│ Step 3: Treatment-Control pairing                       │
│  Pair each treatment well with its matched control well │
│  x_pert = treatment expression (978 genes)              │
│  x_ctl  = control expression (978 genes)                │
│  Control states: CTRL, CTRL-DMSO, CTRL-EtOH             │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│ Step 4: Compute DEG                                    │
│  deg = x_pert - x_ctl  (log2 fold change, 978 genes)   │
│  Output: xpert/finetune/data/dtoxs.h5ad                 │
└─────────────────────────────────────────────────────────┘
```

处理的执行代码位于 `xpert/` 项目中。`dtoxs/` 目录仅包含原始数据和文档，无可执行代码。

输出文件 `dtoxs.h5ad` 为 AnnData 格式，包含：
- `X`: 药物扰动 DEG 矩阵（treatment - control, log2 fold change）
- 相关元数据：药物信息、细胞供体、实验批次等

## 与 XPert 模型的关系

XPert 是一个双分支 Transformer 模型（约 25M 参数），用于从小分子化学结构预测其在目标细胞类型上的 978 个 landmark 基因表达扰动。

| 阶段 | 数据集 | 数据空间 | 用途 |
|------|--------|----------|------|
| 预训练 (Pretrain) | L1000 (LINCS Phase I/II) | log2 Q2NORM | 通用转录组扰动预测 |
| 微调 (Finetune) | DToxS (LDS-1593) | log2 RMA | 心肌细胞特异性适配 |

两个数据集均在 log2 表达空间，无需 z-score 转换。

### 微调策略

DToxS 支持五种交叉验证策略：

1. **Warm（温热微调）**: 7-fold drug-stratified 交叉验证。按药物分层划分训练/验证集，同一药物可同时出现在训练和验证集中（不同样本）。含验证集用于早停。

2. **Warm-start（温热评价微调）**: 10-fold drug-stratified 交叉验证。无验证集（与 cold-* 一致），同一药物可同时出现在训练和测试集中。**仅用于模型评价**，不用于最终预测。每个 DTOXS 样本至少在一个 fold 的测试集中出现。

3. **Cold-drug（按药物冷启动）**: 10-fold leave-one-drug-out。每种药物整体放入一个 fold，验证时模型未曾见过该药物的任何样本。

4. **Cold-cate（按类别冷启动）**: 6-fold leave-one-category-out。36 种药物按化学母核分为 6 类，每次留出一整类药物作为验证集，测试模型对未见结构类别的泛化能力。

5. **DTOXS（全量训练）**: 单 fold，所有样本作为训练集，无验证/测试集。**用于最终预测模型训练**，不参与模型评价。

新增的 warm-start_fold-1..10 和 dtoxs_full 划分列由 `xpert/finetune/create_dtoxs_splits.py` 脚本生成。

## 关键注意事项

### 基因符号映射

DToxS 数据集的基因符号与 L1000 landmark 基因因不同项目的基因符号更新时间差而存在差异。映射涉及 **27 个基因**：
- **18 个**为基因重命名（如 ELP1 原称 IKBKAP, COQ8A 原称 ADCK3）
- **9 个**为直接匹配（符号一致但需要确认）

映射依据 `Metadata/UCSC-hg38-Gene-HGNC-Entrez-Map.tsv` 文件。在实际处理中，需以该映射文件为准进行基因符号转换，否则会导致约 2.8% 的 landmark 基因丢失。

### 数据显示空间

- DTOXS RMA 归一化后的 log2 表达值为**非负数**，典型范围 [0, 15]
- 与 L1000 Q2NORM（分位数归一化 log2）表达空间一致，两者可直接互操作
- 此为 log2 表达值，**非 z-score**，不得额外进行标准化

### 药物 SMILES 可用性

- **672 个**小分子给药孔具有有效的 canonical SMILES（来自 `Small_Molecule_Metadata.txt` 的 `sm_smiles_parent` 或 `sm_smiles_batch` 字段）
- **49 个**孔为抗体或重组蛋白扰动（TRS, IGF, TNF），无 SMILES，仅用于下游 MOA 分析
- **29 个**孔无任何元数据匹配（`END`, `SUN+`, `DIC`），在大部分分析中被排除

### 实验批次效应

5 个测序实验（Exp 1, 2, 4, 5）和 Exp 5 内部的 5 个子实验可能存在批次效应。数据分析时建议将实验批次作为协变量考虑。

### 对照设定

三种对照状态（CTRL, CTRL-DMSO, CTRL-EtOH）均视为对照类，共 206 孔。对照样本分布在各子实验中。XPert 处理管线将给药孔与其对应批次的对照孔配对计算 DEG。

### 数据文件格式

- Level3 目录下每个子实验包含两个配套 TSV 文件：
  - `PromoCell-DGE-Experiment-Design*.tsv`：实验设计表，每行一个孔
  - `PromoCell-DGE-Read-Counts*.tsv`：read count 矩阵，列为孔位 Well
- 总共 200 个 TSV 文件（设计表 + 读数表各约 100 个）
- 孔位 Well 是连接设计表和读数表的主键
