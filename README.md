# arc — AI-Driven MIRI Drug Discovery Pipeline

利用人工智能与生物信息学方法，从小分子化学结构出发，预测可用于治疗**心肌缺血再灌注损伤（MIRI）**的药物。

## Pipeline 总览

```
┌─────────────────────────────────────────────────────────────┐
│                     Drug SMILES (小分子结构)                   │
└─────────────────────────────────────────────────────────────┘
                              │
                    ┌─────────┴─────────┐
                    ▼                   ▼
             ┌──────────────┐   ┌──────────────┐
             │   HG (hg/)   │   │    UniMol    │
             │ 异质图GNN嵌入  │   │  3D分子构象特征 │
             │  InfoNCE对比  │   │  unimol_tools│
             └──────────────┘   └──────────────┘
                    │                   │
                    └─────────┬─────────┘
                              │
                    药物嵌入向量拼接 (Drug Embeddings)
                    11104 药物 × (HG + UniMol)
                              │
                              ▼
                    ┌──────────────────┐
                    │  XPert (xpert/)  │
                    │  双分支Transformer │
                    │  ~25M parameters │
                    │                  │
                    │  Input: 药物嵌入  │
                    │       + 细胞状态  │
                    │  Output: 978基因 │
                    │     转录组扰动预测 │
                    └──────────────────┘
                              │
                    3127 FDA 药物 × 978 基因
                    差异表达预测 (DEG)
                              │
                              ▼
                    ┌──────────────────┐
                    │   fda_miri_vs/ (药物排名)  │
                    │                  │
                    │  MIRI 疗效评分    │
                    │  ├─ 表型基因集    │
                    │  └─ 25条通路     │
                    │                  │
                    │  心脏毒性评分     │
                    │  ├─ 11类药物签名  │
                    │  └─ 25条通路     │
                    │                  │
                    │  整合排名 + 文献验证│
                    └──────────────────┘
                              │
                              ▼
                    fda-miri_score.csv
                    (3127 药物最终排名)
```

## 子项目索引

| 目录 | 角色 | 核心技术 | 状态 | 文档 |
|------|------|----------|:----:|------|
| `xpert/` | 药物转录组扰动预测模型 | 双分支 Transformer (~25M params) | ✅ | [CLAUDE.md](xpert/CLAUDE.md) / [README.md](xpert/README.md) |
| `dtoxs/` | 心肌细胞药物扰动组数据集 | DToxS LINCS LDS-1593 (914样本) | ✅ | [CLAUDE.md](dtoxs/CLAUDE.md) / [README.md](dtoxs/README.md) |
| `hg/` | 异质图神经网络 | HeteroGNN + InfoNCE 对比学习 | ✅ | [CLAUDE.md](hg/CLAUDE.md) / [README.md](hg/README.md) |
| `fda_miri_vs/` | 药物排名 + 临床验证 | GSEA/RGES + PRISMA 文献审查 | ✅ | [CLAUDE.md](fda_miri_vs/CLAUDE.md) / [README.md](fda_miri_vs/README.md) |
| `references/` | 参考文献存档 | PDF | — | [CLAUDE.md](references/CLAUDE.md) / [README.md](references/README.md) |

### 子项目间数据流

| 上游 | 下游 | 数据文件 | 说明 |
|------|------|----------|------|
| hg/ | xpert/ | `all_drugs_idx2HG.npy` | 生物嵌入 (11104×256) |
| — | xpert/ | `all_drugs_unimol_arr.npy` | UniMol 3D 特征 (11104×122×514) |
| xpert/ | fda_miri_vs/ | `xpert_cm-fda_deg.csv` | DEG 预测 (3127×978 genes) |

特征以预计算 `.npy` 文件交换，无共享运行时依赖。全部子项目使用 conda env `xpert` (Python 3.13)。

## 核心方法学

### 1. 药物特征提取

两种互补的药物表示被拼接为 XPert 的输入：

- **HG（生物学嵌入）**：在药物-基因-疾病异质图上通过 InfoNCE 对比学习训练 HeteroGNN，生成反映药物生物学功能的嵌入向量
- **UniMol（3D 构象特征）**：基于分子 3D 构象的预训练模型，捕获三维空间结构信息

### 2. 转录组扰动预测 (XPert)

双分支 Transformer 架构：
- **Base Encoder**（4 层 Self-Attention）：编码未扰动细胞状态和基因-基因互作关系
- **Pert Encoder**（CA+SA+SA+CA）：通过 Cross-Attention 将药物信息与细胞上下文融合
- **输出**：扰动表达 (x_pert)、对照表达 (x_ctl)、差异表达 DEG (x_deg = x_pert − x_ctl)
- **Loss**：0.2·L_pert + 0.003·L_ctl + 0.2·L_deg + 1.0·L_pcc

模型在 L1000 数据集（8,981 种药物，978 个 landmark 基因）上预训练，在 DToxS 心肌细胞数据集（914 样本）上微调。

### 3. 药物排名 (fda_miri_vs/)

采用 **GSEA (Gene Set Enrichment Analysis)** 方法学，分为四个评分模块：

**MIRI 疗效评分**（越高越好）：
- **表型评分 (E_pheno)**：基于 MIRI 疾病差异表达基因（41 上调 + 28 下调）计算 RGES 逆转评分
- **通路评分 (E_pathway)**：25 条 MIRI 相关通路的富集分数均值

**心脏毒性评分**（越高越毒）：
- **药物类别评分 (T_drug)**：与 11 类已知心脏毒性药物的转录组签名进行 RGES 相似度比对
- **通路评分 (T_pathway)**：25 条心脏毒性相关通路的富集分数均值

**最终整合**（v3.0 统一 [0,1] 架构）：
- `Score = (Efficacy + (1 − Toxicity)) / 2`
- `ParetoRank`：在 (Efficacy, 1−Toxicity) 上的非支配排序

### 4. 心脏毒性参考药物类别

| # | 类别 | 参考药物 | 数据来源 |
|:--:|------|----------|----------|
| 1 | Anthracycline | Doxorubicin | GEO GSE181517 |
| 2 | Anti-HER2 TKI | Lapatinib | DToxS |
| 3 | BCR-ABL TKI | Ponatinib | DToxS |
| 4 | VEGF/Multi-KI | Sunitinib | DToxS |
| 5 | Class III AAD | Amiodarone | DToxS |
| 6 | Class I AAD | Flecainide | DToxS |
| 7 | Alkylating Agent | Cyclophosphamide | GEO GSE194073 |
| 8 | Antimetabolite | 5-Fluorouracil | GEO GSE166957 |
| 9 | BTK Inhibitor | Ibrutinib | GEO GSE128688 |
| 10 | Dopamine Antagonist | Domperidone | DToxS |
| 11 | Proteasome Inhibitor | Bortezomib | XPert predicted |

## 数据空间说明

- **L1000**：LINCS Phase I/II, Level 3 Q2NORM（分位数归一化 log2 表达），978 landmark genes，非负值
- **DTOXS**：人源心肌细胞 3'-DGE RNA-seq，RMA 归一化 log2 表达，非负值，范围 [0, ~15]
- **关键兼容性**：预训练和微调数据都在 **log2 表达空间** — 无需 z-score 转换
- **药物嵌入**：覆盖 11,104 药物（indices 0–11103）。前 8,981 个用于 L1000 预训练；8,981–11,103 为 FDA 扩展药物。HG 嵌入：训练行 = 作者参考，FDA 行 = HG 微调 + 分位数范数映射 (FDA/Train std=1.000)。UniMol 特征：全部 11,104 药物通过标准 UniMolRepr 生成

## 当前进度 (2026-06-10)

| 阶段 | 状态 | 关键输出 |
|------|:----:|----------|
| HG — 异质图预训练 + 嵌入生成 | ✅ | 11,104×256 嵌入 (微调+范数映射, FDA/Train std=1.000) |
| XPert — L1000 预训练 | ✅ | 8,981 药物模型权重 |
| XPert — DTOXS 微调 (cold-drug/cold-cate) | ✅ | 16 折模型权重 |
| XPert — DTOXS 全量训练划分 (dtoxs) | ✅ | create_dtoxs_splits.py, dtoxs.yaml |
| XPert — FDA 药物 DEG 预测 | ✅ | 3127×978 DEG 矩阵 |
| Drug Ranking — MIRI 疗效评分 (phenotype + pathway) | ✅ | E_pheno + E_pathway |
| Drug Ranking — 心脏毒性评分 (drug-class + pathway) | ✅ | T_drug + T_pathway |
| Drug Ranking — 最终整合排名 | ✅ | `fda-miri_score.csv` |
| GPT Researcher — PRISMA 系统文献审查 | ✅ | Top candidates evidence report |
| 候选药物实验验证优先级排序 | 🔜 | — |

### XPert 模型评估

模型在 DTOXS 数据集上通过 4 组模型进行评价（`evaluate_dtoxs.py`）：
- **L1000**: 预训练模型在 warm-start 测试集上的表现
- **DToxS warm-start**: warm-start 微调模型 (10 folds)
- **DToxS cold-drug**: cold-drug 微调模型 (10 folds)
- **DToxS cold-cate**: cold-cate 微调模型 (6 folds)

每组模型报告 14 个 DEG-only 指标（MSE, RMSE, MAE, Pearson, Spearman, R2, Precision@K Positive/Negative for K=10,20,50,100），结果以 CSV + PNG 柱状图保存至 `metrics/dtoxs/`。

当前 XPert 预测结果药物间方差极低（97.2% 为下调 DEG），影响下游药效评分的区分能力。

详见 [xpert/CLAUDE.md](xpert/CLAUDE.md) 和 [xpert/README.md](xpert/README.md)。

## 快速开始

### 预测 FDA 药物转录组扰动

```bash
cd xpert/
bash predict/predict.sh \
    --prediction_root predict/fda-cm \
    --ctl_csv cm_ctl_landmark.csv \
    --drug_csv fda_drugs.csv \
    --model_path finetune/output/dtoxs_dtoxs/
```

### 运行药物排名

```bash
cd fda_miri_vs/
python run_ranking.py                    # 主排名 pipeline
python plots/aggregate_and_plot.py       # 评分聚合 + 可视化
```

### 训练模型

```bash
# XPert 预训练
cd xpert/ && bash pretrain/run.sh [dataset] [split_method] [config]

# XPert 微调
cd xpert/ && bash finetune/run.sh <pretrained_path> [dataset] [split_method] [config]

# HG 训练 (从头)
cd hg/ && bash run.sh --device cuda:0 --epochs 500

# HG FDA 药物嵌入生成 (微调 + 范数映射)
cd hg/ && python -m src.finetune_hg --device cuda:0 --lr 5e-5 --epochs 50 --update-xpert-embeddings
```

## 文档约定

每个子目录包含两个文档文件：
- **CLAUDE.md** — LLM 智能体上下文索引（英文，无表格，简洁）
- **README.md** — 人类可读项目说明书（中文叙述，含表格和流程图）

详见根目录 [CLAUDE.md](CLAUDE.md) 中的完整约定说明。

### Agent Skills 约定

仓库已接入工程技能套件（Matt Pocock's Skills）。约定如下：

- **Issue tracker** — 本地 Markdown，issue/spec 存放于 `.scratch/<feature>/`。见 [docs/agents/issue-tracker.md](docs/agents/issue-tracker.md)
- **Triage labels** — 五个标准角色标签：`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`。见 [docs/agents/triage-labels.md](docs/agents/triage-labels.md)
- **Domain docs** — 单上下文布局：根目录 `CONTEXT.md` + `docs/adr/`。见 [docs/agents/domain.md](docs/agents/domain.md)
