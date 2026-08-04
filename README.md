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
                    ┌──────────────────────────┐
                    │  miri_reposition/xpert-fda-v1/ │
                    │  (药物排名工作流, stages 01-10) │
                    │                              │
                    │  MIRI 疗效评分 (3 个评分器)    │
                    │  ├─ WCS + RGES               │
                    │  ├─ Zhang / ssCMap           │
                    │  └─ Pathway RGES             │
                    │                              │
                    │  4 个 Z 分量 → 整合排名         │
                    │  过滤 (WCS<0 + 关键通路)        │
                    └──────────────────────────────┘
                              │
                              ▼
                    data/ranked/ranked_drugs.csv
                    (最终排名药物列表)
 ```

## 子项目索引

| 目录 | 角色 | 核心技术 | 状态 | 文档 |
|------|------|----------|:----:|------|
| `pert_model/xpert/` | 药物转录组扰动预测模型（库） | 双分支 Transformer (~25M params) | ✅ | [CLAUDE.md](pert_model/xpert/CLAUDE.md) / [README.md](pert_model/xpert/README.md) |
| `pert_data/dtoxs/` | 心肌细胞药物扰动组数据集 | DToxS LINCS LDS-1593 (914样本) | ✅ | [CLAUDE.md](pert_data/dtoxs/CLAUDE.md) / [README.md](pert_data/dtoxs/README.md) |
| `drug_emb_model/hg/` | 异质图神经网络 | HeteroGNN + InfoNCE 对比学习 | ✅ | [CLAUDE.md](drug_emb_model/hg/CLAUDE.md) / [README.md](drug_emb_model/hg/README.md) |
| `miri_reposition/xpert-fda-v1/` | 药物排名工作流（活跃） | WCS/RGES + Zhang + Pathway 评分 → 整合排名 | ✅ | [CLAUDE.md](miri_reposition/xpert-fda-v1/CLAUDE.md) |
| `miri_reposition/xpert-fda-v1_history/` | 冻结的旧版 fda_miri_vs 快照 | 历史参考（gitignored） | — | — |
| `disease_data/miri/` | MIRI 疾病数据（per-GSE） | GEO 数据 + DEG 列表 | ✅ | [CLAUDE.md](disease_data/miri/CLAUDE.md) |
| `reference/` | 参考文献存档 | PDF | — | [CLAUDE.md](reference/CLAUDE.md) / [README.md](reference/README.md) |

### 子项目间数据流

| 上游 | 下游 | 数据文件 | 说明 |
|------|------|----------|------|
| drug_emb_model/hg/ | xpert-fda-v1 | `data/embeddings/all_drugs_idx2HG.npy` | 生物嵌入 (11104×256)，stage 02 生成 |
| — | xpert-fda-v1 | `data/embeddings/all_drugs_unimol_arr.npy` | UniMol 3D 特征 (11104×122×514)，冻结输入 |
| pert_model/xpert/ | xpert-fda-v1 | `data/deg/xpert_fda-cm_deg.csv` | DEG 预测 (3127×978 genes)，stage 05 |
| disease_data/miri/ | xpert-fda-v1 | `disease_data/miri/{GSE}/deg/` | per-GSE DEG 列表 → stage 07 参考构建 |

特征以预计算 `.npy` 文件交换，无共享运行时依赖。全部子项目使用 conda env `arc` (Python 3.13)，由根目录 `environment.yml` 创建。

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

### 3. 药物排名 (miri_reposition/xpert-fda-v1/)

工作流为编号阶段 01-10，由 `screen.ipynb` 编排。评分采用三种互补方法（沿用冻结旧管线的逻辑，已逐字节验证一致）：

- **WCS + RGES**（基因集）：Weighted Connectivity Score + Reverse Gene Expression Score，基于 MIRI 上下调基因集
- **Zhang Score / ssCMap**（基于排名）：签名基因排名乘积，保留疾病效应幅度
- **Pathway RGES**（机制）：40 条 MIRI 通路（20 UP + 20 DOWN，8 条关键）的逐通路 RGES

**整合**：`integrated = −0.3·RGES_zs − 0.2·Zhang_zs + 0.3·pathway_reversal_zs + 0.025·Σcritical_efficacy_zs`，然后过滤（WCS<0 + ≥6/8 关键通路）并排序。心脏毒性评分从 v1 推迟（见 `docs/adr/0004`）。

### 4. 心脏毒性评分

心脏毒性评分从 v1 推迟（当前代码库中不存在该评分器），不设占位阶段。参见 `docs/adr/0004`。

## 数据空间说明

- **L1000**：LINCS Phase I/II, Level 3 Q2NORM（分位数归一化 log2 表达），978 landmark genes，非负值
- **DTOXS**：人源心肌细胞 3'-DGE RNA-seq，RMA 归一化 log2 表达，非负值，范围 [0, ~15]
- **关键兼容性**：预训练和微调数据都在 **log2 表达空间** — 训练无需 z-score 转换；但评分前必须对预测 DEG 做 per-gene z-score（工作流 stage 06）
- **药物嵌入**：覆盖 11,104 药物（indices 0–11103）。前 8,981 个用于 L1000 预训练；8,981–11,103 为 FDA 扩展药物。HG 嵌入：训练行 = 作者参考，FDA 行 = HG 微调 + 分位数范数映射 (FDA/Train std=1.000)。UniMol 特征：全部 11,104 药物通过标准 UniMolRepr 生成

## 当前进度 (2026-08-04)

| 阶段 | 状态 | 关键输出 |
|------|:----:|----------|
| 工作区重构为分类目录结构 (ADR 0001) | ✅ | pert_data/ drug_emb_model/ pert_model/ miri_reposition/ disease_data/ |
| 组件即库 (ADR 0002) | ✅ | 组件只产出代码+权重；工作流脚本生成所有产物 |
| conda env `arc` 规格 (ADR 0003) | ✅ | 根目录 environment.yml |
| 旧 fda_miri_vs 冻结 + 重建 (ADR 0004) | ✅ | xpert-fda-v1_history/ (gitignored) + xpert-fda-v1/ |
| 疾病数据迁移 disease_data/miri/{GSE} | ✅ | 每 GSE data/ + deg/ + build_deg.R |
| xpert-fda-v1 工作流 stages 01-10 | ✅ | screen.ipynb 编排 |
| MIRI 疗效评分（3 评分器） | ✅ | 4 Z 分量 + 整合分数 |
| 评分逻辑与冻结管线一致性 | ✅ | 相同输入下逐字节一致 |
| 心脏毒性评分 | 🔜 | v1 推迟 (ADR 0004) |

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

### 运行完整筛选工作流

```bash
conda activate arc
cd miri_reposition/xpert-fda-v1/
jupyter notebook screen.ipynb      # 按顺序运行 stages 01-10
# 或逐阶段运行:
python scripts/01_drug_list.py
python scripts/06_zscore.py
python scripts/07_build_miri_deg_reference.py
python scripts/08_efficacy_scores.py
python scripts/10_filter_rank.py
```

### 训练模型

```bash
# XPert 预训练
cd pert_model/xpert/ && bash pretrain/run.sh [dataset] [split_method] [config]

# HG 训练 (从头)
cd drug_emb_model/hg/ && bash run.sh --device cuda:0 --epochs 500

# HG FDA 药物嵌入生成 (微调 + 范数映射)
cd drug_emb_model/hg/ && python -m src.finetune_hg --device cuda:0 --lr 5e-5 --epochs 50 \
    --update-xpert-embeddings --embeddings-output embeddings/all_drugs_idx2HG.npy
```

XPert 微调、FDA 预测、baseline 构建等步骤已移入工作流
`miri_reposition/xpert-fda-v1/`（stages 03-05），详见其 CLAUDE.md。

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
