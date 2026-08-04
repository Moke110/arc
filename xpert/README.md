# XPert: 药物转录组扰动预测模型

## 项目目标

XPert 是一个双分支 Transformer 模型，用于从小分子药物的化学结构和生物特征出发，预测药物处理后在目标细胞中引起的转录组扰动效应（978 个 L1000 landmark 基因的表达变化）。模型输出三个核心量：处理组表达 (x_pert)、对照组表达 (x_ctl)、以及差异表达 (x_deg = x_pert - x_ctl)。

XPert 是本项目 `arc/` 的核心预测引擎 —— 它为下游的药物排名系统（`drug_ranking/`）提供 FDA 药物的预测 DEG 矩阵，从而支持心肌缺血再灌注损伤 (MIRI) 的候选药物筛选。

## 工作流位置

XPert 在整个项目 pipeline 中处于**中心桥梁**位置：

```
Drug SMILES
    │
    ├──→ UniMol                → 3D 分子特征 (atom-level features, 122×512)
    └──→ HeteroGNN (hg/)      → 生物学嵌入 (biological embedding, 256-dim)
            │
            ▼
        XPert (xpert/)         ← 本模块
            │
            ▼
    预测 DEG 矩阵 (N_drugs × 978 genes)
            │
            ▼
    Drug Ranking (drug_ranking/) → MIRI 疗效评分 + 心脏毒性评分
```

## 目录结构

| 目录 / 文件 | 说明 |
|---|---|
| `src/train.py` | CLI 入口：参数解析、训练/验证/推理/测试循环、主函数 |
| `src/model_xpert.py` | 核心模型：`XPertNet`（主模型）+ `AttnEncoder`（双分支 Transformer 堆叠） |
| `src/model_utils.py` | 基础模块：自注意力 (SA)、交叉注意力 (CA)、各类 Embedding、Encoder/Decoder、FFN |
| `src/dataset.py` | 数据加载：`MyDataset` 将 AnnData 转为 PyTorch 张量（惰性加载，预提取 numpy 数组） |
| `src/utils.py` | 工具函数：DataLoader 构建、损失函数 (MSE/PCC/Cosine)、早停、UniMol 特征生成 |
| `src/config.py` | 配置系统：YAML 加载、相对路径解析、输出目录规范 |
| `src/metrics.py` | 评估指标：MSE、PCC、Precision@K，支持 raw 和 deg 两种变体 |
| `src/eval_metrics.py` | 扩展评估：R²、MMD、Wasserstein 距离 |
| `src/predict.py` | FDA 预测模块：无 DataLoader 依赖，直接构建 10-元组批次，集成推理 |
| `src/vis_gene_attn.py` | 基因注意力可视化：热图、箱线图、散点图 |
| `src/vis_molecule.py` | 分子注意力可视化：基于 RDKit 渲染 |
| `pretrain/` | 阶段一：在 L1000 大规模数据集上预训练 |
| `pretrain/config/l1000.yaml` | 预训练配置：batch=64, lr=0.004, hidden=256 |
| `pretrain/run.sh` | 预训练启动脚本（环境变量驱动） |
| `finetune/` | 阶段二：在目标数据集上微调 |
| `finetune/config/dtoxs.yaml` | DTOXS 心肌细胞数据集微调配置，含五种数据划分方案 |
| `finetune/create_dtoxs_splits.py` | 生成 warm-start 和 dtoxs_full 数据划分列 |
| `finetune/run.sh` | 微调启动脚本 |
| `predict/` | 阶段三：对任意药物（如 FDA 药物库）进行扰动预测 |
| `predict/predict.sh` | FDA 预测 CLI 封装 |
| `predict/build_baseline.py` | 从 DTOXS 原始数据构建心肌细胞共识基线 |
| `metrics/` | 评估脚本和结果 |
| `metrics/evaluate_dtoxs.py` | 标准 DTOXS 评估脚本 |
| `metrics/evaluate_dtoxs_zscore.py` | Z-score 后处理评估：4 组模型 (deg vs deg-zscore) 双栏比较 |
| `drug_embeddings/` | 共享药物特征文件（11,104 种药物，索引 0-11103） |
| `metadata/` | 基因/药物参考信息（L1000 基因符号、Entrez ID） |

### drug_embeddings/ 详细说明

| 文件 | 形状 | 来源 | 说明 |
|---|---|---|---|
| `all_drugs_unimol_arr.npy` | 11104×122×514 | UniMol (unimol_tools) | 3D 分子构象特征，全部 11,104 药物标准 UniMolRepr 生成 |
| `all_drugs_idx2HG.npy` | 11104×256 | hg/ | HG 生物嵌入。训练行 = 作者参考，FDA 行 = HG 微调 + 分位数范数映射 (FDA/Train std=1.000) |
| `all_drugs_idx2smi.npy` | 11104 (dict) | ChEMBL | 药物索引到 SMILES 的映射 (索引 0–11103) |
| `PPI_gene_vector_128d.npy` | 978×128 | STRING | PPI 网络的基因位置编码 |

## 模型架构

### 整体架构图

```
                          XPertNet Forward Pass
                          ====================

  Drug Side                                     Cell Side
  ---------                                     ---------
  UniMol Feat (122×512)                        Control Expr (978-dim)
      │                                              │
      ▼                                              ▼
  Linear(512→256)                              Digitize into n_bins
      │                                              │
      ├─── HG Embed (256-dim) ← drug_idx lookup      ▼
      │    (replaces token-0)                  cell_Embeddings
      │                                       (bin embed + PPI position embed)
      ├─── Dose Embed (learnable)                     │
      ├─── Time Embed (learnable)                     │
      │                                               │
      ▼                                               ▼
  unimol_Embeddings                            Cell Tokens (978×256)
  (projection + learnable pos encoding)              │
      │                                               │
      ▼                                               │
  Drug Tokens (N_atoms×256)              ┌────────────┘
      │                                   │
      │         ┌─────────────────────────┘
      │         │
      ▼         ▼
  ┌──────────────────────────────────────────────┐
  │            Dual-Branch Transformer           │
  │                                              │
  │   Pert Branch                   Ctl Branch   │
  │   CA + SA + SA + CA           SA + SA + SA + SA │
  │   (cross-attn drug→cell)     (self-attn cell only) │
  │                                              │
  └──────────────────────────────────────────────┘
      │                                   │
      ▼                                   ▼
  trt_fc(256→64→1)                   ctl_fc(256→64→1)
      │                                   │
      ▼                                   ▼
  x_pert (978-dim)                   x_ctl (978-dim)
      │                                   │
      └────────────┬──────────────────────┘
                   │
                   ▼
         x_deg = x_pert - x_ctl
              (978-dim)
```

### 关键设计说明

**药物侧 (Drug Side)**：
- UniMol 原子特征：每个药物分子的原子级 3D 特征（最多 122 个原子，每个原子 512 维），通过线性层投影到 256 维
- HG 生物嵌入：HeteroGNN 学习得到的药物全局特征（256 维），替换 UniMol 序列的第一个 token
- 剂量/时间嵌入：可学习的嵌入向量，编码药物处理的时间和剂量信息
- 药物 token 序列 = [dose_token, time_token, HG_token, atom_1, ..., atom_N]

**细胞侧 (Cell Side)**：
- 表达值分箱 (Binning)：将连续的表达值离散化为 n_bins=128 个区间，每个基因获得一个离散 token
- PPI 基因位置编码：基于 STRING PPI 网络预训练的基因向量（128 维），提供基因间的拓扑关系先验
- 细胞 token = bin_embedding[expr_discrete] + PPI_position_embedding[gene_index]

**双分支设计 (Dual-Branch)**：
- Pert Branch (trt_structure = CA+SA+SA+CA)：通过交叉注意力让药物信息与细胞状态交互，模拟药物处理后细胞的状态
- Ctl Branch (ctl_structure = SA+SA+SA+SA)：仅基于细胞自注意力，建模未处理细胞的基线表达
- 为什么需要两个分支？Ctl branch 提供对照基线，其输出 x_ctl 与 Pert branch 输出 x_pert 共同计算 DEG，这样的分解使得模型能分别学习基线表达模式和处理诱导的变化

**输出头 (Output Heads)**：
- 三个独立 MLP (256→64→1)，分别预测 x_pert、x_ctl、x_deg
- x_deg 的 MLP 输入为 (trt_cell_embed - ctl_cell_embed)，显式建模差异
- 这是一种多任务设计：模型同时学习对照表达、处理表达和差异表达，互相约束，提高预测精度

### 注意力机制细节

```
crossEncoder (每个 CA 层):
  drug: Self-Attention → drug_SA_embed
  cell: Self-Attention → cell_SA_embed
  cell: Cross-Attention(Q=cell_SA_embed, K=drug_SA_embed, V=drug_SA_embed) → cell_CA_embed
  cell: Feed-Forward → cell_output
  返回 (cell_output, drug_SA_embed)

Encoder (每个 SA 层):
  input: Self-Attention → attention_output
  residual + Feed-Forward → layer_output
  返回 (layer_output, attention_weights)
```

Pert branch (CA+SA+SA+CA) 中，drug-specific gene embedding 会在每个 CA 层之后注入 cell embedding（残差连接），用于特定模式下的基因-药物交互增强。

## 配置系统

### 配置文件结构

配置文件按阶段组织：

| 阶段 | 配置文件 | 关键差异 |
|---|---|---|
| 预训练 (Pretrain) | `pretrain/config/l1000.yaml` | max_value=15.002, min_value=0.0, 剂量=10类, 时间=3类, 细胞类型=240种 |
| 微调-心肌 (DTOXS) | `finetune/config/dtoxs.yaml` | max_value=15.993, min_value=1.854, 剂量=1类, 时间=3类, 细胞类型=1种 |
| 微调-临床 (CDSDB) | `finetune/config/cdsdb.yaml` | 32 batch, 特定临床数据集参数 |
| 微调-药效 (PANACEA) | `finetune/config/panacea.yaml` | 32 batch, 特定药物筛选数据集参数 |

### YAML 配置核心字段

| 字段路径 | 含义 | 预训练默认值 | DTOXS 默认值 |
|---|---|---|---|
| `model.ATTN.hidden_size` | Transformer 隐藏层维度 | 256 | 256 |
| `model.ATTN.n_heads` | 注意力头数 | 8 | 8 |
| `model.ATTN.trt_structure` | Pert branch 层结构 | CA+SA+SA+CA | CA+SA+SA+CA |
| `model.ATTN.ctl_structure` | Ctl branch 层结构 | SA+SA+SA+SA | SA+SA+SA+SA |
| `model.ATTN.topk_cell` | 细胞稀疏注意力 top-k | 128 | 128 |
| `model.ATTN.topk_drug` | 药物稀疏注意力 top-k | 32 | 32 |
| `model.ATTN.ppi_gene_vector_path` | PPI 基因向量路径 | drug_embeddings/PPI_gene_vector_128d.npy | 同上 |
| `dataset.n_bins` | 表达值分箱数 | 128 | 128 |
| `dataset.gene_num` | 基因数量 | 978 | 978 |
| `dataset.atom_num` | 原子类型数 | 24 | 24 |
| `dataset.max_atom_size` | 最大原子数 | 122 | 122 |
| `dataset.max_value` | 表达值上限 | 15.002 | 15.993 |
| `dataset.min_value` | 表达值下限 | 0.0 | 1.854 |

### 路径解析机制

`src/config.py` 中的 `resolve_paths()` 会将以下相对路径转换为绝对路径：

- `dataset.*_data_root` → 相对于 `xpert/` 项目根目录
- `model.ATTN.ppi_gene_vector_path` → 相对于 `xpert/` 项目根目录
- `drug_embeddings.*` (unimol, smi, HG) → 相对于 `xpert/` 项目根目录

### 数据划分方案 (Split Methods)

DTOXS 微调支持五种数据划分，分别在 h5ad 的 `.obs` 中存储为列：

| 方案 | 折数 | 划分策略 | 应用场景 |
|---|---|---|---|
| `warm-start` | 10 | 药物分层随机划分，无验证集，每个样本至少在一个 fold 的测试集中出现 | **仅用于模型评价** |
| `cold-drug` | 10 | 每次留出一个药物类别的所有样本作为测试集 | 评估模型对新药物的泛化能力（零样本药物预测） |
| `cold-cate` | 6 | 每次留出一个药物类别的所有样本作为测试集（更大类别） | 更严格的药物类别级泛化评估 |
| `dtoxs` | 1 | 所有样本作为训练集，无验证/测试集 | **最终预测模型**，不参与模型评价 |

新增的 `warm-start` 和 `dtoxs` 划分由 `finetune/create_dtoxs_splits.py` 脚本生成。

## 训练流程

### 损失函数

```
Loss = α·L_trt + β·L_ctl + γ·L_deg + δ·L_pcc

其中：
  L_trt = MSE(x_pert_pred, x_pert_true)       # 处理组表达重构损失
  L_ctl = MSE(x_ctl_pred, x_ctl_true)          # 对照组表达重构损失
  L_deg = MSE(x_deg_pred, x_deg_true)          # 差异表达损失
  L_pcc = (1 - PCC(x_deg_pred, x_deg_true))     # Pearson 相关系数损失

默认权重 (α, β, γ, δ) = (0.2, 0.003, 0.2, 1.0)
```

为什么 PCC 损失权重最高？PCC 衡量的是预测值与真实值在基因维度上的排序一致性，这对于下游 GSEA 分析和药物排序至关重要 —— 药物效应评分的可靠性依赖于基因相对变化方向而非绝对表达值。

### 训练阶段

```
预训练 (Pretrain):
  ┌─────────────────────────────────────────────┐
  │  数据集: L1000 (336,852 样本, 8,981 药物)    │
  │  学习率: 0.004 (前40 epoch 不变, 后降为 0.5) │
  │  批次: 64                                    │
  │  最大 Epoch: 2500                             │
  │  早停: patience=50, monitor=val PCC loss      │
  │  预热: init_epoch=70 (70 epoch 内使用         │
  │        batch-normalized loss 回传)            │
  │  输出: pretrain/output/<dataset>_<split>/    │
  └─────────────────────────────────────────────┘
                        │
                        ▼
微调 (Finetune):
  ┌─────────────────────────────────────────────┐
  │  数据集: DTOXS (914 心肌细胞样本)            │
  │  加载预训练权重 (排除特定 key 可配置)         │
  │  学习率: 0.002 (微调)                        │
  │  批次: 32-64                                 │
  │  最大 Epoch: 2-3                              │
  │  可冻结部分层 (finetune_params)               │
  │  输出: finetune/output/<dataset>_<split>/    │
  └─────────────────────────────────────────────┘
```

### 启动命令

**预训练**：
```bash
cd xpert/
# Shell 脚本方式（推荐，环境变量可覆盖参数）
bash pretrain/run.sh l1000_mdmt_full pretrain_split pretrain/config/l1000

# 直接 CLI 方式
python -m src.train \
    --mode train \
    --config pretrain/config/l1000 \
    --dataset l1000_mdmt_full \
    --device cuda:0
```

**微调**（10 epochs, 使用新 HG 嵌入）：
```bash
PRETRAINED=pretrain/output/pretrain_mdmt_full_50_epoch.pth

# DTOXS cold-drug 微调 (10-fold, leave-one-drug-out)
bash finetune/run.sh $PRETRAINED dtoxs cold-drug finetune/config/dtoxs l1000_mdmt_full

# DTOXS cold-cate 微调 (6-fold, leave-one-category-out)
bash finetune/run.sh $PRETRAINED dtoxs cold-cate finetune/config/dtoxs l1000_mdmt_full

# DTOXS warm-start 微调 (10-fold, 仅用于评价)
bash finetune/run.sh $PRETRAINED dtoxs warm-start finetune/config/dtoxs l1000_mdmt_full

# DTOXS dtoxs 微调 (1-fold, 最终预测模型)
bash finetune/run.sh $PRETRAINED dtoxs dtoxs finetune/config/dtoxs l1000_mdmt_full

# CDSDB 临床数据微调
bash finetune/run.sh \
    pretrain/output/.../best.pth \
    cdsdb random finetune/config/cdsdb l1000_mdmt_full
```

## FDA 预测流程

### 工作原理

```
  ┌─────────────────┐     ┌──────────────┐     ┌──────────────────┐
  │ FDA Drug List    │     │ Baseline      │     │ Trained Models    │
  │ (3127 drugs)     │     │ (978 genes)   │     │ (7-fold ensemble) │
  └────────┬────────┘     └──────┬───────┘     └────────┬─────────┘
           │                     │                        │
           ▼                     ▼                        ▼
  ┌─────────────────────────────────────────────────────────────┐
  │                    predict.py                                │
  │                                                              │
  │  1. SMILES → UniMol features (UniMolRepr + cache)           │
  │  2. SMILES → HG lookup index (via all_drugs_idx2smi.npy)    │
  │  3. Build batch: ctl_expr tiled N times + drug feats         │
  │  4. Ensemble: forward all models, average deg predictions    │
  │  5. Output: N_drugs × 978 genes DEG matrix                   │
  └─────────────────────────────────────────────────────────────┘
           │
           ▼
  ┌──────────────────────────────────────────────┐
  │  输出: xpert_cm-fda_deg.csv                  │
  │  - drug_name: 药物名                          │
  │  - smiles: 规范 SMILES                        │
  │  - drug_idx: 全局药物索引 (0-11103)           │
  │  - 978 个基因符号列 (各基因的预测 DEG 值)      │
  └──────────────────────────────────────────────┘
```

### 关键设计决策

1. **不依赖 DataLoader**：`predict.py` 直接构建 10-元组张量批次，不依赖 h5ad 文件和 `MyDataset`。这使得预测流程完全独立，可以在任何环境中运行。

2. **trt_raw_data = ctl_raw_data**：在预测时没有真实的处理组表达值，因此将 `trt_raw_data` 设为 `ctl_raw_data`。这不会影响前向传播结果（因为模型不从 `trt_raw_data` 获取输入信息，仅用于损失计算）。

3. **集成平均**：扫描 model_path 目录下所有 `.pth` 文件，对每个模型的 deg 预测取平均。dtoxs 全量训练集成利用了不同数据划分训练的模型的互补性。

4. **共识基线**：`cm_ctl_landmark.csv` 是 206 个 DTOXS 对照孔的均值，提供了稳定、低噪声的心肌细胞基线转录组。

### 使用方法

```bash
cd xpert/

# 预测 FDA 药物在心肌细胞上的转录组扰动
bash predict/predict.sh \
    --prediction_root predict/fda-cm \
    --ctl_csv cm_ctl_landmark.csv \
    --drug_csv fda_drugs.csv \
    --output_csv deg-perturb.csv \
    --dose_idx 0 --time_idx 2 \
    --model_path finetune/output/dtoxs_dtoxs/ \
    --device cuda:0

# dose_idx: 剂量指数 (DTOXS 固定为 0)
# time_idx: 时间指数 (0=3h, 1=6h, 2=24h，默认 24h)
```

### 中间文件

| 文件 | 说明 |
|---|---|
| `{prediction_root}/process/unimol_feats.npy` | UniMol 特征缓存（pickle 格式），首次生成后复用 |
| `{prediction_root}/process/skipped_drugs.csv` | 因无 SMILES 映射或 UniMol 生成失败而跳过的药物列表 |
| `{prediction_root}/deg-perturb.csv` | 最终预测结果：N 个药物 × (1 + 978) 列 |

## 评估指标

所有指标同时评估原始表达值 (raw) 和差异表达值 (deg)：

| 指标 | 缩写 | 含义 | 重要程度 |
|---|---|---|---|
| Mean Squared Error | MSE | 均方误差，整体预测精度 | 中 |
| Root Mean Squared Error | RMSE | 均方根误差，与原始单位一致 | 中 |
| Mean Absolute Error | MAE | 平均绝对误差，对异常值鲁棒 | 中 |
| Pearson Correlation | Pearson | 基因水平线性相关性 | **高**（下游 GSEA 依赖） |
| Spearman Correlation | Spearman | 基因水平排序相关性 | **高**（排序一致性） |
| R² Score | R² | 决定系数 | 中 |
| Precision@K Positive | P@K Pos | Top-K 上调基因的预测精度 | 高（评估方向预测） |
| Precision@K Negative | P@K Neg | Top-K 下调基因的预测精度 | 高（评估方向预测） |

只有 `_deg` 变体的指标被保留在最终评估结果中，因为差异表达才是下游应用关注的核心。

### Z-Score 后处理评估

`evaluate_dtoxs_zscore.py` 对 4 组模型（L1000 pretrain, DToxS warm-start, DToxS cold-drug, DToxS cold-cate）进行双栏比较：
- 每组两个柱子：deg（纯色填充）和 deg-zscore（同色+条纹填充）
- 误差棒：±1 std（跨 fold）
- 输出 PNG 和 CSV，保存到 `metrics/dtoxs/`

```bash
cd xpert/
python metrics/evaluate_dtoxs_zscore.py
```

## 关键设计决策

### 1. 为什么使用双分支架构？

单分支模型只能预测最终的表达状态，无法区分基线表达和药物诱导的变化。双分支设计让 Ctl branch 专门学习对照表达模式（仅依赖细胞自身特征），Pert branch 则融合药物信息来学习处理后的表达。两者的差值建模使得模型能够更准确地归因于药物处理的效果。

### 2. 为什么预训练 + 微调？

L1000 数据集包含 336,852 个样本涵盖了 8,981 种药物在多种细胞系上的扰动实验，提供了丰富的药物-细胞交互模式。预训练让模型学习通用的转录组扰动规律，然后在 DTOXS 心肌细胞数据（914 样本）上进行微调，以适配心肌细胞特定的生物学背景。

### 3. 为什么 PCC 损失权重远大于其他损失？

下游药物排名依赖于 GSEA（基因集富集分析）和 RGES（反向基因表达评分），这些方法关注的是基因的相对排序而非绝对表达值。高 PCC 权重迫使模型更关注基因间的排序一致性，直接优化了下游任务关心的指标。

### 4. 为什么使用 UniMol + HG 双路药物特征？

UniMol 提供原子级 3D 结构信息（512 维），捕捉分子的化学和构象特性。HG 生物嵌入（256 维）编码了药物在 PPI/DTI/DDS 网络中的拓扑关系，提供了生物学先验。两者的组合使得模型既能从化学结构预测药效，又能利用已知的生物交互信息。

### 5. 为什么使用分箱表达而非连续值？

将连续表达值离散化（128 bins）后使用 Embedding 表征，使得模型可以学习上下文相关的表达模式，而非简单的线性插值。这种离散化对于建模非线性基因表达关系至关重要，也让 Transformer 能够更好地捕捉基因间的交互。

### 6. 为什么 DTOXS 需要三种数据划分？

- **cold-drug** 划分验证模型对全新药物的预测能力（最严格的测试，但实际应用中重要）
- **cold-cate** 划分验证模型对全新药物类别的预测能力（评估化学空间泛化性）

三种划分提供了从宽松到严格的多层次泛化评估。

### 7. DTOXS 预测结果是以下调为主的

FDA 药物的 DEG 预测均值约 -0.38，表明大部分药物在心肌细胞中引起整体基因下调。这与激酶抑制剂类药物的已知效应一致 —— 它们通常抑制细胞周期基因和增殖相关通路，导致转录组水平下降。

### 8. POSTPROCESS_ZSCORE：预测后处理 Z-score 标准化

**动机**：XPert 预测的 DEG 值药物间方差极低（97.2% 为下调，均值 ≈ −0.38），导致基于排序的指标（Precision@K、Spearman）表现不佳，影响下游 MIRI 药物评分的区分能力。

**方法**：在推理阶段，对每个药物的预测 DEG 向量（978 个基因）进行逐行 Z-score 标准化：

```
deg_zscore[g] = (deg[g] - mean(deg)) / std(deg)   # 对每个药物单独计算
```

标准化后，每个药物扰动组内的基因相对差异被放大，改善了排序相关指标的评估。

**重要特性**：
- 仅影响推理阶段的后处理，不改变模型训练
- 保持模型参数和训练流程不变
- 仅改善排序指标（Precision@K、Spearman），不改变绝对误差指标（MSE、MAE）
- 通过 `--postprocess_zscore` 布尔标志控制，支持在 train.py 的 validate/infer 模式、predict.py 和 evaluate_dtoxs_zscore.py 中使用

**评估方式**：
- 使用 `metrics/evaluate_dtoxs_zscore.py` 对 4 组模型（pretrain, warm-start, cold-drug, cold-cate）进行 deg vs deg-zscore 双栏可视化比较
- 每组模型有两个柱子：deg（纯色填充）和 deg-zscore（同色+条纹填充）
- 重点评估 Precision@K Positive_deg、Precision@K Negative_deg 和 Spearman_deg

## 使用方法

### 环境准备

```bash
conda activate arc   # Python 3.13, PyTorch 2.0+

# 确保以下包已安装
pip install torch scanpy numpy scipy pandas scikit-learn matplotlib unimol-tools pyyaml tqdm
```

### 第零步：生成数据划分（仅首次需要）

```bash
cd xpert/
# 生成 warm-start 和 dtoxs_full 划分列到 dtoxs.h5ad
python finetune/create_dtoxs_splits.py
```

### 第一步：预训练

```bash
cd xpert/
bash pretrain/run.sh l1000_mdmt_full pretrain_split pretrain/config/l1000
```

输出位于 `pretrain/output/l1000_mdmt_full_pretrain_split/{timestamp}/`。

### 第二步：DTOXS 微调

```bash
cd xpert/

# 找到预训练的最佳模型
PRETRAINED=$(ls pretrain/output/l1000_mdmt_full_pretrain_split/*/l1000_mdmt_full_pretrain_split_fold-1_best.pth | tail -1)

# Warm 微调

# Cold-drug 微调
bash finetune/run.sh "$PRETRAINED" dtoxs cold-drug finetune/config/dtoxs l1000_mdmt_full
```

输出位于 `finetune/output/l1000_mdmt_full-dtoxs_dtoxs/{timestamp}/`。

### 第三步：FDA 药物预测

```bash
cd xpert/
bash predict/predict.sh \
    --prediction_root predict/fda-cm \
    --ctl_csv cm_ctl_landmark.csv \
    --drug_csv fda_drugs.csv \
    --output_csv xpert_cm-fda_deg.csv \
    --dose_idx 0 --time_idx 2 \
    --model_path finetune/output/l1000_mdmt_full-dtoxs_dtoxs/ \
    --device cuda:0
```

输出位于 `predict/fda-cm/xpert_cm-fda_deg.csv`。

### 第四步：结果送入药物排名

生成的 DEG 矩阵直接用于 `drug_ranking/` 模块的 MIRI 疗效和毒性评分：

```bash
cp predict/fda-cm/xpert_cm-fda_deg.csv drug_ranking/xpert/
cd drug_ranking/xpert/
python run_ranking.py
```

## 模型技术规格

| 参数 | 值 |
|---|---|
| 总参数量 | ~25M |
| 隐藏层维度 | 256 |
| 每层注意力头数 | 8 |
| Pert branch 层数 | 4 (CA+SA+SA+CA) |
| Ctl branch 层数 | 4 (SA+SA+SA+SA) |
| FFN 中间维度 | 512 (hidden * 2) |
| 输出 MLP 结构 | 256 → 64 → 1 |
| Dropout | 0.1 (attention, hidden, input) |
| 注意力实现 | F.scaled_dot_product_attention (PyTorch 2.0+) |
| 混合精度 | torch.amp GradScaler + autocast |
| 优化器 | Adam |
| LR 调度器 | LambdaLR (1.0 × 40 epochs, 0.5 after) |
| 预热机制 | init_epoch=70 (前70 epoch用batch-normalized loss) |

## 依赖关系

```
xpert/ 依赖：
├── hg/   → all_drugs_idx2HG.npy (必需)
├── UniMol → all_drugs_unimol_arr.npy (必需)
├── STING PPI → PPI_gene_vector_128d.npy (必需)
└── ChEMBL → all_drugs_idx2smi.npy (必需)

其他模块依赖 xpert/：
└── drug_ranking/ → xpert_cm-fda_deg.csv (FDA 药物 DEG 预测结果)
```
