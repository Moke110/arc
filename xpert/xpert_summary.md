# XPert 项目说明书

## 项目概述

XPert 是一个双分支 Transformer 模型，用于预测药物诱导的转录组扰动响应——即给定某个细胞在某种药物处理下，其 978 个 landmark 基因的表达变化。模型支持剂量-时间依赖的响应建模，并支持从临床前数据（L1000）到临床数据（CDSDB）的迁移学习。

- **论文**: Nature Machine Intelligence (2026)
- **模型参数量**: ~25M
- **框架**: PyTorch ≥ 2.5.1

> **历史项目文档**: 见 [.claude/xpert-history_summary.md](.claude/xpert-history_summary.md) — 原始代码仓库 file-by-file 详细总结，包含已移除的 HG 模型训练、论文图表复现等完整内容。

---

## 目录结构

```
xpert/
├── run.ipynb                  # 交互式 Notebook 启动器（3 步流程）
├── src/                        # 源代码
│   ├── train.py                # 命令行入口，参数解析，训练/验证/推理循环
│   ├── model_xpert.py          # XPertNet 主网络 + AttnEncoder
│   ├── model_utils.py          # 注意力、Embedding、Encoder 组件库
│   ├── dataset.py              # MyDataset：AnnData → PyTorch 训练样本
│   ├── utils.py                # 数据加载、损失函数、早停、UniMol 特征生成
│   ├── config.py               # YAML 配置加载与路径解析
│   ├── metrics.py              # 核心指标（MSE/PCC/Precision@K）+ get_metrics_infer
│   ├── eval_metrics.py         # 扩展指标（R²/MMD/Wasserstein）+ get_metrics_new
│   ├── vis_gene_attn.py        # 基因注意力热力图/箱线图/散点图
│   ├── vis_molecule.py         # 分子注意力可视化（RDKit 渲染）
│   └── __init__.py             # 包入口（lazy import）
├── pretrain/                   # 步骤 1：预训练
│   ├── config/l1000.yaml       # 预训练配置（batch 64, lr 0.004）
│   ├── run.sh                  # 预训练启动脚本
│   ├── data/                    # 预训练输入数据
│   └── output/                  # 预训练模型权重
├── finetune/                    # 步骤 2：微调
│   ├── config/cdsdb.yaml       # 临床 CDSDB 配置（batch 32）
│   ├── config/panacea.yaml     # PANACEA 独立验证配置
│   ├── run.sh                   # 微调启动脚本
│   ├── data/                     # 目标数据集 h5ad 文件
│   └── output/                   # 微调模型权重
├── metrics/                      # 步骤 3：评估
│   ├── evaluate.ipynb            # 评估 Notebook
│   └── output/                   # 评估结果和参考数据
├── drug_embeddings/              # 共享药物特征文件
└── metadata/                     # 基因/药物参考信息
```

---

## 一、工作流说明书

### 运行前准备

#### 环境

```bash
conda create -n xpert python=3.10
conda activate xpert
pip install -r requirements.txt
```

关键依赖：`torch>=2.5.1`, `scanpy>=1.10.0`, `flash-attn>=2.7.0`（可选，加速训练），`unimol-tools`（UniMol 分子特征提取）。

#### 数据文件总览

| 文件 | 路径 | 大小 | 说明 |
|------|------|------|------|
| L1000 MDMT 全集 | [pretrain/data/l1000_mdmt_full_336852.h5ad](pretrain/data/l1000_mdmt_full_336852.h5ad) | ~GB | 预训练主数据集，336,852 样本 |
| L1000 MDMT 子集 | [pretrain/data/l1000_mdmt_268022_pretrain.h5ad](pretrain/data/l1000_mdmt_268022_pretrain.h5ad) | ~GB | 预训练备选子集 |
| UniMol 分子特征 | [pretrain/data/all_drugs_unimol_arr.npy](pretrain/data/all_drugs_unimol_arr.npy) | 4.2 GB | 8981 药物 × 122 原子 × 770 维 |
| 药物 SMILES 映射 | [pretrain/data/all_drugs_idx2smi_8981.npy](pretrain/data/all_drugs_idx2smi_8981.npy) | — | 药物索引 → SMILES 字典，用于在线生成 UniMol 特征 |
| KPGT 指纹 | [pretrain/data/all_drugs_idx2KPGT.npy](pretrain/data/all_drugs_idx2KPGT.npy) | — | 备选药物特征 |
| Morgan 指纹 | [pretrain/data/all_drugs_idx2morgan.npy](pretrain/data/all_drugs_idx2morgan.npy) | — | 备选药物特征 |
| PPI 基因向量 | [pretrain/data/PPI_gene_vector_128d.npy](pretrain/data/PPI_gene_vector_128d.npy) | 979 KB | 978×128，基因相互作用网络嵌入 |
| HG 药物嵌入表 | [drug_embeddings/all_drugs_idx2HG.npy](drug_embeddings/all_drugs_idx2HG.npy) | 8.8 MB | 8981×256，异构图预训练的药物嵌入 |
| CDSDB 临床数据 | [finetune/data/CDSDB_mdmt_613.h5ad](finetune/data/CDSDB_mdmt_613.h5ad) | 11 MB | 613 样本 × 978 基因 |
| PANACEA 数据 | [finetune/data/PANACEA_sdst_352.h5ad](finetune/data/PANACEA_sdst_352.h5ad) | 5.9 MB | 352 样本 × 978 基因 |
| L1000 SDST | [finetune/data/l1000_sdst_78453.h5ad](finetune/data/l1000_sdst_78453.h5ad) | 1.2 GB | 78,453 样本 |
| L1000 MDMT 子集 | [finetune/data/l1000_mdmt_68830_subset.h5ad](finetune/data/l1000_mdmt_68830_subset.h5ad) | 1.1 GB | 68,830 样本 |

#### h5ad 数据结构

每个 `.h5ad` 文件是 AnnData 对象，包含：

- **`.X`**: 扰动后表达谱矩阵（样本 × 978 基因），可能是稀疏矩阵
- **`.obsm['X_ctl']`**: 对照组表达谱矩阵（样本 × 978 基因）
- **`.obs`**: 样本元数据 DataFrame，必须包含以下列：
  - `pert_idx`: 药物索引（整数，对应 UniMol 特征的第一维）
  - `pert_id`: 药物标识符（字符串，如药物名/ID）
  - `pert_dose`: 药物剂量（连续值 μM）
  - `pert_time`: 处理时间（小时，通常 3/6/24）
  - `cell_idx`: 细胞系索引（整数）
  - `cell_iname`: 细胞系名称（CDSDB 则用 `Cancer subtype`）
  - `tissue_idx`: 组织索引（整数）
  - `split` / 自定义划分列: 标记 train/valid/test 的列

---

### 步骤 1：预训练

**目标**：在大规模 L1000 数据集上预训练 XPert 模型，学习通用的药物-细胞扰动表示。

**运行命令**：

```bash
bash pretrain/run.sh [dataset] [nfold] [config]
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `dataset` | `l1000_mdmt_full` | 数据集名，对应配置中的数据路径键。可选：`l1000_mdmt_full`（336K 样本）、`l1000_mdmt_pretrain`（268K 样本） |
| `nfold` | `pretrain_split` | 数据划分列名，从 `.obs[nfold]` 读取 train/valid/test |
| `config` | `pretrain/config/l1000` | 配置文件路径（不含 `.yaml`），相对于项目根目录 |

**脚本内部调用**（等价 CLI）：

```bash
python -m src.train \
    --mode train \
    --model XPert \
    --config pretrain/config/l1000 \
    --drug_feat unimol \
    --drug_embeds unimol,HG \
    --nfold pretrain_split \
    --dataset l1000_mdmt_full \
    --use_gradscaler True \
    --include_cell_idx True
```

**数据流**：

1. [`src/config.py:load_config()`](#src-config-py) 加载 `pretrain/config/l1000.yaml`，将所有相对路径解析为绝对路径
2. [`src/utils.py:_load_data_and_drug_feat()`](#src-utils-py) 调用 scanpy 读取 `pretrain/data/l1000_mdmt_full_336852.h5ad`，加载 `drug_embeddings/all_drugs_unimol_arr.npy`
3. [`src/utils.py:load_dataloader()`](#src-utils-py) 按 `nfold` 列划分数据，创建 [`src/dataset.py:MyDataset`](#src-dataset-py) 实例和 PyTorch DataLoader
4. [`src/model_xpert.py:XPertNet`](#src-model-xpert-py) 初始化模型：
   - 加载 PPI 基因向量 → `cell_Embeddings`
   - 加载 HG 嵌入表 `drug_embeddings/all_drugs_idx2HG.npy` → 冻结的查找表
   - 初始化双分支编码器 `AttnEncoder`（对照 SA×4，处理 CA+SA+SA+CA）
5. [`src/train.py:train()`](#src-train-py) 循环：前向 → 四损失加权 → 反向 → 早停监控
6. 输出：`pretrain/output/<dataset>/<timestamp>/<dataset>_best.pth`

**关键配置参数**（[`pretrain/config/l1000.yaml`](#pretrain-config-l1000-yaml)）：

| 参数 | 值 | 说明 |
|------|-----|------|
| `train.batch_size` | 64 | 批次大小 |
| `train.train_lr` | 0.004 | 预训练学习率（Adam） |
| `train.num_epochs` | 2500 | 最大训练轮数 |
| `train.init_epoch` | 70 | 预热轮数（期间使用 batch 级 RMS 损失） |
| `train.patience` | 50 | 早停耐心值（监控 val_loss4） |
| `train.loss_weight` | `[0.2, 0.003, 0.2, 1]` | L_pert, L_ctl, L_deg, L_pcc 权重 |
| `model.ATTN.hidden_size` | 256 | 隐藏层维度 |
| `model.ATTN.n_heads` | 8 | 注意力头数 |

---

### 步骤 2：微调

**目标**：将预训练模型在目标数据集上进行微调，适配临床或特定实验场景。

**运行命令**：

```bash
bash finetune/run.sh <pretrained_model_path> [dataset] [nfold] [config]
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `pretrained_model_path` | **必填** | 步骤 1 输出的预训练权重路径，如 `pretrain/output/pretrain_mdmt_new.pth` |
| `dataset` | `cdsdb_mdmt` | 目标数据集名。可选：`cdsdb_mdmt`（临床 613 样本）、`panacea_mdmt`（352 样本）、`l1000_sdst`（78K 样本）、`l1000_mdmt`（68K 样本） |
| `nfold` | `split_1,split_2,...,split_5` | 逗号分隔的多 fold 列名。支持单 fold（如 `split_1`）或多 fold |
| `config` | `finetune/config/cdsdb` | 配置文件路径。可选：`finetune/config/cdsdb`、`finetune/config/panacea`、`pretrain/config/l1000` |

**多 fold 示例**：

```bash
# CDSDB 5 折交叉验证微调
bash finetune/run.sh pretrain/output/pretrain_mdmt_new.pth cdsdb_mdmt split_1,split_2,split_3,split_4,split_5 finetune/config/cdsdb

# PANACEA 微调
bash finetune/run.sh pretrain/output/pretrain_mdmt_new.pth panacea_mdmt split_1,split_2,split_3,split_4,split_5 finetune/config/panacea

# L1000 SDST 微调
bash finetune/run.sh pretrain/output/pretrain_mdmt_new.pth l1000_sdst split_1 pretrain/config/l1000

# L1000 冷剂量-时间少样本微调
bash finetune/run.sh pretrain/output/pretrain_mdmt_new.pth l1000_mdmt split_cold_dose\&time_random_1_one_shot pretrain/config/l1000
```

**脚本内部调用**（等价 CLI）：

```bash
python -m src.train \
    --mode train \
    --model XPert \
    --config finetune/config/cdsdb \
    --drug_feat unimol \
    --drug_embeds unimol,HG \
    --nfold split_1,split_2,split_3,split_4,split_5 \
    --dataset cdsdb_mdmt \
    --use_gradscaler True \
    --include_cell_idx True \
    --pretrained_model pretrain/output/pretrain_mdmt_new.pth
```

**数据流**：

1. 加载配置（如 `finetune/config/cdsdb.yaml`）和目标数据集 h5ad 文件
2. 加载预训练权重，自动跳过 shape 不匹配的参数（由 `--exclude_keys` CLI 参数控制）
3. 调用 [`src/utils.py:set_requires_grad()`](#src-utils-py) 控制可训练参数范围（`finetune_params` 设为 `"all"` 则全参数训练）
4. 使用微调学习率（`finetune_lr`）初始化优化器
5. 训练循环同步骤 1，但切换为微调配色的 init_epoch 等参数
6. 输出：`finetune/output/<pretrain_dataset>-<dataset>-<nfold>/<timestamp>/`（每个 fold 一个）

**预训练权重加载逻辑**（详见 [`src/train.py:main()`](#src-train-py)）：
- 遍历预训练 checkpoint 的 `model_state_dict`
- 只加载 key 存在且 shape 匹配的参数
- 跳过 `--exclude_keys` 中指定的模块子串（如输出头的 `class_fc`、`ctl_fc`、`trt_fc`、`deg_fc`）

---

### 步骤 3：评估

**目标**：在测试集上评估微调后的模型，计算指标，提取预测谱和注意力矩阵。

**3a. CLI 推理模式**：

```bash
python -m src.train \
    --mode infer \
    --model XPert \
    --config finetune/config/cdsdb \
    --dataset cdsdb_mdmt \
    --drug_feat unimol \
    --drug_embeds unimol,HG \
    --saved_model finetune/output/cdsdb_warm_split.pth \
    --output_profile True \
    --output_attention True \
    --output_cls_embed True
```

**3b. 交互式评估**：打开 [`metrics/evaluate.ipynb`](#metrics-evaluate-ipynb) 按 cell 执行：

1. 加载预测 profile 文件（`.npy`）
2. 调用 `get_metrics_new` 计算完整指标集（含 R²、MMD、Wasserstein）
3. 调用 `get_metrics_distribution` 计算分布指标
4. 保存结果为 CSV

**输出文件**（保存在 `metrics/output/<timestamp>/`）：

| 文件 | 内容 |
|------|------|
| `<dataset>_<fold>_test_samples_result.csv` | 逐样本指标（PCC、R²、RMSE 等） |
| `<dataset>_<fold>_test_samples_prediction_profile.npy` | `{y_true, y_pred, ctl_true, deg_true, deg_pred}` |
| `<dataset>_<fold>_test_samples_attention_scores.pkl` | 逐层注意力矩阵（trt + ctl） |
| `<dataset>_<fold>_test_samples_cls_embed.npy` | CLS token 嵌入 `(N, 2, hidden_size)` |

---

### 数据流总图

```
Step 1: XPert 预训练
  输入:  pretrain/data/l1000_mdmt_full_336852.h5ad
        drug_embeddings/all_drugs_unimol_arr.npy
        drug_embeddings/PPI_gene_vector_128d.npy
        drug_embeddings/all_drugs_idx2HG.npy (预计算 HG 嵌入)
  调用:  pretrain/run.sh → python -m src.train --mode train
  输出:  pretrain/output/<dataset>/<timestamp>/*.pth
            │
            ▼
Step 2: XPert 微调
  输入:  finetune/data/{cdsdb,panacea,l1000}/*.h5ad
        drug_embeddings/* (共享特征)
        pretrain/output/*.pth (预训练权重)
  调用:  finetune/run.sh → python -m src.train --mode train --pretrained_model <path>
  输出:  finetune/output/<pretrain_dataset>-<dataset>-<nfold>/<timestamp>/*.pth
            │
            ▼
Step 3: 评估
  输入:  finetune/output/*.pth + 数据集 h5ad
  调用:  python -m src.train --mode infer 或 metrics/evaluate.ipynb
  输出:  metrics/output/<timestamp>/*.csv, *.npy, *.pkl
```

---

## 二、文件详情

### `src/train.py` {#src-train-py}

**作用**：命令行统一入口，支持 3 种运行模式（`train`/`test`/`infer`），包含完整的训练循环、验证逻辑和推理逻辑。

**`arg_parse()` 完整参数表**：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--mode` | str | `train` | 运行模式：`train`（训练+验证+测试）、`test`（仅测试）、`infer`（推理，输出预测谱和注意力） |
| `--model` | str | `XPert` | 模型名，当前仅支持 `XPert` |
| `--config` | str | `pretrain/config/l1000` | YAML 配置文件路径（不含 `.yaml`），相对于项目根目录 |
| `--dataset` | str | `l1000_sdst` | 数据集名：`l1000_sdst`、`l1000_mdmt`、`l1000_mdmt_full`、`l1000_mdmt_pretrain`、`l1000_mdmt_example`、`l1000_example`、`cdsdb_mdmt`、`panacea_mdmt` |
| `--drug_feat` | str | `unimol` | 药物特征类型：`unimol`（UniMol 3D 分子特征）、`KPGT`（KPGT 指纹）、`morgan`（Morgan 指纹）、`smi`（SMILES 字符串） |
| `--nfold` | str | `split` | 数据划分列名，支持逗号分隔多 fold（如 `split_1,split_2,split_3`） |
| `--device` | str | `cuda:0` | 训练设备 |
| `--seed` | int | `2024` | 随机种子 |
| `--pretrained_model` | str | `None` | 预训练权重路径（微调时使用） |
| `--saved_model` | str | `None` | 测试/推理时加载的模型权重路径 |
| `--saved_model_path` | str | `None` | 同上但仅指定目录，自动拼接 `{fold}_fold_early_stop.pth` |
| `--resume_from` | str | `None` | 从 checkpoint 恢复训练 |
| `--use_gradscaler` | bool | `False` | 启用混合精度训练（`torch.amp.GradScaler` + Flash Attention） |
| `--lr_scheduler` | bool | `True` | 启用学习率调度（epoch 40 后 lr×0.5） |
| `--pretrained_mode` | str | `global` | HG 嵌入模式：`global`（通用 HG 嵌入）、`specific`（药物特异性基因嵌入） |
| `--include_cell_idx` | bool | `False` | 启用 CLS token + 细胞系分类辅助任务 |
| `--wo_HG` | bool | `False` | 消融：去除 HG token |
| `--wo_atom` | bool | `False` | 消融：仅保留前 2 个 token（HG + molecule） |
| `--wo_atom_HG` | bool | `False` | 消融：仅保留 HG token |
| `--wo_unimol` | bool | `False` | 消融：仅使用 HG token |
| `--wo_ppi` | bool | `False` | 消融：去除 PPI 基因向量 |
| `--use_gene_pos_emed` | bool | `False` | 单元格外的可学习位置编码 |
| `--output_profile` | bool | `False` | 推理时输出预测谱 `.npy` 文件 |
| `--output_attention` | bool | `False` | 推理时输出注意力矩阵 `.pkl` 文件（强制 batch_size=1） |
| `--output_cls_embed` | bool | `False` | 推理时输出 CLS token 嵌入 `.npy` 文件 |
| `--weighted_loss` | bool | `False` | 备用加权损失开关 |
| `--kl_loss` | bool | `False` | 备用 KL 散度损失开关 |

**`main()` 函数流程**（`mode='train'`）：
1. 解析参数 → 加载配置 → 创建带时间戳的输出目录 → 配置 logger
2. 遍历每个 fold：
   - 加载 train/valid/test DataLoader
   - 初始化 XPertNet，kaiming 初始化权重
   - 若有 `--pretrained_model`：加载预训练权重（跳过 shape 不匹配的 key），设置微调参数范围和学习率
   - 训练循环：`train()` → `validate()` → `EarlyStopping.step(val_loss4)` → `lr_scheduler.step()`
   - 早停后加载最佳模型，在 train/valid/test 上评估并保存 CSV
   - 保存损失曲线 PNG

**损失函数（4 项加权）**：
```python
Total = 0.2 × L_pert + 0.003 × L_ctl + 0.2 × L_deg + 1.0 × L_pcc
```
- `L_pert`: MSE(预测扰动表达, 真实扰动表达)，按样本求和
- `L_ctl`: MSE(预测对照表达, 真实对照表达) 或 CrossEntropy(CLS 分类)
- `L_deg`: MSE(预测表达变化, 真实表达变化)
- `L_pcc`: Σ(1 − PCC(预测变化, 真实变化))

---

### `src/model_xpert.py` {#src-model-xpert-py}

**作用**：XPert 双分支 Transformer 的完整实现，是项目的核心模型文件。

**`get_unimol_drug_feat(input_features)`**：从 UniMol 特征数组中分离 `atom_feat`（原子嵌入）、`atom_symbol`（原子类型 ID）、`atom_musk`（注意力掩码）。

**`AttnEncoder` 类**：通用注意力编码器
- 根据 `structure` 字符串（如 `"CA+SA+SA+CA"`）动态构建编码器层序列
- `CA` → 实例化 `crossEncoder`（交叉注意力层），`SA` → 实例化 `Encoder`（自注意力层）
- 每层按序执行，可选记录注意力矩阵
- 支持在 CA 层后注入 `drug_specific_gene_embedding`（药物特异性基因嵌入）

**`XPertNet` 类**：主网络

初始化：
1. `cell_emb`: 细胞 token 嵌入（表达值 bin embedding + PPI 基因向量作为位置编码）
2. `drug_emb`: 药物 token 嵌入（UniMol 512→256 线性投影 + HG 嵌入填充首 token + MDMT 场景前插剂量/时间 token）
3. `drug_HG_embed`: 预训练 HG 嵌入表（`torch.tensor`，冻结参数）
4. 双分支编码器：`attnEncoder_ctl`（SA×4，对照分支）、`attnEncoder_trt`（CA+SA+SA+CA，处理分支）
5. 三个输出头（hidden→64→1 的 Sequential）：`ctl_fc`、`trt_fc`、`deg_fc`
6. 可选：CLS token + `class_fc` 细胞分类头

前向传播：
```
cell_tokens → ┬─→ attnEncoder_ctl (SA×4) ──→ ctl_fc ──→ x_ctl_pred
              │
              └─→ attnEncoder_trt (CA+SA+SA+CA) ──→ trt_fc ──→ x_pert_pred
                    ↑                                   │
                  drug_tokens                     deg_fc ──→ x_deg_pred
                  (UniMol+HG+Dose+Time)
```

---

### `src/model_utils.py` {#src-model-utils-py}

**作用**：XPert 网络的所有基础组件库。

| 组件 | 功能 |
|------|------|
| `LayerNorm` | 自定义 LayerNorm（`γ·(x-μ)/σ + β`） |
| `Embeddings` | 标准 Embedding：word + position embedding → LayerNorm → Dropout |
| `cell_Embeddings` | 细胞 embedding：表达值 bin embedding + PPI 基因向量（作为位置编码）。维度不匹配时自动线性投影。支持消融开关 `wo_ppi` 和额外位置编码 `use_gene_pos_emed` |
| `unimol_Embeddings` | 药物 embedding：UniMol 特征线性投影（512/2304/1024 → hidden_size）→ 首 token 填充 HG 嵌入 → MDMT 模式前插剂量/时间 learnable embedding → 加 position embedding。支持消融参数 `wo_HG`、`wo_atom`、`wo_atom_HG`、`wo_unimol` |
| `SelfAttention` | 自注意力：双路径——`output_attention=True` 用标准实现返回注意力矩阵，否则用 `F.scaled_dot_product_attention`（Flash Attention） |
| `CrossAttention` | 交叉注意力：cell 作 query，drug 作 key/value，同样支持双路径 |
| `FeedForward` | 两层 MLP：Linear→ReLU→Dropout→Linear |
| `SublayerConnection` | 残差连接：`LayerNorm(x) + input_tensor` |
| `SelfOutput` | FeedForward + 残差连接 |
| `Encoder` | 自注意力编码器：SelfAttention → 残差 → FeedForward → 残差 |
| `crossEncoder` | 交叉注意力编码器：Drug SA → Cell SA → Cell CA（cell 对 drug 交叉注意力）→ 残差 → FeedForward → 残差 |

---

### `src/dataset.py` {#src-dataset-py}

**作用**：定义 `MyDataset` 类，将 AnnData 对象转换为 PyTorch 可迭代的训练样本。

**辅助函数**：

- **`_digitize(x, bins, side)`**：将连续表达值转换为离散 bin 索引。`side="both"` 时在两个相邻 bin 之间加随机抖动实现数据增强（训练时），本项目中实际使用 `side="one"`。
- **`assign_pert_dose_corrected(pert_dose)`**：将连续剂量值映射到 10 个离散区间（0~9）。区间边界：`[0, 0.21, 0.41, 0.71, 1.01, 1.51, 3.1, 4.1, 7.1, 12.1, ∞)`。

**`MyDataset(Dataset)` 类**：

- **`__init__`**：
  1. 从 `.obs` 计算/补全 `pert_dose_idx`（调用 `assign_pert_dose_corrected`）和 `pert_time_idx`（映射：3h→0, 6h→1, 24h→2）
  2. 从 `.X` 获取扰动后表达谱，`.obsm['X_ctl']` 获取对照表达谱
  3. 根据 `min_value`/`max_value` 计算 128 个分位数 bins，`_digitize` 分箱
  4. 预提取所有元数据列为 numpy 数组（`_pert_id`, `_pert_idx`, `_pert_dose`, `_cell_id` 等）加速 `__getitem__`

- **`__getitem__`**：返回一个元组：
  ```python
  # 训练/验证/测试模式（10 个元素）：
  (trt_raw_data, ctl_raw_data, trt_raw_binned, ctl_raw_binned,
   drug_feat, pert_dose_idx, pert_time_idx, pert_idx, cell_idx, tissue_idx)

  # 推理模式（额外追加 4 个元素）：
  (..., pert_id, cell_id, pert_dose, pert_time)
  ```

---

### `src/utils.py` {#src-utils-py}

**作用**：数据加载、药物特征生成、损失函数、早停机制、训练辅助函数。

**核心组件**：

| 函数/类 | 功能 |
|---------|------|
| `atom2id` | 原子符号到 ID 的映射字典（25 种原子 + `hg_molecule` + `molecule` 特殊 token） |
| `generate_unimol_feat(smiles_dict, max_atoms)` | 调用 `unimol_tools.UniMolRepr` 为 SMILES 列表在线生成 3D 分子特征。输出格式：`{pert_idx: [mask, atom_symbol_id, 512-dim repr]}`，pad/truncate 到 `max_atoms`（默认 122） |
| `_load_data_and_drug_feat(args, config, logger)` | 根据 `dataset` 名加载对应 h5ad 文件 + 药物特征。若 `drug_unimol_path=null`，自动调用 `generate_unimol_feat` 在线生成 UniMol 特征 |
| `load_dataloader(args, config, logger, nfold)` | 按 `nfold` 列划分数据，创建 train/valid/test DataLoader（`num_workers=2`） |
| `load_test_dataloader(args, config, logger, nfold)` | 仅加载 test DataLoader |
| `load_infer_dataloader(args, config, logger, nfold)` | 推理模式 DataLoader，`output_attention=True` 时强制 batch_size=1 |
| `set_random_seed(seed)` | 设置 Python/numpy/torch/cuda 随机种子 |
| `log_nested_dict(raw_dict, logger)` | 递归日志输出嵌套字典 |
| `mse_loss_ls_sum(pred, target)` | 按基因维度求 MSE 后按样本求和（论文公式） |
| `pcc_loss_sum(pred, target)` | `sum(1 - PCC(pred_i, target_i))`，逐样本 PCC 损失求和 |
| `cos_loss_sum(pred, target)` | `sum(1 - cos_sim(pred_i, target_i))` |
| `train/validate/infer(model, ...)` | 简化版训练/验证/推理函数（2 损失），用于非 XPert 基线模型 |
| `EarlyStopping` | 早停类：监控指定指标，耐心值计数，保存/加载最佳模型 checkpoint |
| `save_loss_fig(tr_loss, val_loss, n_fold, folder)` | 绘制 train/val 损失曲线保存为 PNG |
| `set_requires_grad(model, finetune_params)` | 控制微调可训练参数范围：`"all"` 全部训练，或指定模块名列表仅训练匹配的层 |

---

### `src/config.py` {#src-config-py}

**作用**：YAML 配置加载器，将相对路径解析为绝对路径，以及输出目录生成。

| 函数 | 功能 |
|------|------|
| `load_config(config_name)` | 从 `configs/<name>.yaml` 加载配置，调用 `resolve_paths()` 处理路径，返回完整配置 dict |
| `resolve_paths(config)` | 遍历配置中的路径键：`model.ATTN.ppi_gene_vector_path` 和 dataset 路径（`l1000_sdst_data_root` 等）→ 相对 `PROJECT_ROOT`；`model.HG.*` 路径 → 相对 `WORKSPACE_ROOT`（`xpert/../`） |
| `get_output_dir(stage, timestamp)` | 根据阶段返回输出目录：`pretrain` → `pretrain/output/logs/<ts>`，`finetune` → `finetune/output/logs/<ts>`，`infer` → `metrics/output/<ts>` |

**路径解析规则**：
- `PROJECT_ROOT` = `xpert/`（`config.py` 的上上级目录）
- `WORKSPACE_ROOT` = `xpert/../`（PROJECT_ROOT 的父目录）
- xpert 项目内路径（数据集、PPI 向量、UniMol 特征）→ 相对 PROJECT_ROOT
- HG 嵌入路径（`drug_hg_pretrained_embed_path`、`hg_path`）→ 相对 WORKSPACE_ROOT

---

### `src/eval_metrics.py` {#src-eval-metrics-py}

**作用**：包含 R²、MMD、Wasserstein 等论文所需的分布级别指标，以及完整的 `get_metrics_new` 函数。

| 函数 | 功能 |
|------|------|
| `precision_k(y_true, y_pred, k)` | Top-K 精确度：预测的前 K 个上/下调基因与真实值的重叠比例 |
| `mae/mse/rmse(y_true, y_pred)` | 标准误差指标，`flag=True` 返回逐样本列表 |
| `pearson/spearman(y_true, y_pred)` | 逐样本计算相关系数后取均值 |
| `compute_mmd(pred, true)` | 最大均值差异（RBF 核），逐基因计算后取均值 |
| `compute_wasserstein(pred, true)` | Wasserstein 距离 |
| `get_metrics(y, f)` | 基础指标集：MSE, RMSE, MAE, Pearson, Spearman, Precision@K（4 种 k 值） |
| `get_metrics_new(y, f, ctl)` | **主评估函数**：所有 x_pert 指标（含 R²、MMD、Wasserstein）+ 所有 x_deg 指标（R²_deg, Pearson_deg, MMD_deg, Wasserstein_deg, Precision@K_deg） |
| `get_metrics_distribution(y, f, ctl)` | 仅分布指标：MMD, Wasserstein（用于 x_pert 和 x_deg） |

**关键细节**：
- x_deg = y - ctl（真实表达变化 = 扰动后 − 对照组）
- MMD/Wasserstein 是逐基因计算后取均值（非全局分布）
- `get_metrics_new` 是论文报告指标的来源

---

### `src/metrics.py` {#src-metrics-py}

**作用**：不含分布指标的评估函数集。函数命名和功能与 `eval_metrics.py` 的前半部分相同（独立可用，无依赖关系），额外包含逐样本版本的 `get_metrics_infer`。

| 函数 | 功能 |
|------|------|
| `precision_k`, `mae`, `rmse`, `mse`, `pearson`, `spearman` | 同 `eval_metrics.py` 的对应函数 |
| `get_metrics(y, f)` | 基础指标（不含 deg，不含 R²/MMD/Wasserstein） |
| `get_metrics_new(y, f, ctl)` | x_pert + x_deg 的误差+相关性+Precision@K 指标（**不含** R²/MMD/Wasserstein），比 `eval_metrics.py` 版本更轻量 |
| `get_metrics_infer(y, f, ctl)` | **推理专用**：逐样本返回全部指标列表（`flag=True` 全部启用），返回 `(mean_metrics_dict, per_sample_metrics_dict)` |

**与 `eval_metrics.py` 的分工**：
- [metrics.py](#src-metrics-py)：误差 + 相关性 + Precision@K，用于训练/验证快速评估和推理逐样本分析
- [eval_metrics.py](#src-eval-metrics-py)：上述全部 + R² + MMD + Wasserstein，用于论文级别完整报告

---

### `src/vis_gene_attn.py` {#src-vis-gene-attn-py}

**作用**：基因级注意力矩阵可视化。

| 函数 | 功能 |
|------|------|
| `plot_attention_heatmap(attention_matrix, ...)` | 基因注意力热力图：样本×基因矩阵，seaborn heatmap。>30 基因时自动采样 x 轴标签 |
| `plot_attention_boxplot(attention_matrix, ...)` | 基因注意力箱线图：每基因一个箱线，展示跨样本注意力分布，Spectral 调色板 |
| `plot_attention_scatter(attention_vectors, ...)` | 多 fold 注意力散点图：叠加不同 fold 的注意力向量，计算并标注余弦相似度和变异系数（CV） |

---

### `src/vis_molecule.py` {#src-vis-molecule-py}

**作用**：药物分子级注意力可视化，将注意力权重映射到 2D 分子结构图上。

| 函数 | 功能 |
|------|------|
| `visualize_molecule_attention(smiles, attn_weights, ...)` | **核心函数**：SMILES → RDKit Mol → 2D 坐标 → Cairo 渲染 → PIL 在每原子坐标处叠加注意力渐变圆。支持 colorbar 和归一化 |
| `plot_attention_heatmap(attention_matrix, ...)` | 原子注意力热力图：样本×原子矩阵 |
| `plot_attention_scatter(attention_vectors, ...)` | 多 fold 原子注意力散点图 |

---

### `src/__init__.py` {#src-init-py}

**作用**：使用 `__getattr__` 实现 lazy import，支持 `from src import XPertNet` 等简洁用法。

导出的符号：`XPertNet`, `MyDataset`, `get_metrics_new`, `get_metrics_infer`, `get_metrics_distribution`, `EarlyStopping`, `set_random_seed`。

---

### `pretrain/config/l1000.yaml` {#pretrain-config-l1000-yaml}

**用途**：L1000 数据集主配置，适用于预训练和 L1000 相关的微调。

**关键参数**：

| Section | 参数 | 值 | 说明 |
|---------|------|-----|------|
| ATTN | `hidden_size` | 256 | 隐藏层维度 |
| ATTN | `n_heads` | 8 | 注意力头数 |
| ATTN | `trt_structure` | `CA+SA+SA+CA` | 处理组编码器结构 |
| ATTN | `ctl_structure` | `SA+SA+SA+SA` | 对照组编码器结构 |
| dataset | `max_value` / `min_value` | 15.002 / 0.0 | L1000 表达值范围 |
| dataset | `n_bins` | 128 | 表达值分箱数 |
| dataset | `gene_num` | 978 | 基因数 |
| dataset | `max_atom_size` | 122 | 药物最大原子数 |
| dataset | `num_pert_dose` | 10 | L1000 剂量区间数 |
| dataset | `num_pert_time` | 3 | L1000 时间区间数（3h/6h/24h） |
| dataset | `num_cell_id` | 240 | L1000 细胞系数 |
| train | `train_lr` | 0.004 | 预训练学习率 |
| train | `batch_size` | 64 | 批次大小 |
| train | `num_epochs` | 2500 | 最大训练轮数 |
| train | `init_epoch` | 70 | 预热轮数 |
| train | `patience` | 50 | 早停耐心值 |
| train | `loss_weight` | `[0.2, 0.003, 0.2, 1]` | 四损失权重 |
| finetune | `finetune_lr` | 0.001 | 微调学习率 |
| finetune | `finetune_params` | `all` | 全参数微调 |

---

### `finetune/config/cdsdb.yaml` {#finetune-config-cdsdb-yaml}

**用途**：临床 CDSDB 数据集微调配置。与 `pretrain/config/l1000.yaml` 的关键差异：

| 参数 | CDSDB | L1000 | 原因 |
|------|-------|-------|------|
| `max_value` / `min_value` | 15.993 / 1.854 | 15.002 / 0.0 | 临床数据范围不同 |
| `num_pert_dose` | 11 | 10 | CDSDB 剂量区间更多 |
| `num_pert_time` | 8 | 3 | CDSDB 时间点多 |
| `num_cell_id` | 14 | 240 | 癌种类别数 |
| `batch_size` | 32 | 64 | 小数据集减 batch |
| `init_epoch` | 1000 | 70 | 极高预热防止过拟合 |
| `finetune_lr` | 0.002 | 0.001 | 稍大的微调学习率 |

---

### `finetune/config/panacea.yaml` {#finetune-config-panacea-yaml}

**用途**：PANACEA 独立验证数据集配置。

与 L1000 的关键差异：`batch_size: 32`、`max_value: 14.305`、`init_epoch: 500`、`finetune_lr: 0.002`。

---

### 冷剂量-时间（CDT）和 Example 场景

CDT 少样本微调和快速冒烟测试均复用 `pretrain/config/l1000.yaml`，通过 CLI 环境变量覆盖超参数（如 `NUM_EPOCHS=4`、`BATCH_SIZE=2`），不再需要独立的 YAML 配置。

---

### `pretrain/run.sh` {#pretrain-run-sh}

**用法**：`bash pretrain/run.sh [dataset] [nfold] [config]`

**固定参数**：`--mode train --model XPert --drug_feat unimol --use_gradscaler True --include_cell_idx True`

### `finetune/run.sh` {#finetune-run-sh}

**用法**：`bash finetune/run.sh <pretrained_model_path> [dataset] [nfold] [config]`

**固定参数**：同 pretrain，额外追加 `--pretrained_model` 传入预训练权重路径。

### `metrics/evaluate.ipynb` {#metrics-evaluate-ipynb}

**作用**：交互式评估 Notebook。加载预测输出 `.npy` 文件，调用 `get_metrics_new` 和 `get_metrics_distribution` 计算完整指标，展示并保存结果。

Cell 结构：
1. 导入 `src.eval_metrics` 和 `src.metrics`
2. 加载预测 profile（`*_predict_profile.npy`），提取 y_true/y_pred/ctl_true/deg_true/deg_pred
3. 调用 `get_metrics_new` → DataFrame 展示
4. 调用 `get_metrics_distribution` → DataFrame 展示
5. 保存结果为 `output/csv/evaluation_results.csv`
