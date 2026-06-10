# XPert-Rank: 引入排序损失的药物扰动预测模型

> 状态：设计方案，待实施。所有修改需要重新预训练，暂不部署。

## 1. 概述

XPert-Rank 是对现有 XPert 模型的改进版本，**架构不变**，仅修改训练损失函数：在原有的 PCC_deg loss 之外新增 Spearman_deg（排序相关性）损失项，使模型在训练时直接优化基因 DEG 的排序一致性。

### 1.1 改动动机

下游 MIRI 药效评分（`xpert_fda-miri_scores/`）主要依赖基因 DEG 的排序关系而非绝对数值：

| 评分方法 | 排序依赖 | 幅值依赖 |
|---------|---------|---------|
| Zhang Score / ssCMap | **纯排名驱动**：仅使用 `signed_rank_by_abs_value(DEG)`，不直接使用 DEG 原值 | 仅在按 \|DEG\| 赋 signed rank 时有间接影响 |
| WCS + RGES | **排名为主**：ES 由 hit gene 在 DEG 降序列表中的位置决定 | \|DEG\|¹ 作为 GSEA running sum 权重，起辅助作用 |
| Pathway RGES | 同 WCS/RGES | 同 WCS/RGES |

因此 Spearman（衡量排序一致性）比 Pearson（衡量线性一致性，对幅值偏差敏感）更直接地对齐下游任务目标。当前模型用 PCC_deg (weight=1.0) 作为最重的损失项，将大量 capacity 花在拟合 DEG 精确幅值上，而这些幅值在下游 z-score 归一化和 rank 操作中被消除。

### 1.2 与现有 XPert 的关系

- **模型架构**：100% 复用 `xpert/src/model_xpert.py`（XPertNet + AttnEncoder），不做任何修改
- **训练流程**：复用 `xpert/src/train.py` 的 train/validate/infer 框架，仅在损失计算中添加 Spearman 项
- **预训练权重**：**不兼容**。损失函数变了，需要从 L1000 预训练开始重新训练
- **下游使用**：FDA 预测和药效评分的输入格式不变（仍是 978 基因 DEG CSV）

## 2. 模型架构

与 XPert 完全相同（~25M 参数），在此简要回顾：

```
输入: Drug(SMILES) + Cell_State(978-gene binned expression)
  │
  ├─ Drug Branch:
  │   UniMol(122×512) → Linear(512→256)
  │   + HG_embed(256) 替换首 token
  │   + dose/time learnable embeddings
  │
  ├─ Cell Branch:
  │   Binned expression embed(128 vocab × 256)
  │   + PPI gene positional embed(128→256)
  │
  ├─ Pert Encoder (trt_structure): CA + SA + SA + CA
  │   Cross-attention 将 drug 信息融入 cell 表征
  │
  ├─ Ctl Encoder (ctl_structure): SA + SA + SA + SA
  │   Self-attention over cell state only
  │
  └─ Output Heads (三个独立 MLP, 256→64→1):
      ├─ x_pert (预测扰动后表达)
      ├─ x_ctl  (预测对照表达)
      └─ x_deg = x_pert - x_ctl (表达变化)
```

## 3. 损失函数设计（核心改动）

### 3.1 当前 XPert 损失

```
L = 0.2·MSE(x_pert, y_pert)        # L1: 扰动表达 MSE
  + 0.003·MSE(x_ctl, y_ctl)        # L2: 对照表达 MSE
  + 0.2·MSE(x_deg, y_deg)          # L3: DEG MSE
  + 1.0·(1 - PCC(x_deg, y_deg))    # L4: DEG Pearson 损失（最重项）
```

### 3.2 XPert-Rank 损失

```
L = 0.2·MSE(x_pert, y_pert)        # L1: 不变
  + 0.003·MSE(x_ctl, y_ctl)        # L2: 不变
  + 0.2·MSE(x_deg, y_deg)          # L3: 不变
  + 0.7·(1 - PCC(x_deg, y_deg))    # L4: PCC 权重从 1.0 降至 0.7
  + 0.5·(1 - Spearman(x_deg, y_deg)) # L5: 新增 Spearman 损失
```

### 3.3 Spearman Loss 的可导实现

Spearman = Pearson(rank(pred), rank(target))。排序操作不可导，采用 **sigmoid 逐对比较软排序**（参考 Blondel et al., ICML 2020）：

```
软排序: rank_i ≈ 1 + Σ_{j≠i} sigmoid((pred_i - pred_j) / τ)

其中 τ 为温度参数（默认 1.0），越小越接近硬排序但梯度越稀疏。
```

具体实现：

```python
def spearman_loss_sum(pred, target, temperature=1.0):
    """
    可导 Spearman 排序相关性损失（逐样本求和）。

    Args:
        pred:   (B, N) 预测 DEG 值, B=batch_size, N=978 genes
        target: (B, N) 真实 DEG 值
        temperature: sigmoid 温度参数

    Returns:
        loss: 标量 = Σ(1 - soft_spearman_i)
    """
    batch_size, n_genes = pred.shape

    # 1. 目标值硬排序（detached，不需要梯度）
    target_ranks = torch.argsort(
        torch.argsort(target, dim=1), dim=1
    ).float()

    # 2. 预测值软排序（可导）
    pred_diff = pred.unsqueeze(1) - pred.unsqueeze(2)   # (B, N, N)
    pred_soft_compare = torch.sigmoid(pred_diff / temperature)
    pred_soft_ranks = 1.0 + pred_soft_compare.sum(dim=2) # (B, N)

    # 3. 在排序空间计算 Pearson = Spearman
    pred_rank_mean = pred_soft_ranks.mean(dim=1, keepdim=True)
    target_rank_mean = target_ranks.mean(dim=1, keepdim=True)

    pred_centered = pred_soft_ranks - pred_rank_mean
    target_centered = target_ranks - target_rank_mean

    numerator = (pred_centered * target_centered).sum(dim=1)
    denominator = torch.sqrt(
        (pred_centered ** 2).sum(dim=1) *
        (target_centered ** 2).sum(dim=1) + 1e-8
    )

    spearman = numerator / denominator
    return (1.0 - spearman).sum()
```

### 3.4 计算开销

- 逐对比较矩阵 (B, N, N) = (64, 978, 978)
- ~61M 次 sigmoid 计算/step
- GPU 开销远小于模型 forward pass，训练时间增加预计 <5%

### 3.5 权重设计依据

| 损失项 | 旧权重 | 新权重 | 说明 |
|-------|--------|--------|------|
| MSE(trt) | 0.2 | 0.2 | 不变 |
| MSE(ctl) | 0.003 | 0.003 | 不变 |
| MSE(deg) | 0.2 | 0.2 | 不变 |
| 1-PCC(deg) | 1.0 | 0.7 | 略降，为 Spearman 腾出权重空间。PCC 仍重要：WCS/RGES 的 \|DEG\| 加权和 Zhang 的 signed rank 赋权都依赖 DEG 相对幅值 |
| 1-Spearman(deg) | — | 0.5 | **新增**。直接优化排序，与下游评分目标对齐 |

总相关性 loss 权重：0.7 + 0.5 = 1.2（略高于原来的 1.0）。

## 4. 训练配置

### 4.1 预训练（L1000）

```
数据集:     l1000_mdmt_full
Split:      pretrain_split
Epochs:     2500
Batch size: 64
Optimizer:  Adam
LR:         0.004
LR schedule: 1.0 (epoch 0-39), 0.5 (epoch 40+)
Weight decay: 1e-5
Loss weight: 0.2,0.003,0.2,0.7,0.5    # 新增第5个值
Warmup:     70 epochs (batch-normalized loss → full weighted loss)
Early stop: patience=50, metric=PCC (监控 val loss4)
GradScaler: 可选
```

### 4.2 微调（DTOXS + CDSDB + PANACEA）

```
LR:         0.001
Weight decay: 1e-5
Finetune params: all
其他同预训练
```

### 4.3 CLI 参数变更

`--loss_weight` 参数从 4 个值扩展为 5 个值：

```bash
# 旧（4 个值）
--loss_weight 0.2,0.003,0.2,1

# 新（5 个值）
--loss_weight 0.2,0.003,0.2,0.7,0.5
#                              ^^^  ^^^
#                              PCC  Spearman
```

## 5. 实施清单

所有修改集中在两个文件：

| 文件 | 修改内容 |
|------|---------|
| `xpert/src/utils.py` | 在 `pcc_loss_sum()` （L282）之后添加 `spearman_loss_sum()` |
| `xpert/src/train.py` | (1) import `spearman_loss_sum`; (2) `--loss_weight` 默认值增加第 5 个值; (3) `train()`/`validate()` 解包 5 个权重，添加 `loss5` 计算和 `total_loss5` 追踪; (4) 更新日志输出包含 `Loss5` |

模型架构文件 `xpert/src/model_xpert.py` 无需修改。

## 6. 预期效果

1. **Spearman_deg 指标提升**：直接优化目标，应有显著提升
2. **PCC_deg 指标可能略降**：权重从 1.0 降到 0.7，且 Spearman 和 PCC 优化方向不完全一致
3. **下游 Zhang Score 区分度提升**：Zhang 是纯排名方法，与 Spearman 直接对齐
4. **跨方法排名一致性可能改善**：当基因排序更准确时，不同评分方法的排名应该更一致

## 7. 延后原因

- **预训练开销大**：L1000 预训练（2500 epochs × ~9000 drugs）是最大开销，损失函数任何改动都需重新预训练
- **连锁效应**：预训练 → 3 个数据集的微调 → FDA 预测 → 药效评分，全流程需重跑
- **策略选择**：将所有需要重新预训练的改进方案汇总到 `new-model/` 中，统一评估后一次性实施

## 8. 相关文件

- 现有模型损失函数：[xpert/src/train.py](../xpert/src/train.py) (L113-114), [xpert/src/utils.py](../xpert/src/utils.py) (L266-282)
- 现有模型架构：[xpert/src/model_xpert.py](../xpert/src/model_xpert.py)
- 下游评分方法：[xpert_fda-miri_scores/scripts/](../xpert_fda-miri_scores/scripts/)
- 详细分析（原始记录）：[xpert/docs/spearman_loss_analysis.md](../xpert/docs/spearman_loss_analysis.md)（待删除）

## 参考文献

- Blondel, M., et al. (2020). Fast Differentiable Sorting and Ranking. *ICML 2020*.
- Zhang, S. & Gant, T. (2008). A simple and robust method for connecting small-molecule drugs using gene-expression signatures. *BMC Bioinformatics*, 9, 258.
- Subramanian, A., et al. (2017). A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*, 171(6), 1437–1452.
- Chen, B., et al. (2017). Reversal of cancer gene expression correlates with drug efficacy and reveals therapeutic targets. *Nature Communications*, 8, 16022.