# miri_reposition/xpert-fda-v1 — MIRI 药物重定位工作流

对 3,127 个 FDA 批准小分子药物进行 MIRI（心肌缺血再灌注损伤）疗效筛选与排序的工作流。这是旧版 `fda_miri_vs/` 管线的重构版（旧版已冻结于 `miri_reposition/xpert-fda-v1_history/`）。

## 管线总览

```
miri_reposition/xpert-fda-v1/
├── screen.ipynb              # 唯一编排器：按顺序运行 stages 01-10
├── scripts/                  # 编号阶段脚本 (01-10)
├── finetune/                 # 工作流拥有的 DTOXS 微调 (configs, 权重)
├── predict/                  # 扰动转录组预测实现
├── data/                     # 各阶段产物
├── reference-genesets/       # MIRI 参考基因集 (07 阶段可再生成)
├── prompts/                  # gpt-researcher 提示词
└── reports/, plots/          # 生成报告
```

## 十个阶段

| 阶段 | 输入 | 输出 | 说明 |
|:----:|------|------|------|
| 01 | predict/fda-cm/fda_drugs.csv | data/drug_list.csv | 3,127 药物列表 + 嵌入索引 |
| 02 | 冻结嵌入 + HG | data/embeddings/ | 幂等：完整 HG 表存在则跳过 |
| 03 | DTOXS 对照孔 | data/baseline/ | 工作流拥有的共识基线 |
| 04 | dtoxs.h5ad + 预训练权重 | finetune/output/ | DTOXS 微调 |
| 05 | 基线 + 药物列表 + 微调模型 | data/deg/ | XPert 预测 (索引预计算 UniMol) |
| 06 | data/deg/ | data/deg-zs/ | per-gene z-score |
| 07 | disease_data/miri/{GSE}/deg/ | reference-genesets/ | 集成 MIRI 参考 |
| 08 | deg-zs + 参考基因集 | data/scores/1-scores.csv | 3 评分器 → 4 Z 分量 |
| 09 | scores + ranked | data/reports/ | gpt-researcher 报告 (agent 驱动) |
| 10 | 1-scores.csv | data/ranked/ranked_drugs.csv | WCS<0 + 关键通路过滤 → 排名 |

## 评分方法

- **WCS + RGES**：基于 MIRI 上下调基因集的加权连通性/逆转基因表达评分
- **Zhang / ssCMap**：基于签名基因排名的连通性评分
- **Pathway RGES**：40 条 MIRI 通路（20 UP + 20 DOWN，8 条关键）的逐通路 RGES

整合公式：`integrated = −0.3·RGES_zs − 0.2·Zhang_zs + 0.3·pathway_reversal_zs + 0.025·Σcritical_efficacy_zs`，随后过滤（WCS<0 + ≥6/8 关键通路）并排序。

心脏毒性评分从 v1 推迟（见 `docs/adr/0004`），不设占位阶段。

## 快速开始

```bash
conda activate arc
cd miri_reposition/xpert-fda-v1/
jupyter notebook screen.ipynb
```

## 关键设计

- **组件即库 (ADR 0002)**：`pert_model/xpert/`、`drug_emb_model/hg/` 只提供库代码 + 权重；工作流脚本生成所有产物。
- **嵌入一致性**：finetune 与 predict 均按药物索引消费冻结的预计算 UniMol 数组，避免特征漂移。
- **参考可再生成**：07 阶段从 `disease_data/miri/{GSE}/deg/` 重建集成参考。
- **幂等嵌入**：02 阶段在完整 HG 表存在时跳过昂贵的推理。
