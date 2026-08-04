# reference/ — 参考文献存档

## 目录用途

本目录存放与项目相关的学术论文 PDF 文件，作为领域背景阅读材料。这些论文与主工作流（药物转录组扰动预测 + MIRI 药物排名）**没有代码级依赖关系**，仅供研究人员了解相关领域进展和项目方法论基础。

## 论文列表

| 文件名 | 作者/年份 | 主题 | 与项目的关系 |
|--------|-----------|------|-------------|
| `2026_PerturbHD.pdf` | 待确认 / 2026 | "Are Current AI Virtual Cell Models Useful for Scientific Discovery?" — 评估 AI 虚拟细胞模型在科学发现中的实际效用 | 与 XPert、HG 等模型的评估方向和实际应用价值判断相关 |
| `2026_t world.pdf` | 待确认 / 2026 | 推测为 "Transcriptome World"，转录组建模或扰动空间表示学习 | 可能涉及转录组扰动空间的表示方法，与 XPert 的 gene perturbation prediction 任务相关 |
| `Deep representation learning of chemical-induced transcriptional profile for phenotype-based drug discovery.pdf` | Tong et al. / 2024 (Nature Communications) | 基于深度学习的化学诱导转录谱表示，用于表型驱动的药物发现 | 与 XPert 项目核心任务直接对应：从化学结构预测转录组扰动，再基于扰动进行药物筛选 |

## 与主工作流的关系

本目录中的论文不参与任何代码执行或数据处理流程。主工作流的完整文档见根目录 [CLAUDE.md](../CLAUDE.md)，各子模块的详细说明见对应的 `CLAUDE.md` 文件：

- `xpert/` — 药物转录组扰动预测模型（XPert）
- `kpgt/` — 分子图 Transformer（KPGT）
- `hg/` — 异质图神经网络（HeteroGNN）
- `drug_ranking/` — 药物排名与评分
