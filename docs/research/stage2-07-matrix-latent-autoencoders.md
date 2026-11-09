# Stage 2 Step 07 矩阵潜变量的转录组自编码器

研究日期 2026-11-10。这里把「矩阵输出」理解为每条表达谱产生可直接供解码器和下游处理分支使用的固定形状 `K × D` 潜表示。每行由独立的查询或槽位产生，后续计算沿 token 轴进行。逐基因上下文 `G × D`、长度随可见基因变化的 `L × D`、单个池化向量，以及把向量重新排列成矩阵，分别对应不同的接口。矩阵能缩短 token 序列，标量维数是否减少要看 `K × D` 的乘积。[ADR 0051](../adr/0051-latent-xpert.md) · [ADR 0052](../adr/0052-stage2-autoencoder-reconstruction-boundary.md)

## 当前 AAE 的结果

`AAE_512` 已输出每谱 **512 × 128** 的潜矩阵，共 65,536 个标量，输入为 18,164 个基因的 `log1p(CP10K)` 表达。编码器用 512 个查询读取带固定正弦位置的表达 token，解码器用相同的基因位置查询潜 token。两端各有一次交叉注意力，潜 token 之间没有自注意力或前馈层。冻结配置训练四轮，以 tenv 验证 MSE 选检查点。保留测试板上，相对训练集逐基因均值的平方误差，AAE 在 Tahoe/tenv 分别为 **1.133/0.496**，PCA_512 为 **0.409/0.121**。这个实验只比较已训练配置，无法归因到位置编码、模型深度或训练轮数中的某一项。[Step 07 说明](../../stage2/07_autoencoder/README.md) · [测试分数](../../stage2/07_autoencoder/output/test_scores.json)

逐轮验证记录显示，四轮内 tenv MSE 为 0.186048 → 0.183549 → 0.151560 → 0.129115，Tahoe 为 0.029343 → 0.016502 → 0.011645 → 0.009630。第四轮是最后一轮，也是选中的检查点，拟合耗时 5 分 32 秒。两源误差到最后一轮仍在下降，延长原架构的训练可作为对照，用来区分训练预算与架构差异。[逐轮记录](../../stage2/07_autoencoder/output/models/AAE_512/epochs.json) · [拟合记录](../../stage2/07_autoencoder/execution/aae_record/run.json)

## CellNavi 的论文和公开代码

**论文 Methods 的细胞流形模型。** 预训练对每个细胞的原始计数做二项降深，降深因子从 `[1, 20]` 取样。模型将可学习的基因身份、表达量和基因图信息送入六层 GeneGraph 编码器。训练时最多抽取 2,048 个非零表达基因，推理可使用全部非零基因。编码器产生含 `[CLS]` 的 `(L+1) × 256` 状态，论文定义的细胞坐标则取 `[CLS]` 的 256 维状态并投影到 **2,048 维向量**。Methods 描述的解码器把该坐标与各目标基因嵌入拼接，再用共享的两层 MLP 回归逐基因表达，对完整表达计算 MSE。论文的下游任务使用细胞坐标判断什么基因驱动源细胞到目标细胞的转变，预训练解码器随后丢弃。[CellNavi 原文 Fig. 1 与 Methods](https://www.nature.com/articles/s41556-025-01755-1)

**公开代码的预训练路径。** `PretrainModel.forward` 将完整的 `[B, L+1, 256]` 编码序列与目标基因的可学习嵌入交给两层 `TransformerDecoder`，以输入基因和目标基因两组回归 MSE 与分箱交叉熵训练。这个解码器消费的是逐基因序列，长度随输入基因数变化。公开微调代码则直接用首 token 的 256 维状态接分类头。复现时需要明确采用论文 Methods 的单向量路径，还是公开实现的逐基因解码路径。对 Step 07 可分别借鉴可学习基因身份和降噪训练。降深应当在原始计数上实施，并与 ADR 0052 的完整输入重建分开评估。[预训练代码](https://github.com/DLS5-Omics/CellNavi/blob/main/cellnavi/model/pretrain_model.py) · [微调代码](https://github.com/DLS5-Omics/CellNavi/blob/main/cellnavi/model/finetune_model.py)

## 已报道的矩阵表示

| 模型 | 每谱表示与训练目标 | 与 Step 07 的关系 |
| --- | --- | --- |
| **scLDM，Palla 等，ICML 2026** | 用基因 ID 和计数构造 token，固定数量的可学习查询读取基因集合，潜 token 再经 Transformer 处理。VAE 超参数表搜索 `K ∈ {64, 128, 256}`、每 token 的潜变量通道 `D ∈ {8, 16, 32}`，编码和解码潜块各搜 2 或 4 层。基因查询从潜矩阵还原表达，输出可用 Gaussian 或负二项头。§3.2 也给出只输出均值的确定性 AE 路径。 | 最贴近完整表达的矩阵重建。编码可筛非零基因，解码仍能查询包括实测零在内的完整基因轴。VAE 单独训练，扩散先验属于第二阶段。作者公开实现、配置和权重。观察性单细胞实验的基因轴为 17,002、19,734 或 27,997，Replogle 和 Parse1M 的扰动实验则用 2,000 个 HVG。附录 K.1 显示更深的模型也可能降低 Dentate Gyrus 的表现，仍需在 ARC 的表达尺度上验证。[正式论文](https://proceedings.mlr.press/v306/palla26a.html) · [§3.2 与附录 H](https://arxiv.org/html/2511.02986v2) · [实现](https://github.com/czi-ai/scldm) |
| **ExpressionVAE，Bhattacharya 等，2026 预印本及 ICML GenBio 工作坊** | 编码与解码采用双向的交叉注意力压缩和展开。连续 Gaussian 版本的潜矩阵为 **64 × 8**，即 512 个标量。论文主推 FSQ 版本，输出 64 个离散码索引，词表有 512 项。 | 连续版给出了与 PCA_512 对齐的矩阵预算。论文主实验用 2,000 个基因评估单细胞扰动，结果需在 Step 07 的全轴、伪汇总表达和保留板上重新检验。[预印本](https://www.biorxiv.org/content/10.64898/2026.06.15.732063v1.full) · [附录 A](https://openreview.net/pdf?id=kiYPCGFPO3) · [实现](https://github.com/sanjukta7/sc-evae) |
| **Exceiver，Connell 等，2022 预印本** | 确定性 Perceiver IO，基因嵌入乘表达值，交叉注意力读入固定 latent queries，潜 token 自注意力，再用未乘表达的基因嵌入作解码查询。官方命令示例为 128 个 latent queries。覆盖 CCLE 定义的 19,067 个蛋白编码基因，训练使用离散噪声遮盖和遮盖位置 MSE，可加组织标签辅助分类。 | 提供近全基因轴、无 KL 的架构先例。可借鉴基因身份和潜 token 自注意力。公开模型验证的是遮盖恢复，改为完整输入并对全轴评分属于 ARC 的实验设计。[原文 §6](https://arxiv.org/pdf/2210.14330) · [代码和预训练权重](https://github.com/keiserlab/exceiver) |
| **PETRI，Conrad 等，ICLR 2026** | 转录组 Perceiver 编码器输出 **8 × 768** 的每细胞表示。训练遮盖 75% 输入 token，按共同实验条件组合细胞，在多模态集合中联合重建。HepG2 转录组解码覆盖 18,082 基因。 | 支持少量 token 的近全轴遮盖重建。可参考其单模态转录组 MAE 对照。PETRI 主任务使用跨细胞和跨模态上下文。[正式论文](https://proceedings.iclr.cc/paper_files/paper/2026/hash/f52329d64e0a0110408c6717e333c0b1-Abstract-Conference.html) · [架构与附录](https://openreview.net/pdf?id=Vu8YXDooG5) |
| **Funomics T0，ICLR 2026 工作坊** | 在单细胞和 bulk RNA-seq 上联合训练。固定 Perceiver latent 可读取表达和实验条件，再经基因查询解码，任务包括遮盖计数重建。附录给出的潜隐藏状态为 **512 × 1,024**。 | 是 bulk 与单细胞共用矩阵表示的先例。它缩短约 19,000 基因的 token 序列，潜矩阵仍有 524,288 个标量。[原文 §2 与附录 F](https://openreview.net/pdf?id=exRCSDdksX) |
| **scFoundation，Hao 等，Nature Methods 2024** | 非对称遮盖自编码器。官方 `gene` / `gene_batch` 接口可导出每谱 **19,264 × 512** 的上下文基因矩阵，`cell` 接口才做池化。 | 如果矩阵的行必须对应具体基因，这是可直接核对的先例。逐基因输出的成本与基因数成正比。[论文](https://doi.org/10.1038/s41592-024-02305-7) · [实现与输出形状](https://github.com/biomap-research/scFoundation#gears) |

CZI/Biohub 的 scLDM 与 Rogers 等同名方法是两项工作。Rogers 等在 §2.7 报告的是每谱 100 维的向量 VAE。[Rogers 等原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13661255/)

## 与现有基准的对齐

可先按 scLDM 的确定性路径与 Exceiver 的结构设计一个 `64 × 8` AE。保留可学习的基因 ID、潜 token 之间的自注意力，以及由基因身份发起的全轴解码。内部注意力宽度可以大于 8，但传给解码器和 Step 08 的瓶颈保持 512 个标量，使其与 PCA_512、MLP_512 对齐。沿用 ADR 0052 的完整输入、联合训练板、两源分报和全轴 MSE。若需定位单项改动的影响，可分别对照固定正弦位置、可学习基因 ID 和潜自注意力，并记录标量预算与运行成本。外部论文支持架构可行性，ARC 的重建收益仍待在保留板验证。[ADR 0052](../adr/0052-stage2-autoencoder-reconstruction-boundary.md) · [scLDM §3.2 与附录 H](https://arxiv.org/html/2511.02986v2) · [Exceiver §6](https://arxiv.org/pdf/2210.14330)
