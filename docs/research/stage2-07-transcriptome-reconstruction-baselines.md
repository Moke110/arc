# Stage 2 Step 07 全转录组重建方法与基线

研究日期 2026-11-10。范围为同一表达谱经编码再解码的重建、遮盖或降深后的恢复，以及它们与下游扰动预测的边界。这里的全转录组指预先固定且明确报告的测量基因轴。对 Step 07，这个轴是 18,164 个 CM 基因。本文整理方法和证据，供后续实验设计使用。[Step 07 设计记录](../adr/0051-latent-xpert.md)

## 最贴近 Step 07 的证据

Fu 等的 **ReconEval**（2026 年 bioRxiv 预印本）把三个实验分开：完整表达输入的端到端重建、冻结基础模型编码器后训练独立解码器、潜空间经扰动模型变换后再解码。端到端比较 PCA、确定性 AE 和 scVI 类 VAE。冻结编码器比较 State SE、scConcept、scGPT、SCimilarity，分别搭配 MLP、kNN 或 Transformer 解码器，并与相同标量维数的 PCA 对比。指标包括单细胞 MSE、群体均值 R²、能量距离，以及差异表达、共表达和通路信号等。该研究报告低维时确定性 AE 的重建表现较好，PCA 随维数升高变得有竞争力，冻结表示的重建表现依赖解码器。这些结论来自其数据和任务设置。[ReconEval 原文](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full) · [官方实现](https://github.com/theislab/ReconEval)

ReconEval 先归一化到 10,000 再做 log1p，但模型在各数据集自己的筛选基因集训练和测试。Tahoe-100M 是按板选 HVG 加细胞周期标记，共 **6,087** 个基因。PBMC-10M 为 **5,412** 个，LuCA 为 **5,983** 个。因此可借用它的评估设计，Step 07 仍需在自己的 18,164 基因轴验证重建。ReconEval 的 Tahoe 使用单细胞数据，Step 07 使用伪汇总表达谱。[ReconEval Methods](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full) · [Step 07 输入发布](../../stage2/07_lxpert/README.md#expression-input-release)

## 论文中实际用的输入和目标

| 路线 | 编码器输入和解码目标 | 论文对照及适用范围 |
|---|---|---|
| **SCimilarity**，Nature 2024 | 固定 **28,231** 基因的表达向量，用全连接编码器映射到 **128** 维，再用全连接解码器重建。MSE 是辅助目标，同时使用监督 triplet loss。 | 提供全轴 MLP 自编码的架构先例。论文的细胞相似搜索结果与 Step 07 重建分数属于不同任务。[原文 Methods](https://www.nature.com/articles/s41586-024-08411-y) |
| **xTrimoGene / scFoundation**，NeurIPS 2023 / Nature Methods 2024 | 固定 **19,264** 基因参考轴。xTrimoGene 同时抽取零值和非零值遮盖，重编码器只处理可见非零基因，轻量解码器结合可见表示及其余基因位置，用遮盖位置 MSE 训练。scFoundation 加入读深条件，做表达增强。 | scFoundation 的增强对照包括 SAVER、scImpute 和 MAGIC，服务于低深度输入恢复。可借鉴固定基因身份、连续值嵌入和非对称编码解码。[xTrimoGene §3](https://proceedings.neurips.cc/paper_files/paper/2023/file/db68f1c25678f72561ab7c97ce15d912-Paper-Conference.pdf) · [scFoundation 论文](https://doi.org/10.1038/s41592-024-02305-7) · [增强实验代码说明](https://github.com/biomap-research/scFoundation#read-depth-enhancement) |
| **scGPT**，Nature Methods 2024 | 基因身份和分箱表达，学习从可见基因上下文预测遮盖表达（GEP）以及由细胞嵌入辅助预测（GEPC）。 | Kedzierska 等评估预训练目标时给 scGPT 输入未遮盖基因，用逐基因表达均值作对照，发现 scGPT 的上下文预测可输给该均值。这个先例支持加入无输入的逐基因均值下限。其评分使用分箱表达。[scGPT 原文](https://www.nature.com/articles/s41592-024-02201-0) · [零样本研究 Fig. 4](https://doi.org/10.1186/s13059-025-03574-x) |
| **State Embedding（SE）**，State v2 方法 | 从每细胞的 top-2,048 表达基因构造输入，用细胞嵌入和目标基因身份等信息经 MLP 解码，预训练每细胞抽取 **1,280** 个目标基因，损失兼顾细胞内基因模式与同基因跨细胞变化。 | State 的 **context mean**、**perturbation mean** 用于下游扰动泛化对照，其中一个主要比较在 **2,000 HVG** 评估。SE 的重建表示对照可采用 ReconEval 的同维 PCA。[State 正式论文](https://doi.org/10.1016/j.cell.2026.07.052) · [State v2 §2.2、§4.4](https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2.full) · [ReconEval Fig. 3](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full) |
| **TxFM**，NeurIPS 2026 | 用可见基因的可学习身份和连续表达作为输入，将样本编码到一个向量，再通过 MLP 预测测量基因表达。 | 提供遮盖自编码和全向量解码的近例，论文着重测试表示的迁移能力。[TxFM 原文](https://arxiv.org/html/2605.31562v2) |
| **scVI**，Nature Methods 2018 | 用原始计数做变分推断，论文中的观测模型为 ZINB，显式建模文库大小和批次，并用保留数据似然检查模型拟合。 | 若另做计数重建，可作为计数 VAE 对照。若加入现有图，先把预测变换到与真值相同的 log1p(CP10K) 空间，单独记录其额外使用的计数和协变量。[scVI 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289068/) · [ReconEval 的 VAE 设置](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full) |

## Step 07 可比的评估设计

**主实验采用完整输入的独立重建。** 同一条表达谱输入和作为真值，所有可学习方法在同一训练板上拟合，在保留板用同一表达变换和两源共有的基因轴评分。纳入处理谱和对照谱，无需按 DMSO 可匹配性筛选样本。原发布因 Ensembl ID 差异漏掉了 Tahoe 的 FAM153B，下一版按已核实的跨库映射恢复该基因，实测零保留。此前全量 refit 权重已经见过这些板。可用它做样本内诊断，保留板结果则需用排除了这些板的训练权重，额外划出独立测试板可把调参与最终评估分开。[Step 07 输入与训练](../../stage2/07_lxpert/README.md) · [基因轴决定](../adr/0052-stage2-autoencoder-reconstruction-boundary.md) · [ReconEval 划分与评估](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full)

| 建议对照 | 拟合和推理口径 | 用途 |
|---|---|---|
| 训练集逐基因均值 | 对每个观测基因只从训练样本估计均值。推理输出固定向量。若另报按来源的均值，明确它使用已知来源标签。 | 无输入下限。与旧图的同板 DMSO 加平均 Δ 是两种不同定义。 |
| **PCA → inverse transform** | 只在训练板拟合中心和投影，测试时编码输入谱后反投影。建议扫 32、128、512 个成分，资源允许时增加更高维点。明确分来源拟合或共享拟合，统一基因掩码和聚合口径。 | 主要线性重建对照。512 成分只代表该 PCA 的维数。当前 Step 07 为 512 × 128 个 latent 标量，共 **65,536**，大于输入轴 18,164。它缩短的是 token 序列。对齐表示预算时需固定标量维数或存储开销，并同时报告 token 数。 |
| **全连接 MLP AE** | 在同一 log1p(CP10K) 数据上训练，重建 MSE 只对同一组已测基因计算，报告明确的瓶颈维数和参数规模。可先用与 PCA 相同的几个标量维数。 | 主要非线性对照，可帮助定位固定位置编码和 cross-attention 相比普通 AE 的贡献。 |
| 原样复制 | 完整输入下直接输出输入，PCC 在非恒定样本上为 1，RMSE 为 0。 | 评分管线的正确性检查，单独标注。 |

表中 PCA、确定性 AE 和对齐潜空间维数的思路来自 ReconEval，全连接架构参考 SCimilarity，均值参考 scGPT 的零样本目标评估。这里给出的维数扫点与来源口径是针对 Step 07 的实验建议。[ReconEval](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full) · [SCimilarity](https://www.nature.com/articles/s41586-024-08411-y) · [scGPT 零样本评估](https://doi.org/10.1186/s13059-025-03574-x)

保留逐谱 PCC 和 RMSE，明确按谱等权和恒定向量的处理，并按 Tahoe、tenv 分报。补充跨谱的逐基因 PCC，或相对训练均值的残差误差与重建方差，检查是否只恢复了平均表达轮廓。若另做遮盖，只对共同遮盖的已测位置计分，区分实测零和非零，所有方法看到同一份遮盖输入。该实验可比较逐基因训练均值、带缺失推断的低秩重建和遮盖 AE。若另做读深恢复，固定降深输入与高深真值，单报相应基线。[ReconEval 指标](https://www.biorxiv.org/content/10.64898/2026.06.15.731445v1.full) · [xTrimoGene 遮盖策略](https://proceedings.neurips.cc/paper_files/paper/2023/file/db68f1c25678f72561ab7c97ce15d912-Paper-Conference.pdf) · [scFoundation 增强实验](https://github.com/biomap-research/scFoundation#read-depth-enhancement)

现有图中的 `no_change` 使用同板实测 DMSO，`mean` 使用该 DMSO 加训练平均 Δ。两者都使用了完整表达重建设定之外的样本。它们对应控制到处理的预测问题，不能隔离当前模型的编码解码能力。本地在固定抽取的 96 条 tenv 处理谱上核对过一个更贴近任务的对照。正式 AE 的 PCC/RMSE 为 0.555/0.337，仅使用 3,414 条匹配训练处理谱估计的逐基因均值为 0.825/0.221。这个样本内探针说明，换掉 DMSO 对照后仍有重建差距。正式结论需等同划分、全覆盖的 PCA 和 MLP AE 对照。对架构的优先检验是同数据和同预算的 PCA 与 MLP AE，再决定是否引入可学习基因身份或遮盖训练。[现有图](../../stage2/07_lxpert/output/latent/metrics.png) · [Step 07 的重建与处理分支](../adr/0051-latent-xpert.md)
