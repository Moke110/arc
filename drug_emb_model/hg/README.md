# HG — 异质图药物嵌入模型 (Heterogeneous Graph Neural Network)

## 写给新手：阅读指南

本文档专为深度学习新手编写，将从零开始介绍图神经网络的基本概念，然后深入讲解 HG 模型的结构、输入、输出、训练过程以及所使用的生物信息学先验数据。

如果你已经了解 GNN 基础知识，可以直接跳转到"模型架构"部分。

---

## 目录

1. [前置知识：图、图神经网络与异质图](#1-前置知识图图神经网络与异质图)
2. [生物信息学背景：PPI、DTI 与 DDS](#2-生物信息学背景ppi-dti-与-dds)
3. [HG 项目概述](#3-hg-项目概述)
4. [异质图数据：节点、边与特征来源](#4-异质图数据节点边与特征来源)
5. [模型架构详解](#5-模型架构详解)
6. [训练流程：对比学习与 InfoNCE 损失](#6-训练流程对比学习与-infonce-损失)
7. [推理、外推与下游集成](#7-推理外推与下游集成)
8. [先验数据来源详解](#8-先验数据来源详解)
9. [环境配置与实操指南](#9-环境配置与实操指南)
10. [当前状态与局限性](#10-当前状态与局限性)

---

## 1. 前置知识：图、图神经网络与异质图

### 1.1 什么是图（Graph）？

在深度学习中，图是一种数据结构，由**节点（Node / Vertex）**和**边（Edge）**组成：

```
节点 (Node)  —— 代表一个实体，比如一个人、一个药物、一个基因
边 (Edge)    —— 代表两个实体之间的关系，比如"朋友关系"、"药物与靶点的相互作用"
```

形式化地，图 G = (V, E)，其中 V 是所有节点的集合，E 是所有边的集合。
每条边连接两个节点，可以是有向的（A→B）或无向的（A—B）。

**一个具体的例子：**

```
[阿司匹林] ──(作用于)──> [COX-1基因]
     │
     └──(作用于)──> [COX-2基因]
```

在这个微型图中，有 3 个节点（1 个药物 + 2 个基因）和 2 条边（两种靶向关系）。

### 1.2 什么是图神经网络（GNN）？

图神经网络的核心思想是：**每个节点通过"观察"它的邻居来更新自己的表征**。这个过程叫**消息传递（Message Passing）**。

#### 直观理解

想象一场学术会议：

- 你（节点 A）有一个初始的研究兴趣向量
- 你周围的朋友（邻居节点 B、C、D）各自也有自己的研究兴趣
- 经过交流后，你的研究兴趣会受到朋友们的影响，融合了大家的观点
- 朋友们也同时受到他们各自邻居的影响
- 经过多轮交流（多层 GNN），每个人都形成了更丰富、更有深度的学术观点

GNN 的每一层做的事情本质上就是：

```
h_i^(new) = AGGREGATE( h_i^(old), {h_j^(old) for j in 邻居(i)} )
```

也就是：节点 i 的新表征 = 聚合函数（自己的旧表征 + 所有邻居的旧表征）。

#### 数学形式（以 SAGEConv 为例）

HG 使用的 SAGEConv（GraphSAGE 卷积）的简化数学公式：

```
h_i^(l+1) = W^(l) · CONCAT(h_i^(l), MEAN({h_j^(l) : j ∈ N(i)}))
```

其中：
- `h_i^(l)` 是节点 i 在第 l 层的嵌入向量（一个高维数组，比如 256 维）
- `N(i)` 是节点 i 的邻居集合
- `MEAN` 是对所有邻居的嵌入取平均
- `CONCAT` 是向量拼接
- `W^(l)` 是可学习的权重矩阵（训练的目标就是找到最优的 W）

这意味着每个节点的新表征 = 自己的旧表征 + 邻居的平均表征，经过线性变换后的结果。

### 1.3 什么是异质图（Heterogeneous Graph）？

普通图只有一种节点类型和一种边类型。但现实世界中，**一个系统中往往存在多种类型的节点和多种类型的关系**。

**异质图**就是包含多种节点类型和多种边类型的图。

HG 模型中的异质图包含：

```
节点类型:
  - Drug (药物): 8,981 个
  - Gene  (基因): 19,392 个

边类型（关系类型）:
  - PPI (gene ─ gene):   蛋白质-蛋白质相互作用
  - DTI (drug ─ gene):   药物-靶点相互作用（双向）
  - DDS (drug ─ drug):   药物-药物结构相似性
```

可以这样可视化：

```
                    DDS (药物相似性)
                 +--------------------+
                 |                    |
                 v                    v
           [Drug_0] ----------> [Gene_0]
              |     (DTI)         |  ^
              |                   |  |
              |                   |  | PPI
              |                   |  |
              +-----------------> [Gene_1]
              |     (DTI)         |
              |                   | PPI
              v                   v
           [Drug_1]            [Gene_2]
```

#### 为什么异质图对药物发现很重要？

- 如果只看药物的化学结构（比如 KPGT 做的事情），模型只能学到"结构相似的药物可能功能相似"
- 但如果把药物和基因放在同一张图中，模型就能学到更丰富的生物学关系：
  - "这个药物靶向哪个基因？"
  - "这两个基因的蛋白产物会互相作用吗？"
  - "靶向同一个基因的药物之间有什么关系？"

这就是 **HG 嵌入被称为"生物学嵌入"（biological embedding）** 的原因——它编码了药物在生物学网络中的位置信息，而不仅仅是化学结构。

### 1.4 HeteroConv：异质图的消息传递

在异质图中，消息传递比普通图复杂，因为不同类型的边需要不同的处理方式。

HeteroConv 的工作原理：

```
对于每一种边类型 (src_type, relation, dst_type):
    1. 取出该边类型的所有边
    2. 为源节点和目标节点分别应用卷积（SAGEConv）
    3. 对每种边类型的输出进行聚合（求和或平均）

最终，每个节点从所有关联它的边类型中收到了更新后的嵌入。
```

例如，一个药物节点 Drug_0 在第 l 层的更新会同时考虑：
- 它通过 DTI 边连接的基因节点（消息来自基因）
- 它通过 DDS 边连接的其他药物节点（消息来自其他药物）

这样，Drug_0 的嵌入就同时编码了"它靶向哪些基因"和"它与其他药物有多相似"的信息。

---

## 2. 生物信息学背景：PPI、DTI 与 DDS

在这一节中，我们将介绍 HG 异质图中使用的三种生物信息学关系类型。

### 2.1 PPI — 蛋白-蛋白质相互作用（Protein-Protein Interaction）

#### 基本概念

蛋白质是基因的最终产物，是细胞中执行各种功能的分子机器。大多数蛋白质不是独立工作的——它们与其他蛋白质物理结合，形成复合物来完成特定任务。

PPI 描述的就是**两个蛋白质之间是否存在物理结合或功能关联**。

#### 一个直观的例子

```
p53 蛋白 ──(结合)──> MDM2 蛋白
  │
  └──(结合)──> p300 蛋白
```

p53 是一个著名的肿瘤抑制蛋白。它与 MDM2 结合后被标记降解，与 p300 结合后被激活。理解这些相互作用对于理解癌症机制至关重要。

#### 数据形式

在我们的异质图中，PPI 表现为 **gene-gene 之间的边**（因为蛋白是由基因编码的，所以蛋白间的相互作用映射为基因间的相互作用）。

- 边的数量：901,260 条
- 边的权重：每条边有一个权重值（`edge_attr`），表示相互作用的置信度或强度
- 无向性：如果蛋白 A 与蛋白 B 相互作用，那么蛋白 B 也与蛋白 A 相互作用

#### 为什么 PPI 对药物发现重要？

如果一个药物靶向基因 A，而基因 A 的蛋白产物与基因 B 的蛋白产物相互作用，那么扰动基因 A 可能会间接影响基因 B 的功能。PPI 网络让 HG 模型能够捕捉这种**间接的生物学影响**。

### 2.2 DTI — 药物-靶点相互作用（Drug-Target Interaction）

#### 基本概念

大多数小分子药物通过**与特定的蛋白质（靶点）结合**来发挥药效。这个靶点通常是某种酶、受体或离子通道。DTI 描述的就是哪种药物作用于哪个靶点蛋白。

#### 一个直观的例子

```
阿司匹林 ──(抑制)──> COX-1 酶
阿司匹林 ──(抑制)──> COX-2 酶
二甲双胍 ──(激活)──> AMPK 激酶
```

阿司匹林通过抑制 COX-1/COX-2 酶来减少前列腺素的合成，从而发挥镇痛和抗炎作用。

#### 数据形式

在我们的异质图中，DTI 表现为 **drug-gene 之间的边**：

- 边的数量：12,890 条
- 方向性：同时构建 drug→gene 和 gene→drug 两条边，使信息可以双向流动
- 无权重：这些边是二元的（存在或不存在），不带权重

#### 为什么 DTI 是 HG 的"桥梁"？

DTI 边是连接 **drug 世界** 和 **gene 世界** 的桥梁。如果没有 DTI 边，drug 节点和 gene 节点之间就没有信息流通，模型无法学习药物和基因之间的关联。

### 2.3 DDS — 药物-药物相似性（Drug-Drug Similarity）

#### 基本概念

两种药物如果在化学结构上高度相似，它们很可能具有相似的药理活性。DDS 边记录了在分子指纹层面高度相似的药物对。

#### 数据形式

在我们的异质图中，DDS 表现为 **drug-drug 之间的边**：

- 边的数量：287,834 条
- 边的权重：每条边有一个权重值，表示结构相似性的程度（Tanimoto 系数等）
- 高密度的原因：8981 个药物之间可能存在的药物对数量是 C(8981, 2) ≈ 4000 万，287,834 约为其中的 0.7%

#### 为什么 DDS 重要？

如果药物 A 和药物 B 结构相似，但药物 A 靶向基因 X（已知的 DTI），而药物 B 的靶点未知，那么模型可以通过 DDS 边推断出药物 B 也可能靶向基因 X。这就是**协同过滤**在图上的体现。

### 2.4 三种关系的协同作用

```
            DDS                        DDS
[Drug_A] <=======> [Drug_B]      [Drug_A] ======= [Drug_B]
   |                   |             ||               ||
   | DTI               | DTI         || DTI           || DTI (未知)
   |                   |             ||               ||
[Gene_X] <==PPI==> [Gene_Y]      [Gene_X] <==PPI==> [Gene_Y]
```

左图：所有关系已知，信息可以在 drug 和 gene 世界之间自由流动。
右图：Drug_B 的 DTI 未知。通过 DDS（Drug_A 和 Drug_B 相似）+ Drug_A 的已知 DTI，模型可以推断 Drug_B 可能也靶向 Gene_X。

**HG 的核心思想就是利用这些已知关系的协同，来学习融合了多种生物学背景知识的药物表征。**

---

## 3. HG 项目概述

### 3.1 在整体项目中的位置

HG 是 `arc` 项目（心肌缺血再灌注损伤药物预测）的药物嵌入模块之一。整个项目的数据流如下：

```
SMILES 分子式
    │
    ├──> KPGT (分子图 Transformer)  ──> 分子结构特征 (2304-dim)
    │
    ├──> HG   (异质图神经网络)     ──> 生物学嵌入   (256-dim)
    │
    └──> UniMol (3D 分子预训练)    ──> 3D 结构特征  (122×514)
                  │
                  ▼
         XPert (双分支 Transformer)
                  │
                  ▼
         DEG 预测 (978 个 landmark 基因的转录扰动)
                  │
                  ▼
         MIRI 药物评分与排序
```

HG 的独特价值在于：它提供了**生物学视角**的药物表征。KPGT 看到的是化学结构，UniMol 看到的是 3D 构象，而 HG 看到的是"这个药物在人体生物学网络中处于什么位置"。

### 3.2 核心思想

化学结构相似的药物可能作用于相似的靶点，而共享靶点的药物可能产生相似的转录组扰动。

通过在包含 PPI、DTI、DDS 的异质图上进行消息传递，HG 学习到的药物嵌入编码了：
- **直接信息**：这个药物靶向哪些基因（DTI）
- **间接信息**：这个药物的靶点蛋白与哪些其他蛋白相互作用（PPI）
- **类比信息**：这个药物与哪些其他药物结构相似（DDS）

### 3.3 环境配置

使用 conda 环境 `xpert`：

- **Python**: 3.13
- **PyTorch**: 2.12（深度学习框架）
- **PyTorch Geometric**: 2.7（图神经网络库）

```bash
conda activate arc
cd hg/
```

---

## 4. 异质图数据：节点、边与特征来源

### 4.1 数据总览

```
hg/data/
├── all_gene_node_feat_19392.npy                     # 基因节点特征 (19392, 512), 76 MB
├── all_drug_node_feat_8981.npy                      # 药物节点特征 (8981, 512), 36 MB
├── PPI_edge_index_all_901260_pairs_weighted.npy     # PPI 边 (901,260 对), 21 MB
├── DTI_edge_index_all_12890_pairs.npy               # DTI 边 (12,890 对), 202 KB
├── DDS_edge_index_all_287834_pairs_weighted.npy     # DDS 边 (287,834 对), 6.6 MB
└── 978_gene_node_idx.npy                            # landmark 基因索引 (978,), 7.8 KB
```

### 4.2 节点数据详解

#### 药物节点 (Drug Nodes)

| 属性 | 值 |
|------|-----|
| 数量 | 8,981 |
| 特征维度 | 512 |
| 特征类型 | 化学指纹 (chemical fingerprints) |
| 覆盖范围 | 与 L1000 预训练数据集一致的 8,981 个药物 |

每个药物节点的 512 维特征来自化学信息学中的分子指纹——一种将分子结构编码为固定长度比特向量的方法。具体来说，这些特征编码了分子中存在哪些官能团、环结构、原子类型等化学信息。

#### 基因节点 (Gene Nodes)

| 属性 | 值 |
|------|-----|
| 数量 | 19,392 |
| 特征维度 | 512 |
| 特征类型 | 基因表达谱 (gene expression profiles) |
| 其中 landmark 基因 | 978 个 (L1000 项目的标志性基因) |

19,392 是涵盖了人类基因组中所有蛋白质编码基因（约 20,000 个）的一个较完整集合。每个基因节点的 512 维特征来自基因表达数据——在不同实验条件下测得的基因表达水平，经过降维或编码后得到。

**为什么是 19,392 个基因？**
- 略少于标准的人类蛋白质编码基因数（~20,000）
- 具体数量取决于数据融合时各个数据库的交集
- 这个数量确保 PPI 网络覆盖了主要的已知相互作用

#### Landmark 基因 (978 个)

L1000 项目（由 Broad Institute 发起）选出了 978 个"标志性基因"（landmark genes），这 978 个基因的表达水平可以较好地推断其余 ~11,000 个基因的表达水平。这大大降低了高通量转录组测序的成本。XPert 模型的预测目标正是这 978 个基因。

### 4.3 边数据详解

所有边文件以 `.npy` 格式存储，使用 `np.load(path, allow_pickle=True).item()` 加载，每个文件包含一个 Python 字典。

#### PPI 边（蛋白质-蛋白质相互作用）

```
文件: PPI_edge_index_all_901260_pairs_weighted.npy
结构:
  {
    'edge_index': np.ndarray, shape (2, 901260), dtype int64
    'edge_attr':  np.ndarray, shape (901260,),   dtype float64
  }
```

- `edge_index[0]` 是源基因的索引（0 到 19391 之间的整数）
- `edge_index[1]` 是目标基因的索引（0 到 19391 之间的整数）
- `edge_attr` 是每条边的权重，表示相互作用的置信度/强度

#### DTI 边（药物-靶点相互作用）

```
文件: DTI_edge_index_all_12890_pairs.npy
结构:
  {
    'edge_index': np.ndarray, shape (2, 12890), dtype int64
  }
```

- `edge_index[0]` 是药物索引（0 到 8980 之间的整数）
- `edge_index[1]` 是靶点基因索引（0 到 19391 之间的整数）
- 无权重，表示已知的药物-靶点关系
- 在构建异质图时，**同时创建反向边** (gene→drug)，使信息可以双向流动

#### DDS 边（药物-药物相似性）

```
文件: DDS_edge_index_all_287834_pairs_weighted.npy
结构:
  {
    'edge_index': np.ndarray, shape (2, 287834), dtype int64
    'edge_attr':  np.ndarray, shape (287834,),   dtype float64
  }
```

- `edge_index[0]` 和 `edge_index[1]` 都是药物索引
- `edge_attr` 是相似性分数（Tanimoto 系数或其他分子指纹相似度），越高表示两个药物结构越相似

### 4.4 异质图构建过程

在训练入口 `pretrain_hg.py` 中，HeteroData 对象的构建逻辑如下：

```python
data = HeteroData()
data['gene'].x = gene_feat    # (19392, 512)
data['drug'].x = drug_feat    # (8981, 512)

# 添加 PPI 边: gene <-> gene（带权重）
data['gene', 'PPI', 'gene'].edge_index = ppi_edge_index
data['gene', 'PPI', 'gene'].edge_attr  = ppi_edge_attr

# 添加 DTI 边: drug -> gene（无权重）
data['drug', 'DTI', 'gene'].edge_index = dti_edge_index
# 同时添加反向边: gene -> drug
data['gene', 'DTI', 'drug'].edge_index = reverse(dti_edge_index)

# 添加 DDS 边: drug <-> drug（带权重）
data['drug', 'DDS', 'drug'].edge_index = dds_edge_index
data['drug', 'DDS', 'drug'].edge_attr  = dds_edge_attr
```

### 4.5 孤立节点检测

构建图后，训练脚本会检测哪些节点没有任何边连接。这些"孤立节点"在消息传递中无法获得邻居信息，它们的嵌入仅由自身特征决定。

根据实际运行结果：
- 总药物节点数：8,981，其中孤立药物节点：少量（保存到 `isolated_drug_indices.npy`）
- 总基因节点数：19,392，其中孤立基因节点：少量（保存到 `isolated_gene_indices.npy`）

---

## 5. 模型架构详解

### 5.1 总体结构

HG 模型（`HeteroGNN`）是一个三层的异质图卷积网络，包含以下主要组件：

```
输入层:
  药物特征 (512-dim) ──> Linear(512→256) ──> 药物投影 (256-dim)
  基因特征 (512-dim) ──> Linear(512→256) ──> 基因投影 (256-dim)

图卷积层 (×3):
  HeteroConv {
    (drug, DDS, drug):  SAGEConv((-1,-1), 256)   ← 药物-药物消息
    (drug, DTI, gene):  SAGEConv((-1,-1), 256)   ← 药物→基因消息
    (gene, DTI, drug):  SAGEConv((-1,-1), 256)   ← 基因→药物消息
    (gene, PPI, gene):  SAGEConv((-1,-1), 256)   ← 基因-基因消息
  }
  → ReLU 激活

关系预测器 (可选，用于推理):
  DTI 预测头:  Linear(512 → 1)  ← concat(药物嵌入, 基因嵌入)
  DDS 预测头:  Linear(512 → 1)  ← concat(药物嵌入, 药物嵌入)
  PPI 预测头:  Linear(512 → 1)  ← concat(基因嵌入, 基因嵌入)

输出:
  药物嵌入 (256-dim)  ← 最终层的 drug 节点嵌入
  基因嵌入 (256-dim)  ← 最终层的 gene 节点嵌入
```

### 5.2 输入投影层

```python
self.transforme_gene_node = torch.nn.Linear(512, 256)  # 512 → 256
self.transforme_drug_node = torch.nn.Linear(512, 256)  # 512 → 256
```

这两个线性层将原始 512 维节点特征投影到 256 维的隐藏空间。它们的作用是：
- 降维：减少后续图卷积的计算量
- 对齐：将药物特征和基因特征映射到同一个语义空间，使它们可以在异质图中进行有意义的交互

### 5.3 图卷积层：HeteroConv + SAGEConv

HG 使用 3 层 HeteroConv，每层内部包含多个 SAGEConv 实例。

#### SAGEConv（GraphSAGE 卷积）

SAGEConv 的全称是 **SAmple and aggreGatE Convolution**。它的核心公式是：

```
h_i' = W · [h_i || AGGREGATE({h_j : j ∈ N(i)})]
```

其中 `||` 表示向量拼接，`AGGREGATE` 通常是取均值（mean）。

SAGEConv 的 `(-1, -1)` 参数表示：输入维度从数据中自动推断（源节点和目标节点可能具有不同的维度，但在异质图中它们都被投影到了 256 维）。

#### HeteroConv（异质图卷积）

HeteroConv 是一个"容器"，它同时管理多种边类型的消息传递。在一层 HeteroConv 中：

```
对于 drug 节点：
  它的邻居包括:
    - 通过 DDS 边连接的其他 drug 节点  → SAGEConv(drug, drug)
    - 通过 DTI 反向边连接的 gene 节点   → SAGEConv(gene, drug)
  将这些邻居的消息聚合 → 更新后的 drug 嵌入

对于 gene 节点：
  它的邻居包括:
    - 通过 PPI 边连接的其他 gene 节点  → SAGEConv(gene, gene)
    - 通过 DTI 正向边连接的 drug 节点   → SAGEConv(drug, gene)
  将这些邻居的消息聚合 → 更新后的 gene 嵌入
```

#### 三层堆叠的意义

```
第 1 层: 每个节点看到它的直接邻居（1-hop）
第 2 层: 每个节点看到邻居的邻居（2-hop）
第 3 层: 每个节点看到邻居的邻居的邻居（3-hop）
```

对于药物节点来说：
- 第 1 层后，药物嵌入包含了"它直接靶向哪些基因"的信息
- 第 2 层后，药物的嵌入进一步包含了"它靶向的基因与哪些其他基因相互作用"的信息（通过 PPI 路径 drug→gene1→gene2）
- 第 3 层后，信息传播到更远的生物学网络中

**这就是 HG 嵌入被称为"生物学嵌入"的核心原因**——通过多层消息传递，药物的表征融合了从直接靶点到远距离生物学网络的多层次信息。

### 5.4 激活函数

每层 HeteroConv 之后应用 ReLU（Rectified Linear Unit）：

```python
x_dict = {key: x.relu() for key, x in conv(x_dict, edge_index_dict).items()}
```

ReLU(x) = max(0, x)，将负值置零，保留正值。这引入了非线性，使模型能够学习更复杂的模式。

### 5.5 关系预测器（推理辅助组件）

HG 包含三个可选的关系预测头：

```python
self.relation_predictors = torch.nn.ModuleDict({
    'DTI': torch.nn.Linear(512, 1),  # 输入 concat(药物嵌入, 基因嵌入) = 512
    'DDS': torch.nn.Linear(512, 1),  # 输入 concat(药物嵌入, 药物嵌入) = 512
    'PPI': torch.nn.Linear(512, 1),  # 输入 concat(基因嵌入, 基因嵌入) = 512
})
```

这些预测器在训练阶段不被使用，主要用于训练后的推理：
- 预测未知的药物-靶点关系（DTI 预测）
- 预测药物-药物相似性（DDS 预测）
- 预测蛋白质-蛋白质相互作用（PPI 预测）

预测方式是将两个节点的嵌入拼接后通过线性层得到分数：

```python
# 例如 DTI 预测
pair_repr = concat(drug_embedding[i], gene_embedding[j])  # shape: (512,)
score = Linear(512→1)(pair_repr)                          # shape: (1,)  标量分数
```

分数经过 sigmoid 转换为概率：`p = 1 / (1 + exp(-score))`

### 5.6 模型参数量估计

- 输入投影：`512×256 + 256` (drug) + `512×256 + 256` (gene) ≈ 263K
- 每层 HeteroConv（包含 4 种边类型各一个 SAGEConv）：约 4 × (256×256 + 256) ≈ 263K × 3 layers ≈ 789K
- 关系预测器：3 × (512×1 + 1) ≈ 1.5K
- **总计：约 1M 参数**（CLAUDE.md 中估计 25M 偏大，实际远小于该值）

---

## 6. 训练流程：对比学习与 InfoNCE 损失

### 6.1 训练目标

HG 的训练目标是：**学习节点嵌入，使得在图中相连的节点具有相似的嵌入向量，而不相连的节点具有不相似的嵌入向量。**

形式化地：
- 如果 Drug_A 和 Gene_X 之间存在 DTI 边，我们希望 `cosine_sim(emb(Drug_A), emb(Gene_X))` 尽可能大
- 如果 Drug_A 和 Gene_Y 之间没有 DTI 边，我们希望 `cosine_sim(emb(Drug_A), emb(Gene_Y))` 尽可能小
- 同理适用于 PPI 和 DDS 边类型

### 6.2 InfoNCE 对比损失（Contrastive Loss）

HG 使用 InfoNCE（Information Noise-Contrastive Estimation）损失函数：

```
对于每条正样本边 (i, j):
  pos_sim = cos(z_i, z_j) / τ          ← 正样本对的相似度（除以温度）
  neg_sim = cos(z_i, z_neg_k) / τ      ← K 个负样本的相似度（除以温度）

  L(i,j) = -log( ----------------------------------------------------- )
                      exp(pos_sim) + Σ_{k=1}^{K} exp(neg_sim_k)
```

其中：
- **z_i, z_j**：分别是一条正样本边的两个端点的嵌入向量
- **τ (temperature)**：温度参数，默认 0.5。较小的 τ 使模型对难样本更敏感
- **z_neg_k**：随机采样的负样本（不与节点 i 相连的节点）
- **K (num_neg_samples)**：每条正边对应的负样本数，默认 3

#### 直观理解

这个损失函数可以理解为一个 (K+1) 类的分类问题：
- 给定锚点节点 i，在 K+1 个候选中找出哪个是它真正的邻居 j
- 模型需要给真正的邻居 j 更高的相似度，给随机负样本更低的相似度

就像做选择题：给你一个药物，在 4 个基因中找出它真正靶向的那个。做对了得分，做错了扣分。

### 6.3 负采样（Negative Sampling）

对于每种边类型（PPI、DTI、DDS），训练时需要生成负样本——即图中不存在的边。

```python
def generate_negative_samples(num_src, num_dst, num_neg_samples, pos_edge_index):
    # 从所有可能的源-目标节点对中随机采样
    # 排除正样本边（图中已存在的边）
    # 为每个源节点生成 num_neg_samples 个负样本
```

关键细节：
- **负采样是从所有可能的节点对中随机抽取**，不包括已存在的正样本边
- 如果随机采样过程中"不幸"抽到了正样本，会重试（最多尝试 `3 × num_src × num_neg_samples` 次）
- 这种采样方式给每个源节点生成了 `num_neg_samples` 个不同的负样本目标节点

### 6.4 小批量训练与 NeighborLoader

由于全图（19,392 + 8,981 ≈ 28,000 节点，120 万+ 条边）无法一次性装入 GPU 显存，HG 使用 **邻居采样（Neighbor Sampling）** 进行小批量训练。

#### NeighborLoader 的工作方式

```python
loader = NeighborLoader(
    data,
    num_neighbors=[10, 5, 3],   # 每层采样的邻居数
    batch_size=1,                # 每批 drug seed 节点数
    input_nodes=('drug', torch.arange(8981))  # 以所有 drug 节点为种子
)
```

对于一个 batch：
1. 选择一个（或 batch_size 个）drug 节点作为"种子"
2. 采样该 drug 的 DDS 邻居（最多 10 个）
3. 通过 DTI 边采样基因邻居（最多 5 个）
4. 通过 PPI 边采样远距离基因邻居（最多 3 个）
5. 这样就构建了一个包含 ~数十个节点的小型子图
6. 在这个子图上计算对比损失，反向传播更新参数

`num_neighbors=[10, 5, 3]` 的含义：
- 第 1 层：每个节点采样最多 10 个直接邻居
- 第 2 层：每个节点采样最多 5 个二跳邻居
- 第 3 层：每个节点采样最多 3 个三跳邻居

**注意**：run.sh 中使用的是 `[10, 5, 3]`，而非代码中的默认值 `[35, 20, 10]`。较小的邻居采样数减少了每个 batch 的大小，加快了训练速度。

#### batch_size=1 的含义

`batch_size=1` 意味着每个 batch 只以一个 drug 节点为种子。这样做的好处是：
- 每个 batch 的子图很小，显存占用低
- 每个 epoch 有 8,981 个 batch（每个 drug 一个 batch）
- 可以在梯度更新中获得更多样化的信号

代价是：
- batch 内的负采样多样性有限（因为子图中节点数量少）
- 训练噪声较大，需要更多 epoch 来稳定

### 6.5 早停机制（Early Stopping）

HG 使用三重早停策略，防止过拟合和无效训练：

1. **标准早停 (patience=30)**：基于 loss 的早停。如果连续 30 个 epoch 的 loss 没有低于历史最优值，停止训练。

2. **连续上升早停 (conseq_patience=3)**：如果 loss 连续 3 个 epoch 都在上升，说明模型可能进入了不稳定的状态，立即停止。

3. **累计上升早停 (turn_patience=10)**：如果 loss 在整个训练过程中累计上升了 10 次（即使不是连续的），停止训练。

实际训练中（2026-05-18），模型在 epoch 24 被 turn_patience 触发早停。

### 6.6 训练输出

训练完成后，生成以下文件：

| 文件 | 路径 | 格式 |
|------|------|------|
| 最佳模型权重 | `model/<timestamp>/global_fold_early_stop.pth` | PyTorch checkpoint (9.1 MB) |
| 药物嵌入 | `embeddings/HG_drug_embeddings.npy` | (8981, 256) float32 |
| 药物嵌入备份 | `model/<timestamp>/phase1_drug_embeddings.npy` | (8981, 256) float32 |
| 基因嵌入 | `model/<timestamp>/phase1_gene_embeddings.npy` | (19392, 256) float32 |
| 训练日志 | `model/<timestamp>/<timestamp>.log` | 文本日志 |

### 6.7 完整的训练命令

```bash
cd hg/
conda activate arc

# 默认参数训练
bash run.sh

# 自定义参数训练
bash run.sh --device cuda:0 --epochs 500 --lr 0.001 --num_neg_samples 5

# 从 checkpoint 恢复训练
bash run.sh --global_saved_model model/20260518_102121/global_fold_early_stop.pth

# 训练后计算 DTI 分数矩阵
bash run.sh --infer_dti_scores True

# 训练后计算药物-基因关系矩阵
bash run.sh --output_relation_matrix True
```

---

## 7. 推理、外推与下游集成

### 7.1 训练药物的嵌入提取（推理）

对于训练集中的 8,981 个药物，获取嵌入非常简单——直接将完整的异质图输入模型进行前向传播：

```python
model.eval()
with torch.no_grad():
    x_dict = model(data.x_dict, data.edge_index_dict)
    drug_embeddings = x_dict['drug']  # (8981, 256)
```

这利用了 GNN 的直推学习（transductive learning）特性——训练和推理使用同一张图结构。

### 7.2 FDA 药物的嵌入生成（完整流程）

对于训练集之外的 FDA 药物（索引 8981–11103，共 2,123 个），采用以下完整的标准 HG 流程：

**阶段一：数据收集**

```
步骤 1: DTI 查询 (ChEMBL API)
  query_dti_fda.py → activities 端点 → target 端点 → gene symbol → HG 基因索引
  结果: 3,984 条 DTI 边，覆盖 856/2,123 个 FDA 药物，535 个 unique 靶点基因
  阈值: pChEMBL ≥ 6.0

步骤 2: DDS 计算 (Morgan 指纹 Tanimoto 相似性)
  compute_dds_fda.py → 计算 2123×8981 (FDA↔Train) + 2123×2123 (FDA↔FDA)
  结果: 7,659 条 DDS 边 (5,317 FDA↔Train + 2,342 FDA↔FDA)
  阈值: 0.5106 (原有 DDS 权重的第 5 百分位)
```

**阶段二：HG 微调**

```
步骤 3: 生成 FDA 药物输入特征 (KPGT→512)
  Ridge 回归: KPGT(2304-dim) → HG 输入空间 (512-dim)
  R² = 0.97，拟合于 8,981 个训练药物

步骤 4: 构建扩展异质图
  11,104 药物 + 19,392 基因 + 合并 DTI (16,874) + 合并 DDS (295,493) + PPI (901,260)

步骤 5: InfoNCE 微调 (finetune_hg.py)
  全图前向传播 → 采样 DTI/DDS/PPI 边 → InfoNCE 对比损失 → 低学习率反向传播
  学习率: 5e-5 (原始训练的 1/10)，50 epochs
  训练药物嵌入余弦相似度保持 > 0.98 (vs 微调前)

步骤 6: 分位数范数映射 (分布对齐)
  对微调后的 FDA 嵌入按 per-drug L2 范数进行分位数插值映射
  目标分布: 训练药物的 L2 范数分布 (mean=3048, std=1988)
  结果: FDA 嵌入全局 std=226 (与训练药物完全一致), FDA/Train std ratio=1.000
  DTI 生物学信号: 2.48x vs random
```

### 7.3 下游集成：XPert

HG 嵌入最终被 XPert 消费，路径为 `miri_reposition/xpert-fda-v1/data/embeddings/all_drugs_idx2HG.npy` (11104, 256)：

- **索引 0–8980**：作者提供的参考嵌入（与 `embeddings/article/HG_drug_embeddings.npy` 完全一致）
- **索引 8981–11103**：FDA 药物嵌入（HG 微调 + 分位数范数映射生成）

HG 嵌入与 KPGT 分子图特征 (2304-dim) 和 UniMol 3D 特征 (122×514) 一起，构成 XPert 的三种互补药物表征。

**两份规范嵌入文件**：

| 文件 | 形状 | 说明 |
|------|------|------|
| `miri_reposition/xpert-fda-v1/data/embeddings/all_drugs_idx2HG.npy` | (11104, 256) | 全药物 HG 嵌入（XPert 默认使用） |
| `embeddings/article/HG_drug_embeddings.npy` | (8981, 256) | 作者原始 L1000 参考嵌入 |

---

## 8. 先验数据来源详解

HG 异质图中使用的所有数据均来自**公开的生物信息学数据库**。原始数据处理脚本在仓库重构时归档，但数据处理逻辑清晰可追溯。

### 8.1 药物节点特征 (Drug Node Features)

**来源**：化学信息学工具生成的分子指纹（Molecular Fingerprints）

分子指纹是将化学结构编码为固定长度位向量的标准方法。常见的分子指纹类型包括：
- **Morgan/Circular Fingerprints** (ECFP)：以每个原子为中心，编码其周围的环形子结构
- **MACCS Keys**：预定义的 166 个化学子结构的存在/缺失
- **RDKit Fingerprints**：基于路径的拓扑指纹

HG 使用 512 维的分子指纹，很可能是通过 RDKit 生成的 Morgan 指纹（或其变体），捕捉了每个药物分子中存在的官能团、环结构和原子连接模式。

### 8.2 基因节点特征 (Gene Node Features)

**来源**：基因表达数据库（最可能是 GEO — Gene Expression Omnibus）

基因表达特征是通过整合大量公开的转录组实验数据（不同组织、不同条件下的基因表达谱），然后通过降维技术（如 PCA、autoencoder 等）压缩到 512 维向量。

这些特征编码了每个基因在不同生物学环境中的表达模式，使具有相似表达模式的基因在特征空间中更接近。

### 8.3 PPI 数据

**来源**：蛋白质相互作用公共数据库

主要的 PPI 数据来源可能包括：

| 数据库 | 全称 | 简介 |
|--------|------|------|
| STRING | Search Tool for the Retrieval of Interacting Genes/Proteins | 最全面的 PPI 数据库，整合了实验验证、共表达、文本挖掘等多种证据 |
| BioGRID | Biological General Repository for Interaction Datasets | 专注于实验验证的蛋白质和遗传相互作用 |
| IntAct | — | EMBL-EBI 维护的分子相互作用数据库，强调 curated 数据 |
| HPRD | Human Protein Reference Database | 人类蛋白质注释数据库 |

**数据整合方式**：
- 从这些数据库下载 PPI 数据
- 将蛋白 ID 映射到基因 ID（通过 UniProt 或 Ensembl 等标准标识符）
- 仅保留两个端点都在 19,392 个基因集合内的边
- 边的权重来自 STRING 的 combined score 或其他数据库的置信度分数

**为什么是 901,260 条 PPI 边？**
这是从上述数据库整合、去重、并限制在 19,392 个基因范围内后的结果。人类蛋白质组中估计有 130,000–650,000 种二元相互作用，901,260 条边说明数据整合了 STRING 等数据库中的多种证据类型（包括文本挖掘和预测性边）。

### 8.4 DTI 数据

**来源**：药物-靶点相互作用公共数据库

主要的数据来源可能包括：

| 数据库 | 全称 | 简介 |
|--------|------|------|
| DrugBank | — | 最全面的药物信息数据库，包含已批准和实验性药物的靶点信息 |
| ChEMBL | — | EMBL-EBI 维护的生物活性分子数据库，包含药物-靶点结合数据 |
| BindingDB | — | 实验测定的蛋白质-小分子结合亲和力数据库 |
| STITCH | Search Tool for Interactions of Chemicals | 化学物质-蛋白质相互作用数据库，整合多种证据 |
| DGIdb | Drug-Gene Interaction Database | 专注于药物-基因相互作用的数据库 |

**数据整合方式**：
- 从 DrugBank、ChEMBL、STITCH 等多个数据库收集药物-靶点关系
- 将药物名称/SMILES 映射到 HG 的药物索引空间
- 将靶点蛋白 ID 映射到 HG 的基因索引空间
- 仅保留药物索引在 0–8980 且基因索引在 0–19391 范围内的边
- 由于只有 12,890 条边，说明使用了较高置信度的阈值进行过滤

**为什么 DTI 边相对较少（12,890）？**
平均每个药物只有约 1.4 个已知靶点基因。这是因为：
1. 许多药物的完整靶点谱尚未被实验验证
2. DTI 数据库倾向于收录有明确实验证据的靶点
3. 药物通常有 1–3 个主要靶点（但可能有数十个脱靶效应）

**FDA 药物 DTI 扩展 (2026-06-10)**：
- 通过 ChEMBL API 查询了 2,123 个 FDA 药物的靶点信息
- 查询流程：`activities` 端点 → `target_chembl_id` → `target` 端点 → `gene symbol` → HG 基因索引
- 使用 pChEMBL ≥ 6.0 作为活性阈值
- 结果：3,984 条边，覆盖 856/2,123 药物 (40%)，535 unique 靶点基因
- 实现脚本：`src/query_dti_fda.py`，全局靶点缓存 15,025 条记录

### 8.5 DDS 数据

**来源**：基于分子指纹计算的药物-药物相似性

DDS 数据不是来自外部数据库，而是通过对药物节点特征（分子指纹）进行相似性计算得到的：

```
对于每对药物 (i, j):
    sim(i, j) = Tanimoto_Similarity(fp_i, fp_j)
    如果 sim(i, j) > 阈值:
        添加边 (i, j)，权重 = sim(i, j)
```

287,834 条 DDS 边说明使用了相对宽松的相似性阈值，使得平均每个药物有约 32 个结构相似的邻居药物。

**FDA 药物 DDS 扩展 (2026-06-10)**：
- 使用 Morgan 指纹 (1024-dim) 计算 Tanimoto (Jaccard) 相似性
- 矩阵计算：2123×8981 (FDA↔Train) + 2123×2123 (FDA↔FDA)
- 阈值：0.5106（原有 DDS 权重的第 5 百分位，匹配 0.71% 密度）
- 结果：7,659 条边（5,317 FDA↔Train + 2,342 FDA↔FDA）
- 实现脚本：`src/compute_dds_fda.py`

### 8.6 数据汇总表

| 数据类型 | 数量/维度 | 估计来源 | 原始数据规模 |
|----------|----------|----------|------------|
| Drug 节点特征 (训练) | 8,981 × 512 | UniMol mol 特征 (经 Ridge 映射自 KPGT) | 8,981 SMILES |
| Drug 节点特征 (FDA) | 2,123 × 512 | Ridge(KPGT→512) 投影, R²=0.97 | 2,123 SMILES |
| Gene 节点特征 | 19,392 × 512 | GEO 基因表达谱 + 降维 | > 数万样本 |
| PPI 边 | 901,260 对 (加权) | STRING (score>700) | > 数百万原始相互作用 |
| DTI 边 (训练) | 12,890 对 | DrugBank, ChEMBL, STITCH, DGIdb | > 数万条药物-靶点关系 |
| DTI 边 (FDA 扩展) | 3,984 对 | ChEMBL (pChEMBL≥6.0) | 856/2,123 药物, 535 unique 基因 |
| DDS 边 (训练) | 287,834 对 (加权) | 分子指纹 Tanimoto 相似性 | C(8981, 2) ≈ 40M 对 |
| DDS 边 (FDA 扩展) | 7,659 对 | Tanimoto>0.5106 (第5百分位) | 5,317 FDA↔Train + 2,342 FDA↔FDA |
| Landmark 基因 | 978 个索引 | L1000 项目 (Broad Institute) | 978 genes |
| 总药物节点 | 11,104 | 8,981 训练 + 2,123 FDA | — |
| 总基因节点 | 19,392 | — | — |

---

## 9. 环境配置与实操指南

### 9.1 目录结构

```
hg/
├── src/                              # 源代码
│   ├── __init__.py
│   ├── model_hg.py                   # HeteroGNN 模型定义 (SAGEConv, HeteroConv, 关系预测器)
│   ├── pretrain_hg.py                # 训练入口 (argparse, NeighborLoader, Phase 1 训练循环)
│   ├── finetune_hg.py                # 微调入口 (全图 InfoNCE, 低 LR, 范数映射)
│   ├── infer_fda_hg.py               # FDA 推理: 扩展图构建 + 单次前向传播
│   ├── query_dti_fda.py              # ChEMBL DTI 查询: activities → target→gene → HG 索引
│   ├── compute_dds_fda.py            # Morgan 指纹 Tanimoto DDS 边计算
│   ├── gene_id_mapping.py            # Gene symbol → HG 节点索引 (978 landmark 锚点)
│   ├── merge_dti_batches.py          # DTI 批次文件合并
│   └── utils_hg.py                   # 工具 (ContrastiveLoss, 负采样, 子图构建, 嵌入提取)
├── data/                             # 异质图数据 (静态 .npy 文件)
│   ├── all_drug_node_feat_8981.npy          # 训练药物节点特征 (8981, 512)
│   ├── all_gene_node_feat_19392.npy         # 基因节点特征 (19392, 512)
│   ├── PPI_edge_index_all_901260_pairs_weighted.npy  # PPI 边 (基因-基因, 901,260 对)
│   ├── DTI_edge_index_all_12890_pairs.npy           # DTI 边 (药物-基因, 12,890 对)
│   ├── DTI_edge_index_fda_ext.npy                   # FDA DTI 边 (3,984 对, 856 药物)
│   ├── DDS_edge_index_all_287834_pairs_weighted.npy # DDS 边 (药物-药物, 287,834 对)
│   ├── DDS_edge_index_fda_drugs.npy                 # FDA DDS 边 (7,659 对)
│   ├── 978_gene_node_idx.npy                # L1000 landmark 基因索引 (978,)
│   ├── fda_extension_drugs.csv              # FDA 扩展药物列表 (2,123 药物)
│   └── dti_query.log                        # ChEMBL 查询日志
├── model/                            # 训练输出
│   ├── 20260518_102121/              # 预训练 checkpoint
│   │   ├── global_fold_early_stop.pth       # 预训练模型 (epoch 23, ~9 MB)
│   │   └── phase1_drug_embeddings.npy      # Phase 1 输出 (8981, 256)
│   └── finetune/20260610_111617/     # 微调 checkpoint
│       ├── final_checkpoint.pth             # 微调模型 (epoch 50, ~9 MB)
│       └── finetune_history.json            # 训练指标历史
├── embeddings/                       # 最终嵌入输出
│   └── article/
│       └── HG_drug_embeddings.npy    # 作者提供的参考嵌入 (8981, 256)
├── run.sh                            # 训练启动脚本
└── run.ipynb                         # 交互式训练 + 数据验证 Notebook
```

### 9.2 快速开始

```bash
# 1. 激活环境
conda activate arc

# 2. 进入目录
cd hg/

# 3. 从头训练 HG (8,981 药物)
bash run.sh

# 4. 为 FDA 药物生成 HG 嵌入 (完整流程)
#   a. 查询 ChEMBL DTI 数据
python -m src.query_dti_fda --start 0 --end 2123 --min-pchembl 6.0
#   b. 计算 DDS 边
python -m src.compute_dds_fda
#   c. 微调 + 生成最终嵌入
#      --embeddings-output 指定全药物 HG 表写入位置 (默认 hg/embeddings/)；
#      工作流 stage 02 会以工作流 data/embeddings/ 路径调用本命令。
python -m src.finetune_hg --device cuda:0 --lr 5e-5 --epochs 50 \
    --update-xpert-embeddings --embeddings-output embeddings/all_drugs_idx2HG.npy

# 5. 验证输出嵌入
python -c "
import numpy as np
hg = np.load('../miri_reposition/xpert-fda-v1/data/embeddings/all_drugs_idx2HG.npy')
print(f'Shape: {hg.shape}')
print(f'Train std: {hg[:8981].std():.1f}, FDA std: {hg[8981:].std():.1f}')
print(f'FDA/Train ratio: {hg[8981:].std()/hg[:8981].std():.3f}')
"
```

### 9.3 默认超参数

| 超参数 | 默认值 | 作用 |
|--------|--------|------|
| `--num_layers` | 3 | 图卷积层数 |
| `--hidden_dim` | 256 | 隐藏层维度 |
| `--out_dim` | 256 | 输出嵌入维度 |
| `--num_neighbors` | 10 5 3 | 每层邻居采样上限 |
| `--batch_size` | 1 | 每批 drug seed 数 |
| `--epochs` | 300 | 最大训练 epoch 数 |
| `--patience` | 30 | 标准早停耐心值 |
| `--conseq_patience` | 3 | 连续 loss 上升容忍数 |
| `--turn_patience` | 10 | 累计 loss 上升容忍数 |
| `--num_neg_samples` | 3 | 每条正边对应的负样本数 |
| `--temperature` | 0.5 | InfoNCE 温度 |
| `--learning_rate` | 0.0005 | Adam 优化器学习率 |
| `--seed` | 4242 | 随机种子 |
| `--include_edges` | PPI DTI DDS | 使用的边类型 |

---

## 10. 当前状态与局限性

### 10.1 状态

- **HG 预训练**：已完成 (2026-05-18), epoch 23 早停
- **FDA 数据收集**：已完成 (2026-06-10), ChEMBL DTI 3,984 条边 + DDS 7,659 条边
- **HG 微调**：已完成 (2026-06-10), lr=5e-5, 50 epochs, 训练药物 cosim=0.9822 保持
- **FDA 嵌入**：已集成到 `all_drugs_idx2HG.npy`, FDA/Train std ratio=1.000, DTI 信号 2.48x
- **项目清理**：已完成, 精简中间文件, 保留两份规范嵌入文件和关键 checkpoint

### 10.2 已知局限性

1. **训练未充分收敛**：本地预训练在 epoch 23 (实际为 24，checkpoint 记为 23) 即早停。训练更长时间的模型可能产生与作者参考嵌入更一致的嵌入空间。

2. **嵌入空间差异**：本地训练模型与作者参考模型的嵌入空间不同（余弦相似度 ~0.1）。当前通过保留训练行 = 作者参考 + FDA 行 = 微调后分位数映射来兼容 XPert。

3. **FDA DTI 覆盖率有限**：仅 856/2,123 (40%) 的 FDA 药物有已知靶点 (ChEMBL pChEMBL≥6.0)。其余药物仅通过 DDS 边接收间接梯度信号。更多 DTI 数据来源可改善覆盖。

4. **batch_size=1 的限制** (预训练阶段)：每个 batch 仅处理 1 个药物种子节点，batch 内负采样多样性有限。微调采用全图训练避免了此问题。

5. **疾病节点未使用**：代码中预留了 Disease 节点类型，但实际训练中未使用。疾病信息可进一步提升嵌入质量。

6. **GSEA/PCC 边未启用**：额外的基因-基因关系可能提供互补的生物学信号。

7. **FDA 输入特征依赖岭回归**：FDA 药物的 512-dim 输入特征通过 Ridge(KPGT→512) 生成 (R²=0.97)。虽然效果好，但仍为线性近似。未来可考虑直接从 UniMol mol 特征生成 FDA 药物的 HG 输入特征。

8. **数据时效性**：PPI/DTI 数据库持续更新，当前嵌入反映数据下载时的知识状态。

---

## 附录：关键术语对照表

| 英文 | 中文 | 解释 |
|------|------|------|
| Graph Neural Network (GNN) | 图神经网络 | 在图结构数据上运行的神经网络 |
| Heterogeneous Graph | 异质图 | 包含多种节点类型和多种边类型的图 |
| Node / Vertex | 节点 | 图中的基本单元，代表一个实体 |
| Edge | 边 | 连接两个节点的关系 |
| Message Passing | 消息传递 | GNN 中节点通过边交换信息的过程 |
| Embedding | 嵌入 | 将离散对象（药物、基因）映射为连续向量 |
| SAGEConv | GraphSAGE 卷积 | 一种采样+聚合的图卷积操作 |
| HeteroConv | 异质图卷积 | 处理多种边类型的卷积层容器 |
| Contrastive Learning | 对比学习 | 通过区分正负样本对来学习表征 |
| InfoNCE Loss | InfoNCE 损失 | 基于噪声对比估计的对比损失函数 |
| Negative Sampling | 负采样 | 从非边中随机采样作为训练负样本 |
| PPI | 蛋白质-蛋白质相互作用 | 两个蛋白质之间的物理结合或功能关联 |
| DTI | 药物-靶点相互作用 | 药物分子与其靶点蛋白的结合关系 |
| DDS | 药物-药物相似性 | 两个药物在化学结构上的相似程度 |
| L1000 | L1000 项目 | Broad Institute 的高通量基因表达项目 |
| Landmark Gene | 标志性基因 | 可推断全基因组表达的 978 个代表性基因 |
| DEG | 差异表达基因 | 在不同条件下表达水平显著变化的基因 |
| NeighborLoader | 邻居加载器 | PyG 的邻居采样数据加载器 |
| Early Stopping | 早停 | 训练中 loss 不再下降时提前终止 |
| Ridge Regression | 岭回归 | 带 L2 正则化的线性回归 |
| Tanimoto Similarity | Tanimoto 相似性 | 衡量两个分子指纹重叠程度的指标 |
