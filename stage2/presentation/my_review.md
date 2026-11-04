# 虚拟细胞药物扰动模型：数据基座、模型架构、评价与落地应用

## 摘要

药物效应在本质上表现为化学扰动所驱动的细胞状态改变；而化合物、剂量、时间、细胞情境与组合共同构成的高维实验空间，使穷举式筛选难以实现，药物扰动预测因而成为虚拟细胞研究的重要方向。本文沿“化学扰动→分子响应→机制/靶点→治疗决策”这一药理学链路，从数据基座、模型架构、评价与基准测试、落地应用四个方面，对虚拟细胞药物扰动模型加以系统梳理：在数据层面，归纳化学扰动转录组、药物敏感性与药理先验等资源及其固有偏倚；在架构层面，比较结构驱动、潜变量生成、最优传输、扩散与流匹配、图谱因果以及基础模型等主要技术路线；在评价层面，审视现有指标谱系与基线对照，并提出面向药理学的分层评价协议；在应用层面，按验证级别考察已进入体内与临床前的代表性案例。本文强调，此类模型的价值并不取决于表达量预测的精度，而在于能否产出机制对齐、剂量可信、跨情境可迁移且可经湿实验闭环验证的可行动证据。就现阶段而言，虚拟细胞药物扰动模型更接近药理学假设生成器，有望压缩药物发现早期的试错成本。

**关键词：** 虚拟细胞；药物扰动预测；数据基座；模型架构；基准测试；药物重定位

## 1 引言

自现代药理学诞生以来，为疾病寻找有效药物一直是该领域的主要目标之一。围绕此目标，研究者们发现并总结了大量药物靶点、信号通路与构效关系，并逐步建立起从靶点识别、先导化合物发现到临床前评价与临床转化的完整研究范式。然而，传统的药物发现路径高度依赖实验试错：剂量-反应关系、构效关系以及药代/药效学分析虽能刻画药物作用的宏观规律，却难以说明化学扰动究竟如何重塑细胞的分子状态，而单细胞化学基因组学则揭示出，同一化合物在不同细胞类型与状态下可以引发显著异质的响应[1]。与此同时，“化合物×剂量×时间×细胞情境×组合”所构成的高维实验空间，使穷举式筛选在实践中不可行，这正是药物扰动预测研究的根本动因。

在此背景下，虚拟细胞由数字孪生与世界模型式的长期愿景，逐步收敛为面向药物扰动的可操作预测器[2,3,4]。按输入与输出的差异，药物扰动预测大致可分为三类：药物→扰动组、转录组+药物→扰动组，以及转录组+扰动→扰动组，分别对应药物-效应映射、条件响应生成与扰动推断[5,6]。

这一方向之所以在近年来快速兴起，得益于数据与方法的同步成熟。L1000 连通图谱以 978 个标志基因公开了逾百万条化学扰动谱[7]，sci-Plex 在单细胞分辨率下平行获取了百余种化合物、约 65 万个细胞的转录组[8]，Tahoe-100M 等更大规模资源亦相继出现[9]；变分自编码器、最优传输、扩散模型、流匹配与基础模型，则为扰动响应的生成与表达提供了丰富的建模工具[5,6]。然而，多数方法仍停留于基准指标的比较，真正经剂量-机制-转化闭环完成湿实验验证者依然稀少[3]。

有鉴于此，本文提出“数据—架构—评价—应用”四支柱的分析框架，并以机制对齐、剂量可信、跨情境可迁移与湿实验闭环作为药理学层面的评判尺度：数据决定模型能力的上限，架构决定其可达的假设空间，评价决定其优化方向，应用则检验其真实价值。与既有综述相比，本文以药理学研究流程而非纯算法视角组织材料，单独设立批判性的评价章节，并提出面向药理学的分层评价协议与验证级别分层，强调可行动证据应优先于表达量预测精度。

## 2 数据基座

围绕“化学扰动→分子响应→机制/表型”这一逻辑链，药物扰动模型的数据基座大致可归纳为五类：化学扰动转录组（包括整体与单细胞两个层次）、细胞系组学与药物敏感性标签、形态学表型、参考图谱与对照/标准化资源，以及药理-生物学先验知识源。其中，前两类主要提供模型训练与监督信号，后三类则提供情境基准与机制约束。

### 2.1 化学扰动·整体转录组

L1000 平台以 978 个标志基因近似全转录组，约可恢复 82% 的信息；CMap（Connectivity Map，连通图谱）/LINCS（Library of Integrated Network-based Cellular Signatures，整合网络细胞特征库）在此基础上累积了 1,319,138 个表达谱，其中 19,811 个小分子以三重复、6 h 或 24 h 两个时间点进行处理，Touchstone（基准子集）部分覆盖九个核心细胞系，Discovery（扩展子集）部分则覆盖 3–77 个细胞系[7,10]。该资源已成为多数“结构→谱”模型的基座：DeepCE 将 L1000 中不可靠的实验转化为数据增广信号[11]，ExPO 直接预测 978 个基因的 z 分数谱，并对剂量-时间的连续暴露进行建模[12]。

### 2.2 化学扰动·单细胞转录组

sci-Plex 通过核哈希实现高通量单细胞化学转录组检测，在 A549、K562、MCF7 三个细胞系上处理了 188 个化合物（各设四个剂量），共获得约 650,000 个细胞与约 5,000 个独立样本，并据此揭示了 HDAC（histone deacetylase，组蛋白去乙酰化酶）抑制剂响应的细胞异质性[8]。规模更大的 Tahoe-100M（约 1 亿细胞、约 6 万条件、50 个细胞系）将数据粒度推进至十亿级规模[9]。scPerturb 则收集了 44 个公开单细胞扰动数据集（涵盖转录、蛋白与表观信息），统一质控与特征注释，并以 E 统计量量化扰动效应[13]。

### 2.3 细胞系组学与药物敏感性标签

敏感性读出把转录响应连接到“抑制/杀伤”这一药理学终点。GDSC（Genomics of Drug Sensitivity in Cancer，癌症药物敏感性基因组学）覆盖近 700 个细胞系、138 个抗癌药与约 75,000 次实验，给出了剂量-响应参数[14]；CCLE（Cancer Cell Line Encyclopedia，癌细胞系百科全书）将 947 个细胞系的组学图谱与其中 479 个细胞系对 24 个药物的 8 点剂量-响应曲线相配对[15]；PRISM（Profiling Relative Inhibition Simultaneously in Mixtures，混合物相对抑制同步分析）检测了 4,518 个药物在 578 个细胞系中的生长抑制效应[16]；CTRP（Cancer Therapeutics Response Portal，癌症治疗反应门户）则是另一个常用的整体药敏资源[17]。

### 2.4 形态学表型与参考图谱

JUMP（Joint Undertaking for Morphological Profiling，形态学分析联合项目）细胞绘画（Cell Painting）成像记录了约 136,000 个化学扰动的细胞形态特征，为转录组到形态的跨模态桥接提供了读出[18]。在参考图谱方面，Tabula Sapiens 提供了近 50 万个细胞（来自 24 个组织/器官、超过 400 个细胞类型）的单细胞转录组[19]，TCGA（The Cancer Genome Atlas，癌症基因组图谱）则提供了跨 12 类肿瘤的多组学分子图谱，支撑细胞情境与肿瘤情境的对齐[20]。

### 2.5 药理-生物学先验知识源

先验知识是缓解样本稀疏与情境外推失败的重要约束，通常通过图、本体或预训练的方式注入模型。按类型可大致分为以下几类：

蛋白互作/功能网络：STRING（Search Tool for the Retrieval of Interacting Genes/Proteins，互作基因/蛋白检索工具）被 GEARS、Stable-Shift、PDGrapher、TxPert、GenomeScalePert 用于基因间关系建模[21,22,23,24,25]；

通路与本体：GO（Gene Ontology，基因本体）、Reactome（Reactome 通路数据库）、KEGG（Kyoto Encyclopedia of Genes and Genomes，京都基因与基因组百科全书）、MSigDB（Molecular Signatures Database，分子特征数据库）将预测约束到通路结构，用于 GEARS、scLong、AttentionPert、VEGA 等[21,26,27,28]；

药物-靶点/机制库：DrugBank、ChEMBL、BindingDB、TTD（Therapeutic Target Database，治疗靶点数据库）、CTD（Comparative Toxicogenomics Database，比较毒理基因组学数据库）用于药物-靶关系表征，如 AMVL 综合 DrugBank、TTD、CTD 构建药物视图[29]，OffTargetIDL 借助互作组推断脱靶效应[30]，DeepCE 将 STRING 网络与 DrugBank 合成药物-靶向量[11]；

知识图谱：PrimeKG、PharMeBINet、MAP-KG 整合多源机制关系；MAP 构建的 MAP-KG 统一 14 个公共资源、连接约 18.7 万药物与 2.3 万基因[31]，AMVL 引入 PharMeBINet 的图嵌入[29]；

共表达/基因调控网络：SmallDataDrug、GPO-VAE、MiTCP、VEGA 用共表达或调控网络注入情境先验[28,32,33,34]；

蛋白序列/结构：UniProt、ESM-2（Evolutionary Scale Modeling 2，进化尺度建模模型 2）、AlphaFold 提供蛋白侧表示，见于 MAP 与 UNICORN[31,35]；

文本/本体语义：PubMed、MeSH（Medical Subject Headings，医学主题词）、Disease Ontology 被 AMVL、MAP 用于多模态语义融合[29,31]；

分子预训练：MoleculeSTM、ChemBERTa-2 分别被 MAP 与 ExPO 用作化学结构编码器[12,31]。

### 2.6 药理学视角的质量与偏差

从药理学角度看，当前数据基座仍存在若干值得警惕的质量与偏差问题。其一，剂量与时间覆盖不均，L1000 多固定于 6 h 或 24 h，且浓度梯度稀疏，模型因而容易学到细胞系或批次指纹，而非机制本身。其二，细胞情境偏倚明显，核心训练资源以癌细胞系为主，正常组织与广域剂量-效应覆盖不足，限制了外推能力。其三，溶剂（vehicle）/DMSO（dimethyl sulfoxide，二甲基亚砜）与阴性对照的设计不统一，影响差异谱与表型标签的可靠性。其四，单细胞向整体的伪化会掩盖细胞型异质性，而 MOA（mechanism of action，作用机制）注释的不一致也会带来标签噪声。其五，数据划分若未隔离同一细胞系或同一化合物，将造成信息泄漏，并高估模型跨药物、跨情境的泛化能力。

### 2.7 开放问题

当前数据基座尚有三个开放问题。一是规模与机制覆盖不匹配：L1000 覆盖广但机理解析弱，单细胞数据集机制细但规模窄、成本高。二是转录、形态、药敏三种模态之间的系统对齐不足。三是阴性/安慰剂对照与活性筛选之间的可迁移性尚缺闭环验证。

## 3 模型架构

### 3.1 化学结构驱动的直接预测

最早成熟的一类模型直接从药物化学结构或文本描述（可附加细胞系、剂量等条件）预测基因表达扰动谱，并不显式建模细胞对药物反应的异质性。该路线兴起于药物扰动数据相对稀缺的早期，其训练目标明确、不受单细胞异质性噪声的干扰，也因此在真实药物发现中获得湿实验验证的落地应用最多；其代价则是难以刻画个体或情境（细胞类型、时程）之间的差异。代表工作如 DLEPS，它将 SMILES（Simplified Molecular Input Line Entry System，简化分子线性输入规范）解析为语法树后映射到高维球面隐向量，再经稠密网络预测 L1000 扰动谱（皮尔逊相关系数约 0.74），并以 GSEA（gene set enrichment analysis，基因集富集分析）作为疾病疗效打分，据此在小鼠模型中验证了肥胖、高尿酸血症与非酒精性脂肪性肝炎的候选化合物[36]。DeepCE 以图神经网络与多头注意力同时建模化学子结构-基因与基因-基因关联，并用数据增强利用 L1000 的不可靠批次[11]；CIGER 以基因表达排序学习替代逐基因回归，更贴合扰动谱的排序性质[37]；ExPO 采用神经算子（DeepONet）将扰动谱表达为剂量-时间的连续函数，并结合 ChemBERTa 分子嵌入与符合校准的不确定度[12]。此外，TranSiGen、DeepICER、TensorTrain、CellBox、HerbTransfer、PRnet 等分别以自监督表示学习、双线性注意力、张量分解、显式动力学、迁移学习与深度生成，实现从结构或文本到扰动谱的预测[38,39,40,41,42,43]。

### 3.2 潜变量生成式模型

VAE（variational autoencoder，变分自编码器）家族以可解释的潜变量把“对照状态＋扰动条件”映射为扰动后的状态。scGen 利用潜空间向量算术实现跨细胞类型、跨物种的外推[44]；trVAE 引入最大均值差异正则，缓解非配对条件生成中的 OOD（out-of-distribution，分布外）问题[45]；CPA 与 chemCPA 将药物效应分解为可加的潜变量偏移，天然适配组合扰动与剂量外推，并通过化学编码支持对完全未见药物的预测[46,47]；SAMS-VAE 以稀疏加性机制偏移实现扰动专属潜子空间的解耦与组合推理[48]。该类模型的药理学适配核心在于“加性可解耦”：若扰动能够在潜空间中线性叠加，则组合用药与剂量响应即可由单药结果外推。其局限在于潜空间的解释性依赖先验设计，生成质量受限于训练扰动的覆盖范围。S-VQ-VAE、VEGA、BioInfoVAE、scPreGAN、scPerb、CRADLE-VAE、scNODE、Manatee、Dr.VAE、WAVE 等则分别以向量量化、基因模块接线、生物先验、对抗生成、风格迁移、反事实解耦、神经微分方程与时序建模进一步细化了该框架[28,49,50,51,52,53,54,55,56,57]。

### 3.3 最优传输与分布映射

单细胞测量具有破坏性，扰动前后的样本不成对，OT（optimal transport，最优传输）因而直接建模两个细胞群体分布之间的映射。CellOT 用输入凸神经网络求解 OT 映射，用于药物与蛋白扰动响应[58]；CINEMA-OT 在 OT 之上叠加因果框架，分离混淆因素并得到反事实细胞对，支持效应归因[59]；W1OT 改用 Wasserstein-1 距离以加快求解[60]；CMG 以条件 Monge Gap 使传输映射对药物、剂量、组合等协变量条件化，提升了对未见治疗的泛化[61]。其药理学价值在于用群体分布位移直接对应给药前后的机制性重排，但高维基因空间中的传输代价与可扩展性是主要约束；scPRAM 等则以注意力机制作为更轻量的替代[62]。

### 3.4 扩散与流匹配

扩散模型与条件流匹配通过迭代去噪学习高维条件分布，擅长多模态条件生成与连续剂量-时间场的建模。Squidiff 用扩散模型预测细胞发育与药物响应中的转录组变化[63]；DiffDebias（dbDiffusion）在 VAE 潜空间中以分类器无关引导实现去偏[64]；LDM-Drug 将 VAE 与潜扩散结合，预测条件感知的药物转录响应[65]；MCDP 强调多模态条件扩散以提升跨药物、跨细胞系的泛化[66]；CFM-GP 以条件流匹配统一跨细胞类型的扰动学习[67]；MorphoDiff 以转录组引导的潜扩散生成扰动下的细胞形态变化[68]。此类模型能够精细控制组合与连续暴露条件，代价是采样速度与训练成本较高。

### 3.5 图神经网络与知识图谱先验

把基因互作网络或知识图谱作为归纳偏置引入模型，可使预测受到生物通路的约束。GEARS 结合基因关系知识图预测多基因扰动的转录组[21]；TxPert 用多重知识图谱做潜变量迁移以预测扰动效应[24]；PDGrapher 求解逆问题，直接由疾病状态推荐可逆转表型的组合靶点[23]；OffTargetIDL 基于信号互作组推断药物的脱靶转录效应，并分离靶向（on-target）与脱靶（off-target）效应[30]；GPO-VAE、MiTCP、GenomeScalePert、AttentionPert 等亦以图先验或通路建模增强机制一致性[25,27,33,34]。这类模型的长处在于机制可解释性与靶点/脱靶归因，弱点则在于对图谱完整性与准确性的依赖。

### 3.6 解耦与因果建模

将机制效应与情境/批次效应显式分离，是保证扰动信号不受混淆、可被归因的关键。biolord 解耦已知与未知属性，并支持跨状态的“虚拟迁移”[69]；scDisInFact 显式区分批次效应与条件效应[70]；scCausalVI 分离细胞基础态与处理效应[71]；VCI 以变分因果推断估计干预效应[72]；CellCap 用注意力与稀疏字典学习把扰动响应分解为若干转录程序[73]；scGPA 借助大模型工作流进行方向性虚拟扰动[74]。从药理学角度看，这是从“预测表达”走向“识别可干预机制”的必要环节，但因果识别所需的假设在观测转录组中往往难以严格满足。

### 3.7 单细胞基础模型与表征学习

大规模预训练的单细胞基础模型（scGPT、scFoundation、Geneformer、UCE、CellFM、scLong、GeneCompass、scPRINT、BulkFormer 等）提供了通用的细胞/基因嵌入，以期零样本或少样本迁移到未表征的药物或新情境[26,75,76,77,78,79,80,81,82]；State 则专门以多样本集合表征预测跨情境扰动[83]。此类模型把“数据决定能力上限”推到极致，但它们在药物扰动任务上的实际增益高度依赖下游微调与数据泄漏控制，其有效性与成本之争留待评价章节详述。

### 3.8 多模态/形态学/空间/跨器官

从转录组这一核心读出向外扩展，相关工作也日益增多：IMPA（MorphResponse）以生成式风格迁移预测细胞形态扰动[84]；MorphoDiff 以转录组引导的扩散桥接形态表型[68]；TransTox（BridgingOrgan）用 GAN（generative adversarial network，生成对抗网络）进行肝-肾转录谱的双向翻译，以辅助多器官毒性评估[85]；Celcomen 引入空间因果解耦[86]；MultiPert 以对抗对齐与双注意力进行单细胞多组学扰动预测[87]；GenoME 以 MoE（mixture of experts，混合专家）生成模型整合多模态基因组谱[88]。其药理学意义在于，把虚拟筛选的输出从基因层面延伸至形态表型与跨器官毒性，使预测更贴近可操作的临床前证据。

### 3.9 架构评述与选型

若从机制可解释性、剂量/时间/组合、脱靶、跨情境 OOD 与不确定度等药理学需求综合考量，各类架构各有其适用边界：结构驱动模型适合扰动数据稀疏、需要快速虚拟筛选的场景，但缺乏情境与异质性信息；潜变量生成模型适合组合推理、剂量外推与机制解耦；OT 与扩散模型擅长高保真的分布建模，但计算与可解释性代价较高；图/知识图谱与因果建模提供最强的机制与靶点归因；基础模型在数据充分、需要零样本迁移时占优，而轻量专用模型在目标明确、可解释性优先时更为可控。据此，选型可归纳为：扰动覆盖有限且需优先实现湿实验转化者，宜选用结构驱动或轻量专用模型；需要组合/剂量/异质性推演者，宜选用可解耦的潜变量或条件扩散模型；需要靶点与脱靶归因者，宜选用图谱与因果框架；需要跨模态/跨器官毒性评估者，宜选用多模态配对或翻译模型。需要强调的是，无论选择何种架构，数据泄漏控制、OOD 与不确定度验证都是进入应用阶段前不可省略的环节。

## 4 模型评价与基准测试

### 4.1 评价的问题层级

评价应沿药理学证据链逐级上升：表达预测精度 → 通路/机制保真度 → 剂量、MOA、靶点与组合等可行动性 → 湿实验/临床前决策价值。基因层面的高相关并不自动等于机制可靠，因此须先明确评价目标，再选择与之匹配的指标与数据划分。

### 4.2 指标谱系

指标大致可分为六类。回归/相关类包括 MSE（mean squared error，均方误差）、R²、皮尔逊/斯皮尔曼相关，其核心在于评估扰动效应 Δ（处理组减对照），而非绝对表达；以对照为参考的皮尔逊 Δ 易被系统性变异放大，不依赖参考的 MSE 则更为稳健[89]。差异表达类以 DEG（differentially expressed genes，差异表达基因）重合度、AUPRC（area under the precision–recall curve，精确率-召回率曲线下面积）、Top-k 方向一致性衡量扰动基因的恢复[90,91]。通路/机制类通过 GO/Reactome、GSEA 与 MOA 分类检验生物学一致性[31,92]。剂量-反应类考察单调性与 EC50（half-maximal effective concentration，半数有效浓度）一致性。分布类如 MMD（maximum mean discrepancy，最大均值差异）/能量距离与 Wasserstein 距离，衡量整体位移；不确定度校准则检验输出的可信度。PerturBench 还提出排序指标，以真实扰动能否在候选中排到最前衡量扰动特异性，并借以暴露模式坍缩[91]。

### 4.3 基准测试全景

PerturBench 定义了协变量转移与组合预测任务，并引入排序指标[91]；BenchGenPert 覆盖 29 个数据集、6 类指标，区分扰动泛化与情境泛化，发现情境越不相似预测越差[93]；SysCompPert 用 25 个数据集、13 种方法、4 类指标，发现表达指标利于小效应扰动，而 Δ/DEG 指标在大效应下更为清晰[90]；Systema 揭示常规指标会被系统性变异高估[89]。scArchon 组合了统计与生物指标，Bench-FCM、Bench-LinBaseline、SimpleControls 则强调基线对照[94,95,96,97]。这些平台共同表明，不同指标探测的是不同的生物学信号，而非同一对象的冗余摘要[90]。

### 4.4 评价中的陷阱与争议

其中最一致的发现是，深度模型并未普遍击败简单基线。Anders 等报告深度预测尚未超越线性基线[96]；Moccia 等证明均值与 CRISPR（clustered regularly interspaced short palindromic repeats，成簇规律间隔短回文重复序列）先验对照即可揭示基础模型的有效性上限[97]；Bench-FCM 发现训练均值模型优于 scGPT/scFoundation，而带 GO 特征的随机森林优势显著[95]。其机制在于，扰动态与对照态之间存在系统性差异（选择偏差、应激反应等混杂），参考型指标会把只捕获平均效应的基线放大为“高能力”；Systema 在十个数据集中量化了该系统性变异，发现其与皮尔逊 Δ 的相关达 0.91–0.95，改用不依赖参考的 RMSE（root mean squared error，均方根误差）才能部分缓解[89]。即便定量分数较高，模型也可能丢失扰动特异特征，如 scArchon 所见[94]。更为根本的是，OOD 的定义缺乏统一标准：扰动、细胞类型、剂量与组合的划分难度迥异，使跨研究比较失真[93,98]；而数据泄漏、批次与对照设计缺陷则会进一步加剧偏差[5,6]。此外，报告规范亦普遍缺失。

### 4.5 面向药理学的分层评价协议

本文建议将评价组织为逐级门槛：分子表达 → 通路/机制 → MOA/靶点 → 剂量-反应 → 表型 → 体内，模型须逐级通过。强制报告应包含 Δ 指标、明确的 OOD 划分、至少一个线性/均值/对照基线、不确定度估计与复现脚本；在解读参考型指标之前，应先报告数据集自身的系统性变异程度[89]。最终验证与湿实验闭环对接，由计算评价与实验证实迭代构成“计算→实验”闭环[3]；只有经过这一检验的机制、剂量、组合与靶点结论，才可进入转化决策。

### 4.6 从模型预测推导出药物：下游方法论

基于配体/结构的方法借助指纹、相似性、药效团与对接进行筛选，但缺少细胞机制信息。基于表型/转录组的方法更贴合扰动模型：包括 CMap 连接性与疾病签名逆转、GSEA/ssGSEA（single-sample gene set enrichment analysis，单样本基因集富集分析）富集，其中 DLEPS 由疾病签名反推候选并经实验验证[36]；DeepCE 预测新化合物的差异表达谱以用于重定位[11]；Beyondcell 计算单细胞的药物特征富集分以刻画肿瘤异质性[99]；MAP 将预测聚合为通路评分，在 58 个留出化合物中把 5 个获批肺癌药中的 4 个排入前 15[31]。基于靶点/机制的方法则直指可解释性：OffTargetIDL 同时推断了药物-靶标作用及其下游信号[30]；PDGrapher 预测可逆转疾病表型的组合靶点[23]；MoA-SEN 以稀疏可解释网络预测机制类别并关联剂量-反应[92]；scRank 借靶标扰动网络从无处理数据推断响应细胞类型[100]。排序与决策应辅以多打分聚合、不确定性过滤、方向性与剂量依赖校验以及阳性/阴性对照，并警惕 L1000 仅测 978 个标志基因所带来的平台偏差[11,36]。预测只有被转化为可解释、可排序、可验证的“药物级证据”，方能进入药学的行动空间。

## 5 落地应用

虚拟细胞药物扰动模型的最终价值，在于把“化学扰动→分子响应”的预测能力转化为药理学上可行动的决策证据。典型的落地管线是：以疾病或基因签名定义治疗方向，用模型对候选化合物打分排序，随后完成机制解释，最终由湿实验闭环验证。该管线是否成立，取决于候选能否在体内模型中重现预测效应，因此“验证级别”是判断各案例证据强度的核心标尺[36]。

### 5.1 体内阳性案例

DLEPS 是这条管线中证据最为完整的代表：原始工作以 L1000 扰动谱训练，对未见分子表达谱变化的预测相关系数达 0.74，并在肥胖、高尿酸血症与非酒精性脂肪性肝炎三种小鼠模型中分别验证了紫苏烯、竹节参皂苷 IV 与曲美替尼的疗效，进而提示 MEK–ERK（mitogen-activated protein kinase kinase–extracellular signal-regulated kinase，丝裂原活化蛋白激酶激酶–细胞外信号调节激酶）通路可作为干预靶点[36]。此后，DLEPS 被系统用于多类疾病并取得动物模型阳性结果：坡模酸通过抑制 PTP1B（protein tyrosine phosphatase 1B，蛋白酪氨酸磷酸酶 1B）同时改善衰老与肥胖[101]；阿塔鲁伦、华蟾素与双氢青蒿素分别经 BMP（bone morphogenetic protein，骨形态发生蛋白）-SMAD、Wnt/β-catenin 及组蛋白 H3K9（histone H3 lysine 9，组蛋白 H3 第 9 位赖氨酸）乙酰化通路对抗骨质疏松[102,103,104]；升麻素苷（prim-O-glucosylcimifugin）改善衰老相关的肌腱再生[105]；新型巨核细胞生成诱导剂 DMAG 用于血小板减少[106]；基于白癜风签名的筛查则将氟维司群定位为黑色素瘤免疫治疗的增敏剂[107]。DeepCE 同样贡献了体内阳性案例，如异丁司特用于阿尔茨海默病抗炎[108]、喹唑啉类 LSD1（lysine-specific demethylase 1，赖氨酸特异性去甲基化酶 1）抑制剂用于肿瘤[109]。GPAR 预测并在 ApoE/LDLR（apolipoprotein E/low-density lipoprotein receptor，载脂蛋白 E/低密度脂蛋白受体）双基因缺陷小鼠中验证了复方丹参滴丸对高胆固醇性心衰的保护作用[110]；GPS 针对肝细胞癌与特发性肺纤维化完成了从头筛选与先导优化[111]；scGPT 辅助构建了损伤肾靶向肽库，并筛得可缓解急性肾损伤的靶向递送候选[75,112]。

### 5.2 临床前待核案例

另一类案例停留在体外或早期临床前阶段，证据层级较低：XPert 构建了 TOP1（topoisomerase I，拓扑异构酶 I）/PARP1（poly(ADP-ribose) polymerase 1，聚腺苷二磷酸核糖聚合酶 1）双靶点抑制剂用于膀胱癌[113]；TranSiGen 结合从头生成与专家结构优化，获得了纳摩尔级酪氨酸酶抑制剂[114]；scRank 从椎间盘退变单细胞数据推断 MDSC（myeloid-derived suppressor cells，髓系来源抑制细胞）与脂质运载蛋白-2（Lipocalin-2，LCN2）为治疗靶点[100,115]；MultiDCP 提出了 6 个具性别差异机制的阿尔茨海默病重定位候选[116]；DeepCE 扩展框架筛得可透血脑屏障的自噬增强剂[117]；DLEPS 与 DeepCE 结合 PROTAC（proteolysis-targeting chimera，蛋白降解靶向嵌合体）去卷积，将褪黑素定位为抗衰老与抗 AD（Alzheimer's disease，阿尔茨海默病）的双重候选[118]。

### 5.3 机制与靶点

药物响应预测同时服务于机制解析：OffTargetIDL 借助互作组深度学习推断药物的脱靶信号效应[30]；CellBox、PDGrapher 与 PerturbNet 分别以动力系统、因果图神经网络与深度生成模型预测组合扰动，以设计联合治疗[23,41,119]；CINEMA-OT 与 GPerturb 则提供单细胞层面的协同/拮抗分解与不确定性量化，使组合与剂量设计不再依赖黑箱打分[59,120]。

### 5.4 安全性与跨器官毒性

在安全性方面，BridgingOrgan/TransTox 用生成对抗网络实现肝-肾转录组的双向翻译，使单一器官的毒性数据可外推至多器官评估，契合药物警戒对早期毒性信号的需求[85]。

### 5.5 验证级别分层与落地障碍

综合来看，可将证据强度分为“计算→体外→体内→临床”四级；当前多数阳性结果仍停留在动物模型层面，进入临床的证据依然稀缺。落地障碍并非单纯的工程问题：模型与数据集的可复现性不足、监管对黑箱预测的接受度较低、数据隐私受限、可解释性有限、湿实验验证成本高昂，且“基因/细胞表型→临床终点”的映射仍不清晰。其核心结论是：预测分数高并不等于治疗获益，虚拟细胞模型必须经受从计算到临床的逐级验证，方能兑现其药理学价值。

## 6 讨论与展望

### 6.1 四支柱的相互制约

数据、架构、评价与应用并非彼此并列的独立环节，而是同一闭环中的四个环节。数据决定模型能力的上限：单细胞化学基因组学的兴起，使药物效应得以在细胞亚群与机制层面被刻画[1]，而大规模时间分辨的扰动蛋白组数据进一步表明，更接近功能的读出能够跨细胞系、患者来源类器官与临床活检迁移，从而支撑可操作的疗效与协同预测[121]。架构决定模型可达的假设空间，评价则决定其优化方向；若评价止步于少量指标，模型便容易陷入预测指标看似良好、但药理学上难以行动的困境。闭环验证要求计算评估与湿实验复证交替迭代，从而推动模型从方法新颖性走向机制可行动性[3]。

### 6.2 关键开放问题

其一，是真正未见药物与靶点的 OOD（out-of-distribution，分布外）外推，这也是扰动建模的首要目标[6]；虚拟细胞世界模型主张，其核心优势并不在于基准指标上的增量，而在于输出本质的转变，即支持反事实推演与长时程规划[2]。其二，是组合与剂量-时间外推、个体化与患者来源建模，以及体外至体内的迁移。后者尤为现实：在监管层面，人工智能模型目前仅作为补充证据，独立的新药申请仍须经湿实验验证[3]。此外，预测还须区分关联与因果，并对不确定性做出定量刻画[4,5]。

### 6.3 未来路线图

未来研究应走向多尺度、多模态的虚拟细胞世界模型，以干预条件化的模拟取代孤立的端点预测[2,4]；构建“预测—实验—再训练”的主动闭环，以缩小假设空间[4]；将基础模型与生物先验相结合；并建设 scPerturb 式的扰动图谱与标准化基准，以弥补当前各自比较、缺乏最佳实践的短板[6]。与此同时，监管与伦理框架也须同步推进[3]。

### 6.4 对药理学研究者的实践建议

研究团队宜按“任务—数据—架构”相匹配的原则选题，评价至少应覆盖 OOD、剂量-时间与分布级一致性检验，并与湿实验团队建立固定的闭环协作节律[3,6]。在定位上则须审慎：虚拟细胞是药理学假设生成器，而非湿实验的替代品[4]。

## 7 结论

虚拟细胞药物扰动模型的四个支柱正处于发展不均衡却又相互咬合的阶段：扰动数据与模型架构的积累最为迅速，评价环节则仍是最薄弱、也最具决定性的短板，而真实应用层面的验证目前仍集中在少数起步较早的结构驱动模型。就当前定位而言，这类模型更接近药理学假设生成器，其价值并不在于预测了多少基因的表达，而在于能否产出机制对齐、剂量可信、跨情境可迁移且可被湿实验闭环验证的可行动证据[3]。随着分层评价协议与标准化基准的推广，大规模高质量扰动数据对模型训练的持续支撑[1]，以及虚拟细胞“计算—实验—转化”闭环迭代的逐步普及[2,3]，虚拟细胞药物扰动模型有望实质性压缩药物发现早期的试错成本，并从辅助性工具走向临床前药物评价的核心组件[3]。

## 参考文献

[1] Wyatt A, Hoffer-Hawlik K, Giglio RM, et al. An expanded role for single-cell chemical genomics profiling in drug discovery. Biochemical Journal, 2026. DOI: 10.1042/bcj20253273.

[2] Xing EP, Song L. A world model of the virtual cell. Cell, 2026. DOI: 10.1016/j.cell.2026.08.042.

[3] Ma C, Zhang H, Rao Y, et al. AI-driven virtual cell models in preclinical research: technical pathways, validation mechanisms, and clinical translation potential. npj Digital Medicine, 2025. DOI: 10.1038/s41746-025-02198-6.

[4] Bunne C, Roohani Y, Rosen Y, et al. How to build the virtual cell with artificial intelligence: Priorities and opportunities. Cell, 2024. DOI: 10.1016/j.cell.2024.11.015.

[5] Monfort-Lanzas P, Rungger K, Madersbacher L, et al. Machine learning to dissect perturbations in complex cellular systems. Computational and Structural Biotechnology Journal, 2025. DOI: 10.1016/j.csbj.2025.02.028.

[6] Gavriilidis GI, Vasileiou V, Orfanou A, et al. A mini-review on perturbation modelling across single-cell omic modalities. Computational and Structural Biotechnology Journal, 2024. DOI: 10.1016/j.csbj.2024.04.058.

[7] Subramanian A, Narayan R, Corsello SM, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell, 2017. DOI: 10.1016/j.cell.2017.10.049.

[8] Srivatsan SR, McFaline-Figueroa JL, Ramani V, et al. Massively multiplex chemical transcriptomics at single-cell resolution. Science, 2020. DOI: 10.1126/science.aax6234.

[9] Zhang J, Ubas AA, de Borja R, et al. Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling. bioRxiv, 2025. DOI: 10.1101/2025.02.20.639398.

[10] Koleti A, Terryn R, Stathias V, et al. Data Portal for the Library of Integrated Network-based Cellular Signatures (LINCS) program: integrated access to diverse large-scale cellular perturbation response data. Nucleic Acids Research, 2018. DOI: 10.1093/nar/gkx1063.

[11] Pham T, Qiu Y, Zeng J, et al. A deep learning framework for high-throughput mechanism-driven phenotype compound screening and its application to COVID-19 drug repurposing. Nature Machine Intelligence, 2021. DOI: 10.1038/s42256-020-00285-9.

[12] Spadaro A, Sharma A, Dehzangi I. ExPO: an exposure-conditioned neural operator for L1000 signature prediction. Journal of Cheminformatics, 2026. DOI: 10.1186/s13321-026-01226-1.

[13] Peidli S, Green TD, Shen C, et al. scPerturb: harmonized single-cell perturbation data. Nature Methods, 2024. DOI: 10.1038/s41592-023-02144-y.

[14] Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Research, 2012. DOI: 10.1093/nar/gks1111.

[15] Barretina J, Caponigro G, Stransky N, et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature, 2012. DOI: 10.1038/nature11003.

[16] Corsello SM, Nagari RT, Spangler RD, et al. Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. Nature Cancer, 2020. DOI: 10.1038/s43018-019-0018-6.

[17] Cancer Therapeutics Response Portal (CTRP) [DB/OL]. Broad Institute. URL: https://portals.broadinstitute.org/ctrp/.

[18] Chandrasekaran SN, Ackerman J, Alix E, et al. JUMP Cell Painting dataset: morphological impact of 136,000 chemical and genetic perturbations. bioRxiv, 2023. DOI: 10.1101/2023.03.23.534023.

[19] The Tabula Sapiens Consortium*, Jones RC, Karkanias J, et al. The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. Science, 2022. DOI: 10.1126/science.abl4896.

[20] The Cancer Genome Atlas Research Network, Weinstein JN, Collisson EA, et al. The Cancer Genome Atlas Pan-Cancer analysis project. Nature Genetics, 2013. DOI: 10.1038/ng.2764.

[21] Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. Nature Biotechnology, 2024. DOI: 10.1038/s41587-023-01905-6.

[22] Dip SA, Zhang L. Stable-Shift: Predicting Transcriptional Responses of Unseen Gene Perturbations Using Graph Neural Networks with Biological Priors. Proceedings of the 17th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics, 2026. DOI: 10.1145/3807503.3820871.

[23] Gonzalez G, Lin X, Herath I, et al. Combinatorial prediction of therapeutic perturbations using causally inspired neural networks. Nature Biomedical Engineering, 2025. DOI: 10.1038/s41551-025-01481-x.

[24] Wenkel F, Tu W, Masschelein C, et al. TxPert: using multiple knowledge graphs for prediction of transcriptomic perturbation effects. Nature Biotechnology, 2026. DOI: 10.1038/s41587-026-03113-4.

[25] Zhan L, Wang Y, Wang A, et al. A genome-scale deep learning model to predict gene expression changes of genetic perturbations from multiplex biological networks. Briefings in Bioinformatics, 2024. DOI: 10.1093/bib/bbae433.

[26] Bai D, Mo S, Zhang R, et al. scLong: a billion-parameter foundation model for capturing long-range gene context in single-cell transcriptomics. Nature Communications, 2026. DOI: 10.1038/s41467-026-69102-y.

[27] Bai D, Ellington CN, Mo S, et al. AttentionPert: accurately modeling multiplexed genetic perturbations with multi-scale effects. Bioinformatics, 2024. DOI: 10.1093/bioinformatics/btae244.

[28] Seninge L, Anastopoulos I, Ding H, et al. VEGA is an interpretable generative model for inferring biological network activity in single-cell transcriptomics. Nature Communications, 2021. DOI: 10.1038/s41467-021-26017-0.

[29] Yan Y, Yang Y, Tong Z, et al. Adaptive multi-view learning method for enhanced drug repurposing using chemical-induced transcriptional profiles, knowledge graphs, and large language models. Journal of Pharmaceutical Analysis, 2025. DOI: 10.1016/j.jpha.2025.101275.

[30] Meimetis N, Lauffenburger DA, Nilsson A. Inference of drug off-target effects on cellular signaling using interactome-based deep learning. iScience, 2024. DOI: 10.1016/j.isci.2024.109509.

[31] Feng J, Zhao Z, Zhang X, et al. A knowledge-driven framework for predicting single-cell responses for unprofiled drugs. Nature Machine Intelligence, 2026. DOI: 10.1038/s42256-026-01286-w.

[32] Alsulami R, Lehmann R, Khan SA, et al. Predicting and interpreting cell-type-specific drug responses in the small-data regime using inductive priors. Nature Machine Intelligence, 2026. DOI: 10.1038/s42256-026-01202-2.

[33] Baek S, Park S, Chok YT, et al. GPO-VAE: modeling explainable gene perturbation responses utilizing GRN-aligned parameter optimization. Bioinformatics, 2025. DOI: 10.1093/bioinformatics/btaf256.

[34] Yang K, Cheng J, Cao S, et al. Predicting transcriptional changes induced by molecules with MiTCP. Briefings in Bioinformatics, 2024. DOI: 10.1093/bib/bbaf006.

[35] Liu T, Huang T, Wang L, et al. UNICORN: Towards universal cellular expression prediction with a multi-task learning framework. Nature Communications, 2025. DOI: 10.1038/s41467-025-64506-8.

[36] Zhu J, Wang J, Wang X, et al. Prediction of drug efficacy from transcriptional profiles with deep learning. Nature Biotechnology, 2021. DOI: 10.1038/s41587-021-00946-z.

[37] Pham T, Qiu Y, Liu J, et al. Chemical-induced gene expression ranking and its application to pancreatic cancer drug repurposing. Patterns, 2022. DOI: 10.1016/j.patter.2022.100441.

[38] Tong X, Qu N, Kong X, et al. Deep representation learning of chemical-induced transcriptional profile for phenotype-based drug discovery. Nature Communications, 2024. DOI: 10.1038/s41467-024-49620-3.

[39] Meng F, Wang C, Lin Y, et al. DeepICER: A deep learning framework for predicting compound-induced gene expression profiles. Acta Pharmaceutica Sinica B, 2026. DOI: 10.1016/j.apsb.2026.01.046.

[40] Iwata M, Yuan L, Zhao Q, et al. Predicting drug-induced transcriptome responses of a wide range of human cell lines by a novel tensor-train decomposition algorithm. Bioinformatics, 2019. DOI: 10.1093/bioinformatics/btz313.

[41] Yuan B, Shen C, Luna A, et al. CellBox: Interpretable Machine Learning for Perturbation Biology with Application to the Design of Cancer Combination Therapy. Cell Systems, 2021. DOI: 10.1016/j.cels.2020.11.013.

[42] Liu Q, Wang B, Li S. From chemical compounds to herbal interventions: a transfer learning-based perturbational transcriptome prediction framework for cancer drug discovery. Cancer Biology & Medicine, 2026. DOI: 10.20892/j.issn.2095-3941.2026.0032.

[43] Qi X, Zhao L, Tian C, et al. Predicting transcriptional responses to novel chemical perturbations using deep generative model for drug discovery. Nature Communications, 2024. DOI: 10.1038/s41467-024-53457-1.

[44] Lotfollahi M, Wolf FA, Theis FJ. scGen predicts single-cell perturbation responses. Nature Methods, 2019. DOI: 10.1038/s41592-019-0494-8.

[45] Lotfollahi M, Naghipourfar M, Theis FJ, et al. Conditional out-of-distribution generation for unpaired data using transfer VAE. Bioinformatics, 2020. DOI: 10.1093/bioinformatics/btaa800.

[46] Lotfollahi M, Klimovskaia Susmelj A, De Donno C, et al. Predicting cellular responses to complex perturbations in high-throughput screens. Molecular Systems Biology, 2023. DOI: 10.15252/msb.202211517.

[47] Hetzel L, Boehm S, Kilbertus N, et al. Predicting Cellular Responses to Novel Drug Perturbations at a Single-Cell Resolution. Advances in Neural Information Processing Systems 35, 2022. DOI: 10.52202/068431-1937.

[48] Bereket M, Karaletsos T. Modelling Cellular Perturbations with the Sparse Additive Mechanism Shift Variational Autoencoder. Advances in Neural Information Processing Systems 36, 2023. DOI: 10.52202/075280-0001.

[49] Xue Y, Ding MQ, Lu X. Learning to encode cellular responses to systematic perturbations with deep generative models. npj Systems Biology and Applications, 2020. DOI: 10.1038/s41540-020-00158-2.

[50] Doncevic D, Herrmann C. Biologically informed variational autoencoders allow predictive modeling of genetic and drug-induced perturbations. Bioinformatics, 2023. DOI: 10.1093/bioinformatics/btad387.

[51] Wei X, Dong J, Wang F. scPreGAN, a deep generative model for predicting the response of single-cell expression to perturbation. Bioinformatics, 2022. DOI: 10.1093/bioinformatics/btac357.

[52] Tang Z, Zhou M, Zhang K, et al. scPerb: Predict single-cell perturbation via style transfer-based variational autoencoder. Journal of Advanced Research, 2025. DOI: 10.1016/j.jare.2024.10.035.

[53] Baek S, Park S, Chok YT, et al. Cradle-VAE: Enhancing Single-Cell Gene Perturbation Modeling with Counterfactual Reasoning-based Artifact Disentanglement. Proceedings of the AAAI Conference on Artificial Intelligence, 2025. DOI: 10.1609/aaai.v39i15.33695.

[54] Zhang J, Larschan E, Bigness J, et al. scNODE: generative model for temporal single cell transcriptomic data prediction. Bioinformatics, 2024. DOI: 10.1093/bioinformatics/btae393.

[55] Yang Y, Seninge L, Wang Z, et al. The manatee variational autoencoder model for predicting gene expression alterations caused by transcription factor perturbations. Scientific Reports, 2024. DOI: 10.1038/s41598-024-62620-z.

[56] Rampášek L, Hidru D, Smirnov P, et al. Dr.VAE: improving drug response prediction via modeling of drug perturbation effects. Bioinformatics, 2019. DOI: 10.1093/bioinformatics/btz158.

[57] Lv T, Chen B, Dai X, et al. Learning chemical-induced gene expression perturbations with WAVE. Nature Communications, 2026. DOI: 10.1038/s41467-026-77645-3.

[58] Bunne C, Stark SG, Gut G, et al. Learning single-cell perturbation responses using neural optimal transport. Nature Methods, 2023. DOI: 10.1038/s41592-023-01969-x.

[59] Dong M, Wang B, Wei J, et al. Causal identification of single-cell experimental perturbation effects with CINEMA-OT. Nature Methods, 2023. DOI: 10.1038/s41592-023-02040-5.

[60] Chen Y, Hu Z, Chen W, et al. Fast and scalable Wasserstein-1 neural optimal transport solver for single-cell perturbation prediction. Bioinformatics, 2025. DOI: 10.1093/bioinformatics/btaf253.

[61] Driessen A, Rajwade DA, Harsanyi B, et al. Conditional Monge Gap enables generalizable single-cell perturbation modelling. Nature Machine Intelligence, 2026. DOI: 10.1038/s42256-026-01242-8.

[62] Jiang Q, Chen S, Chen X, et al. scPRAM accurately predicts single-cell gene expression perturbation response based on attention mechanism. Bioinformatics, 2024. DOI: 10.1093/bioinformatics/btae265.

[63] He S, Zhu Y, Tavakol DN, et al. Squidiff: predicting cellular development and responses to perturbations using a diffusion model. Nature Methods, 2026. DOI: 10.1038/s41592-025-02877-y.

[64] Shang E, Wei Y, Roeder K. Predicting the unseen: A diffusion-based debiasing framework for transcriptional response prediction at single-cell resolution. Proceedings of the National Academy of Sciences, 2025. DOI: 10.1073/pnas.2525268122.

[65] Kim C, Yoo S. Predicting condition-aware drug-induced transcriptional responses via a latent diffusion model. Bioinformatics, 2026. DOI: 10.1093/bioinformatics/btag173.

[66] Ma Y, Du K, Li Y, et al. Enhancing cross-context generalization in drug perturbation prediction with a multimodal conditional diffusion framework. Bioinformatics, 2026. DOI: 10.1093/bioinformatics/btag482.

[67] Abir AR, Dip SA, Zhang L. CFM-GP: unified conditional flow matching to learn gene perturbation across cell types. NAR Genomics and Bioinformatics, 2026. DOI: 10.1093/nargab/lqag087.

[68] Wang X, Fan Y, Guo Y, et al. Prediction of cellular morphology changes under perturbations with a transcriptome-guided diffusion model. Nature Communications, 2025. DOI: 10.1038/s41467-025-63478-z.

[69] Piran Z, Cohen N, Hoshen Y, et al. Disentanglement of single-cell data with biolord. Nature Biotechnology, 2024. DOI: 10.1038/s41587-023-02079-x.

[70] Zhang Z, Zhao X, Bindra M, et al. scDisInFact: disentangled learning for integration and prediction of multi-batch multi-condition single-cell RNA-sequencing data. Nature Communications, 2024. DOI: 10.1038/s41467-024-45227-w.

[71] An S, Cho J, Cao K, et al. scCausalVI disentangles single-cell perturbation responses with causality-aware generative model. Cell Systems, 2025. DOI: 10.1016/j.cels.2025.101443.

[72] Liu E, Zhang J, Uhler C. Learning genetic perturbation effects with variational causal inference. PLOS Computational Biology, 2026. DOI: 10.1371/journal.pcbi.1013194.

[73] Xu Y, Fleming S, Tegtmeyer M, et al. Explainable modeling of single-cell perturbation data using attention and sparse dictionary learning. Cell Systems, 2025. DOI: 10.1016/j.cels.2025.101245.

[74] Tang H, Li H, Zhao Q, et al. scGPA: an LLM-assisted workflow for directional virtual gene perturbation analysis from single-cell transcriptomes. BMC Genomics, 2026. DOI: 10.1186/s12864-026-13270-0.

[75] Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods, 2024. DOI: 10.1038/s41592-024-02201-0.

[76] Hao M, Gong J, Zeng X, et al. Large-scale foundation model on single-cell transcriptomics. Nature Methods, 2024. DOI: 10.1038/s41592-024-02305-7.

[77] Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. Nature, 2023. DOI: 10.1038/s41586-023-06139-9.

[78] Rosen Y, Roohani Y, Agrawal A, et al. Universal cell embedding provides a foundation model for cell biology. Nature, 2026. DOI: 10.1038/s41586-026-10689-z.

[79] Zeng Y, Xie J, Shangguan N, et al. CellFM: a large-scale foundation model pre-trained on transcriptomics of 100 million human cells. Nature Communications, 2025. DOI: 10.1038/s41467-025-59926-5.

[80] Yang X, Liu G, Feng G, et al. GeneCompass: deciphering universal gene regulatory mechanisms with a knowledge-informed cross-species foundation model. Cell Research, 2024. DOI: 10.1038/s41422-024-01034-y.

[81] Kalfon J, Samaran J, Peyré G, et al. scPRINT: pre-training on 50 million cells allows robust gene network predictions. Nature Communications, 2025. DOI: 10.1038/s41467-025-58699-1.

[82] Kang B, Fan R, Yi M, et al. BulkFormer: A large-scale foundation model for bulk transcriptomes. Cell Systems, 2026. DOI: 10.1016/j.cels.2026.101657.

[83] Adduri AK, Gautam D, Bevilacqua B, et al. Predicting cellular responses to perturbation across diverse contexts with State. Cell, 2026. DOI: 10.1016/j.cell.2026.07.052.

[84] Palma A, Theis FJ, Lotfollahi M. Predicting cell morphological responses to perturbations using generative modeling. Nature Communications, 2025. DOI: 10.1038/s41467-024-55707-8.

[85] Li T, Chen X, Tong W. Bridging organ transcriptomics for advancing multiple organ toxicity assessment with a generative AI approach. npj Digital Medicine, 2024. DOI: 10.1038/s41746-024-01317-z.

[86] Megas S, Chen DG, Polanski K, et al. Celcomen: spatial causal disentanglement for single-cell and tissue perturbation modeling. Nature Communications, 2026. DOI: 10.1038/s41467-026-69856-5.

[87] Zhao M, Tang X, Li J, et al. MultiPert: An adversarial alignment and dual attention framework for single-cell multi-omics perturbation prediction. PLOS Computational Biology, 2026. DOI: 10.1371/journal.pcbi.1014054.

[88] Wei J, Xue Y, Chai H, et al. GenoME: a MoE-based generative model for individualized, multimodal prediction and perturbation of genomic profiles. Nucleic Acids Research, 2026. DOI: 10.1093/nar/gkag902.

[89] Viñas Torné R, Wiatrak M, Piran Z, et al. Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. Nature Biotechnology, 2026. DOI: 10.1038/s41587-025-02777-8.

[90] Li L, You Y, Fu Y, et al. A systematic comparison of single-cell perturbation response prediction models. Science Advances, 2026. DOI: 10.1126/sciadv.aed3414.

[91] Wu Y, Wershof E, Schmon S, et al. PerturBench: Benchmarking Machine Learning Models for Cellular Perturbation Analysis. Advances in Neural Information Processing Systems 38, 2025. DOI: 10.52202/085713-3225.

[92] Sada Del Real K, Rubio A. Discovering the mechanism of action of drugs with a sparse explainable network. eBioMedicine, 2023. DOI: 10.1016/j.ebiom.2023.104767.

[93] Wei Z, Wang Y, Gao Y, et al. Benchmarking algorithms for generalizable single-cell perturbation response prediction. Nature Methods, 2026. DOI: 10.1038/s41592-025-02980-0.

[94] Radig J, Droit R, Doncevic D, et al. scArchon: a scalable benchmarking framework for assessing single-cell perturbation models. Genome Biology, 2026. DOI: 10.1186/s13059-026-04104-z.

[95] Csendes G, Sanz G, Szalay KZ, et al. Benchmarking foundation cell models for post-perturbation RNA-seq prediction. BMC Genomics, 2025. DOI: 10.1186/s12864-025-11600-2.

[96] Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. Nature Methods, 2025. DOI: 10.1038/s41592-025-02772-6.

[97] Wong DR, Hill AS, Moccia R. Simple controls exceed best deep learning algorithms and reveal foundation model effectiveness for predicting genetic perturbations. Bioinformatics, 2025. DOI: 10.1093/bioinformatics/btaf317.

[98] Kernfeld E, Yang Y, Weinstock JS, et al. A comparison of computational methods for expression forecasting. Genome Biology, 2025. DOI: 10.1186/s13059-025-03840-y.

[99] Fustero-Torre C, Jiménez-Santos MJ, García-Martín S, et al. Beyondcell: targeting cancer therapeutic heterogeneity in single-cell RNA-seq data. Genome Medicine, 2021. DOI: 10.1186/s13073-021-01001-x.

[100] Li C, Shao X, Zhang S, et al. scRank infers drug-responsive cell types from untreated scRNA-seq data using a target-perturbed gene regulatory network. Cell Reports Medicine, 2024. DOI: 10.1016/j.xcrm.2024.101568.

[101] Zhu J, An Y, Wang X, et al. The natural product rotundic acid treats both aging and obesity by inhibiting PTP1B. Life Medicine, 2022. DOI: 10.1093/lifemedi/lnac044.

[102] Zeng L, Gu R, Li W, et al. Ataluren prevented bone loss induced by ovariectomy and aging in mice through the BMP-SMAD signaling pathway. Biomedicine & Pharmacotherapy, 2023. DOI: 10.1016/j.biopha.2023.115332.

[103] Lu D, Zeng L, Li Y, et al. Cinobufotalin prevents bone loss induced by ovariectomy in mice through the BMPs/SMAD and Wnt/β-catenin signaling pathways. Animal Models and Experimental Medicine, 2024. DOI: 10.1002/ame2.12359.

[104] Wang R, Wang Y, Niu Y, et al. Deep Learning-Predicted Dihydroartemisinin Rescues Osteoporosis by Maintaining Mesenchymal Stem Cell Stemness through Activating Histone 3 Lys 9 Acetylation. ACS Central Science, 2023. DOI: 10.1021/acscentsci.3c00794.

[105] Wang Y, Jin S, Luo D, et al. Prim-O-glucosylcimifugin ameliorates aging-impaired endogenous tendon regeneration by rejuvenating senescent tendon stem/progenitor cells. Bone Research, 2023. DOI: 10.1038/s41413-023-00288-3.

[106] Wang L, Liu S, Luo J, et al. Targeting a thrombopoietin-independent strategy in the discovery of a novel inducer of megakaryocytopoiesis, DMAG, for the treatment of thrombocytopenia. Haematologica, 2022. DOI: 10.3324/haematol.2022.282209.

[107] Zhu J, Huang L, Bi S, et al. Vitiligo Signature-Based Drug Screening Identifies Fulvestrant as a Novel Immunotherapy Combination Strategy. Advanced Science, 2025. DOI: 10.1002/advs.202503979.

[108] Oliveros G, Wallace CH, Chaudry O, et al. Repurposing ibudilast to mitigate Alzheimer’s disease by targeting inflammation. Brain, 2023. DOI: 10.1093/brain/awac136.

[109] Li Z, Qin T, Li Z, et al. Discovery of quinazoline derivatives as a novel class of potent and in vivo efficacious LSD1 inhibitors by drug repurposing. European Journal of Medicinal Chemistry, 2021. DOI: 10.1016/j.ejmech.2021.113778.

[110] Yang Y, Feng K, Yuan L, et al. Compound Danshen Dripping Pill inhibits hypercholesterolemia/atherosclerosis-induced heart failure in ApoE and LDLR dual deficient mice via multiple mechanisms. Acta Pharmaceutica Sinica B, 2023. DOI: 10.1016/j.apsb.2022.11.012.

[111] Xing J, Tan M, Leshchiner D, et al. Deep-learning-based de novo discovery and design of therapeutics that reverse disease-associated transcriptional phenotypes. Cell, 2026. DOI: 10.1016/j.cell.2026.02.016.

[112] Xiao Y, Tong Z, Xu H, et al. A rationally designed injury kidney targeting peptide library and its application in rescuing acute kidney injury. Science Advances, 2025. DOI: 10.1126/sciadv.adt3943.

[113] Wang W, Li D, Liu J, et al. Molecular Construction and Anti-Bladder Cancer Study of Novel TOP1/PARP1 Dual-Target Inhibitors. Journal of Medicinal Chemistry, 2026. DOI: 10.1021/acs.jmedchem.6c01714.

[114] Sun Y, Wang J, Chen W, et al. Strategy and efficiency-redefined discovery of novel nanomolar tyrosinase inhibitors: AI de novo molecular generation + expert-guided structural optimization. Journal of Advanced Research, 2026. DOI: 10.1016/j.jare.2025.12.041.

[115] Zhang C, Li H, Wang H, et al. Identifying Myeloid-Derived Suppressor Cells and Lipocalin-2 as Therapeutic Targets for Intervertebral Disc Degeneration. Advanced Science, 2025. DOI: 10.1002/advs.202500505.

[116] Figueiredo-Pereira ME, Serrano PA, Rockwell P. Six promising drug repurposing candidates for Alzheimer’s disease and their sex-specific mechanisms and efficacy. Neural Regeneration Research, 2026. DOI: 10.4103/nrr.nrr-d-25-00256.

[117] Dong Y, Xiao X, Zhuang X, et al. DeepDrugDiscovery identifies blood–brain barrier permeable autophagy enhancers for Alzheimer’s disease. Nature Biomedical Engineering, 2026. DOI: 10.1038/s41551-026-01667-x.

[118] Sun Y, Liu S, Chen L, et al. AI-driven discovery of dual antiaging and anti-AD therapeutics via PROTAC target deconvolution of a super-enhancer–regulated axis. Science Advances, 2025. DOI: 10.1126/sciadv.adz9283.

[119] Yu H, Qian W, Song Y, et al. PerturbNet predicts single-cell responses to unseen chemical and genetic perturbations. Molecular Systems Biology, 2025. DOI: 10.1038/s44320-025-00131-3.

[120] Xing H, Yau C. GPerturb: Gaussian process modelling of single-cell perturbation data. Nature Communications, 2025. DOI: 10.1038/s41467-025-61165-7.

[121] Sun R, Qian L, Li Y, et al. An operational perturbation proteomics-based virtual cell model. Nature, 2026. DOI: 10.1038/s41586-026-11001-9.
