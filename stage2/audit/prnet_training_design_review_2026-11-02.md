# PRnet 正式预训练与心肌细胞微调：事实核查与决策访谈

日期：2026-11-02。状态：访谈已收拢，用户已批准按统一spec拆票，11张本地票据已发布；资格审计、QC与正式训练尚未执行。

## 已确认的研究目标

用户在本轮访谈确认：

1. 首先评估对未见药物的心肌细胞转录响应预测，并检验 Tahoe 预训练是否改善预测。
2. 模型面向人源心肌细胞的一般药物响应，后续用于 MIRI 筛选。

3. 用户选择Q3 B：目标药物在CM训练中未见，允许Tahoe预训练包含该药物；不要求全局cold-drug。实际应用面向尚无CM药物湿实验响应的药物，评估通过隐藏已有CM响应来模拟。
4. 用户选择Q4 A：主评估聚焦审定后的药物；环境化合物是否用于扩展训练另行决定。
5. 用户选择Q5 A：主要交付背景为iPSC-CM。
6. 用户明确收紧Q6：本轮仅24h；其他时间点不训练也不测试，不作为辅助训练或外部测试。
7. 用户选择Q7 A：temp_env药物面板为主实验，ScreenSeq作独立来源验证；合并微调属于后续消融，不是主训练。
8. 用户选择Q8 A：主训练只用审定后的药物；其他环境化合物仅进入单独的扩展训练消融，保持药物测试对象固定。
9. 用户选择Q9 B：主实验纳入人用已批准药物、临床在研药物，以及机制明确的研究工具化合物。用途与机制须有证据，ATC/MoA非空不自动获得资格；不限定FDA批准药物。
10. 用户选择Q10 A：通过实验质量审核的低响应、毒性/应激强响应均保留，增加标签并分层评价；不按DEG数量或响应强弱选择训练样本。实验质量不合格仍可排除，理由须可追溯。
11. 用户按推荐确认Q11：首轮保留18,164基因；新增160另作疾病相关覆盖审查，不更改冻结corpus。
12. Q12的原选择已被用户Q20修订替代：方向AUPRC@200为主指标，Δ RMSE、Δ相关性与响应分层结果为辅助；不是原460化合物历史leaderboard的直接延续。
13. 用户按推荐确认Q13：先用现有hidden128/state latent64结构完成random-init/Tahoe-init比较，随后单独消融潜空间与隐藏宽度；每种结构都要配套两种初始化。
14. 用户选择Q14 B：新药推断只使用训练数据形成的固定CM共识对照，不要求目标培养背景新测对照；不得使用测试药物处理表达形成该输入。共识的训练/推断使用方式及目标参数化需进一步决定。
15. 用户修订Q15：本轮仅Tahoe＋temp_env，ScreenSeq不训练、不验证、不测试。此前ScreenSeq独立来源验证及合并微调消融均从当前执行范围撤出，留作后续候选；五折等其余评估条款不因移除ScreenSeq而自动改变。
16. 用户选择Q16 A：保留原PRnet表达目标与结构；训练输入log尺度匹配对照、预测处理表达，validation/test/最终推断输入固定训练CM共识，并由预测表达减该共识得到Δ。不采用signed Δ输出头改造。
17. 用户按推荐确认Q17：每折训练许可范围的合格DMSO孔先按板求log表达均值，再跨板等权；保存对照membership与权重，最终全量训练后生成交付共识。
18. 用户澄清Q18无需另作剂量选择：测试集只来自审定后的temp_env，自然只含0.2/1/10µM。按药物身份切分，保留入选药物通过QC的全部实测剂量；本轮衡量已测三剂量下的CM未见药物泛化。未来任意剂量的预测交付不属于本轮测试范围。
19. 用户按推荐确认Q19：CM微调更新PRnet全部参数，不冻结预训练层；部分冻结留作后续可选消融。
20. 用户修订Q12/Q20：方向AUPRC@200为验证主指标，CM微调按验证AUPRC@200最高选择checkpoint并早停，Δ RMSE为辅助。训练损失保留Gaussian NLL，validation/test使用固定训练共识对照；不再按Δ RMSE或验证NLL选择CM微调checkpoint。历史bring-up的NLL记录不改写。
21. 用户按推荐确认Q21：正式Tahoe预训练同样按验证方向AUPRC@200最高选择checkpoint并早停，训练损失仍为Gaussian NLL。Tahoe的预测/真值Δ用其plate×cell_line匹配对照，不套用固定CM共识。
22. 用户按推荐确认Q22：完成预先约定的对照实验即可结束首轮，接受正结果或负结果；预训练相对随机初始化的增益与预测器相对简单基线的可用性分别报告，不要求反复调参直到成功。
23. 用户确认Q23 A：3个随机种子，固定相同五折；Tahoe独立预训练3次，对应CM两种初始化各进行五折×三种子评估，不择优报告种子。
24. 用户确认Q24 A：双方各试1e-4/3e-4/1e-3三个学习率，搜索预算和epoch上限规则相同，按验证AUPRC选择并允许独立早停。首轮3个Tahoe模型、90个CM候选训练，选定配置后保留30个CM测试单元。

## Q18前置事实：剂量与孔数（实际元数据复算）

以下均为新增药理资格及实验QC过滤前的数量；物理孔和pseudobulk观察不能混用。

### Tahoe

用`pert_data/Tahoe-100M/source/metadata/sample_metadata.parquet`的plate×sample去重，并与`stage2/04_tahoe_data/output/pseudobulk/manifest.csv`按plate×well核对，剂量0个不一致：

| 剂量µM | 物理处理孔 | pseudobulk观察（孔×细胞系） |
|---:|---:|---:|
| 0.05 | 424 | 21,043 |
| 0.5 | 402 | 20,079 |
| 5 | 490 | 24,496 |
| 0：DMSO_TF对照 | 28 | 1,400 |
| 总计 | 1,344 | 67,018 |

379种化合物均覆盖三个剂量，但重复数不统一：1,137个drug×dose条件中975个1孔、158个2孔、3个3孔、1个16孔。不能把每孔约50系当成50个独立处理重复。

### temp_env

原始`GSE262419_hash.csv.gz`以Plate_ID×Well_ID计数，共4,558唯一孔、0重复；按Chemical_Concentration_uM分类：

| 剂量µM | 类型 | 物理孔 |
|---:|---|---:|
| 0.2 | 小分子处理 | 1,402 |
| 1 | 小分子处理 | 1,402 |
| 10 | 小分子处理 | 1,402 |
| 0.01 | ET-1阳性处理（非小分子） | 48 |
| 0 | VEH/DMSO对照 | 192 |
| 0 | Media对照 | 80 |
| 0 | Blank技术空白 | 32 |
| — | 总计 | 4,558 |

`temp_env_all-smiles.h5ad`中三个处理剂量仍各1,402孔，共4,206处理孔；保留272对照，排除ET-1及Blank，共4,478行。现`atc-moa`113个key视图每剂量349处理孔，不能视为未来审定后的药物数量。

一般 CM 响应预测不能自行升级为缺氧/复氧条件下的疗效预测。

## 事实基线与历史边界

- ADR 0026 冻结 CM corpus/v1：18,164 基因、525 个 parent drug_key。已有五个 datasets-1 视图。
- ADR 0029/0031：Tahoe pseudobulk 为 67,018 个 `(plate, well, cell_line)` 表达观察；CM-schema view 为 18,164 基因。
- 2026-11-01 已完成一次 PRnet bring-up。control 是 raw counts 均值，target 是 log1p_cp10k；尺度不同。该 checkpoint 仅证明管线跑通，尚无迁移增益证据。不得静默把它改称正式迁移预训练。
- 现有 18k 基准没有 entrant 超过 delta_mean；PRnet 的 5k 敏感性结果不能代替 18k 主结果。正式迁移实验应保留随机初始化及常数/线性基线。

依据：[Stage 2 README](../README.md)、[ADR 0026](../../docs/adr/0026-stage2-cm-data-contract.md)、[ADR 0031](../../docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md)。

## 湿实验与当前加工：须逐项区分

下表列的是仓库输入/已加工数据可核查的情况，不能当成完整实验 plate map。物理板型、原始重复与剂量安排仍须逐项回接论文、补充材料与设计文件。

| 来源 | 观察单位与暴露 | 当前输入/视图 | 对照核查发现 |
|---|---|---|---|
| Tahoe-100M | 每孔混合约50种癌细胞系；按plate/well/line聚合单细胞，非50个独立培养孔 | 14板，700个plate×line对照组；1,400对照观察、65,618处理观察 | 预训练view已按plate×line匹配，不能跨line平均；每组两条对照观察不等于两个供体生物重复 |
| dtoxs_primary / LINCS-LDS-1593 | 原代人CM，48h；剂量含文献来源 | 原始可解析输入677处理孔、102对照 | 10 DMSO、2 EtOH、90 VEH；VEH须回接设计含义。原始Primary_Cell_Metadata含A/B/D/E细胞来源，loader的donor=NA不代表上游没有信息 |
| dtoxs_ipsc / GSE217424 | iPSC-CM，48h；series matrix剂量 | 296处理孔、16对照 | DMSO；同一plate6存在MSN02/MSN09，现plate均值可能混供体 |
| screen_ctox / GSE244740 | ScreenSeq hiPSC-CM，24h，8点半对数剂量 | 2,520处理孔、552对照；5块Optimization板+3块Screening板，板号须带assay | 536 DMSO、16 NoTreatment；不同对照类型不应未经审定混合 |
| temp_env / GSE262419 | TempO-Seq hiPSC-CM，24h，3剂量 | 16板；4,206结构可解析处理孔、272对照 | 每板12 DMSO+5 Media；另2 Blank/板被排除。现plate均值混合DMSO和Media |
| GSE247541 | 抗病毒药物CM面板，144h | provenance-only，未纳入现CM corpus | 排除沿用了历史时间映射约束；对本轮一般CM目标是否有增量价值待重新评估 |

更具体的配板证据：Tahoe为14×96孔、0.05/0.5/5 µM，每板2个DMSO_TF孔；temp_env为0.2/1/10 µM，原始hash表有4,558行，每板另有3个ET-1阳性处理。ScreenSeq为384孔、triplicate、0.5% DMSO并使用sunitinib阳性处理；原始temp_env坐标属384孔坐标体系，但板规格尚未从methods直接核实。

已下载不等于已使用：当前DToxS-iPSC只读取GSE217423的DGE分支（原研究41药、2个CM系）。同SuperSeries的GSE217421 Conv分支含54药、6个CM系、48h常规RNA-seq，尚未进入当前Stage 2 corpus。它是“是否需要新数据集”必须核查的候选，而不仅是在现有三种temp_env tier之间选择。

原始设计锚点：`pert_data/GSE262419/source/GSE262419_hash.csv.gz`；`pert_data/GSE244740/source/GSE244740_series_matrix.txt.gz`；DToxS两分支各自series matrix及DGE实验设计表。官方来源：[GSE262419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262419)、[GSE244740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244740)、[GSE217423 DGE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217423)、[GSE217421 Conv](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421)。

关键缺口：确认原始供体/细胞批次、溶剂、板位置、培养与技术重复的关系；将原始重复编号和后处理cumcount区分。temp_env原始hash表含Replicate_Number与Solvent，现obs.replicate是生成索引，不能据此证明生物重复数。

本地依据：`pert_data/*/source/` 原始元数据；[02 README](../02_cm_data/README.md)；[views.py](../02_cm_data/execution/cmdata/views.py) 的control_group及control_mean计算；[04 README](../04_tahoe_data/README.md)。Tahoe论文：[Cell原文](https://www.sciencedirect.com/science/article/pii/S0092867426010081)。论文报告时间与现view缺失time字段需通过明确出处补证，不能直接填常数。

## temp_env：化合物数、药理类别数、训练价值数

对已冻结 [compounds.csv](../02_cm_data/output/datasets-1/compounds.csv) 选 `temp_env_all-smiles=True`，按drug_key去重：467条名称记录对应460个key。

| 视图 | 化合物key | 处理观察 | 对照观察 |
|---|---:|---:|---:|
| atc-moa | 113 | 1,047 | 272 |
| atc-moa-mesh | 252 | 2,310 | 272 |
| all-smiles | 460 | 4,206 | 272 |

106个key有ATC，84个有非空MoA，77个两者都有；29个ATC-only、7个MoA-only。提取完整ATC编码的正则为 `\b[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}\b`，按前缀长度1/3/4/5/7统计，L1—L5分别为14/58/88/124/172个不同编码。这是多标签注释覆盖，不是互斥类别或药物数量。建议药理类别先报告L3的88类，并同时报告106个key的覆盖。

MoA有60种不同完整字段字符串；分号拆分得到72个片段，但分号也用于连接target/action，故不能称72种规范机制。

后续校正：上述ATC的14/58/88/124/172仅由完整7位编码向上投影。冻结字段另含22个仅到L4的编码（如Cisapride A03FA）；把这些部分编码也计入其可支持的层级后，L1—L5为14/60/92/134/172。正式类别表须规范化并分别报告完整/部分编码口径，不能继续将88称作所有注释覆盖的L3总数。

现有 `bool(atc or moa)` 筛选把 `Unknown` 也视为有机制。MoA-only含Indoxacarb（External insecticide）、Methyl salicylate与Coumarin（Unknown）等。因此113是自动tier数量，尚不能回答“多少种真实药物”。参见 [drugs.py](../02_cm_data/execution/cmdata/drugs.py)。

训练价值合格数目前**未知**。结构与非零计数QA不等于可学习性：

- 每key有3剂量；442个key各9孔、16个各12孔、2个各18孔。
- 459个key只出现于一板，1个跨两板；药物与板高度耦合。
- 处理孔lib_size最小6,423、中位887,947；检测基因最小2,573、中位10,533，不能只用非零判定高质量。
- 应审定结构、同载体同生物背景对照、读数/活性QC、重复一致性与剂量响应；可重复的低响应也是训练信息。只保留强DEG会改变目标分布。

建议：保留冻结corpus与datasets-1；在新训练输入中增加可追溯资格/质量标记。药物聚焦以审定后的atc-moa为起点，all-smiles作为覆盖消融；MeSH视图不能直接解释为更多药物。

按已确认的研究工具纳入规则，资格审核必须覆盖全部460种结构身份，不能只审核原113。其余347种中出现Cyclopamine、Genistein、Candoxatril、Azimilide、Sematilide、Oxytetracycline hydrochloride、PK11195、Pirinixic acid等须复核候选；旧自动tier存在假阴性。这些是审核候选，不是已完成逐项证据审核的最终名单。

## 统一h5ad：已有契约与待审定部分

已有公共obs 24列：

```text
well, dataset_id, compound_name, drug_key, inchikey, smiles,
dose_um, dose_source, time_h, is_control, vehicle,
plate, control_group, cell_type, donor, replicate,
atc, moa, mesh_classes, pharm_atc_moa, pharm_mesh,
annotation_source, lib_size, n_genes_detected
```

- var：18,164个固定顺序HGNC symbol；entrez_id、ensembl_id、match_type。
- X：原始int32 counts；layers['log1p_cp10k']：同一保留基因轴上按行归一化的float32表达。
- uns['control_mean']：现契约为raw counts均值，不能静默改变语义。
- Tahoe额外保留cell_line_id、cell_name、n_cells、moa_broad、moa_fine、pubchem_cid、has_structure、has_matched_control；作者兼容列SMILES/dose/paired_control_index属于适配层。

推荐待确认：一行统一为可追溯表达观察，另保留物理孔身份与细胞背景；Tahoe唯一观察键须含line，不能以well单列要求唯一。训练对照取匹配对照行的 `mean(log1p_cp10k)`，与target同尺度；不把它偷换成 `log1p_cp10k(mean(raw))`。若缓存，使用带尺度的新字段。

对照应按同plate、同载体、同供体/细胞系、同时间等实际设计定义；无法回接时列为缺失/降级，不捏造。缺失基因应有测量掩码：Tahoe的FAM153B全零槽是缺测，不能解释为实测零表达。split、质量、原始重复与来源字段的最终位置待确定。

历史边界补充：CM h5ad的raw control_mean混合对照，并不意味着03.1历史benchmark也使用混合对照。其execution/script/build_data.py的build_control_bundle明确只取每板12个DMSO行，并在log层重算均值；新协议须区分h5ad缓存、03基准实际输入和05 bring-up实际输入。

## PRnet基因空间与容量

### 剂量输入机制（用户追问核查）

作者`pert_model/PRnet/source/data/_utils.py`的Drug_dose_encoder先算1024维radius=2/useFeatures=True的FCFP4二进制指纹f，再形成`rFCFP=f*log10(d+1)`，经药物适配器进入网络；不是剂量分箱或独立dose embedding。当前03 adapter与05 pretrain都沿用此式，并将d解释为dose_um；编码函数本身不自动识别/换算单位。

按µM输入，0.05/0.2/0.5/1/5/10的缩放分别约0.02119/0.07918/0.17609/0.30103/0.77815/1.04139。同药不同剂量沿同一指纹方向改变幅度，后续非线性网络预测响应；这不施加响应单调、线性、饱和或EC50规律。

因此Tahoe与temp_env剂量不重合不构成输入接口障碍，不需强行分箱对齐；但支持输入任意非负剂量，不等于在未测剂量上已验证可靠。尤其10µM超过Tahoe预训练最高5µM，却是temp_env已测微调剂量，不能把它称为整个流程都未见的剂量外推。用户已澄清测试剂量由temp_env实测范围自然确定，不再把任意剂量交付作为本轮待决策问题。

论文描述FCFP随剂量缩放：[PRnet正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11513139/)；精确log10(d+1)公式以作者代码为准。

24h范围更新后的基因事实：对source_gene_sets.csv按hgnc_id求交，screen_ctox为31,880，temp_env为19,613，两者共同18,324；旧18,164全部包含，新增160（0.881%）。新增160经HGNC→Ensembl连接均被Tahoe测量，新轴仍仅FAM153B缺测。新增160含112 protein-coding、9 non-coding RNA、5 pseudogene、34 other。保留旧轴有基准兼容收益；这些数量不足以证明新增基因没有疾病相关价值，轴选择仍待确认。

[PRnet官方论文](https://www.nature.com/articles/s41467-024-53457-1)及[PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11513139/)：L1000直接建模978 landmark，sci-Plex3建模5,000 HVG；12,328基因来自978输出的线性推断。Table1为hidden128、drug latent64、state latent64、noise10；constructor默认值与论文表不完全相同。

本地同结构只改变基因维度时，参数量为 `386*g+173696`：978→551,204；5,000→2,103,696；18,164→7,185,000。它们是本地结构计算值，不是论文报告的参数量。18k下只改state latent128→7,201,448；hidden256→14,194,280。

建议先固定当前128/64结构，在修正配对与表达尺度后比较随机初始化/Tahoe初始化，再独立消融state latent128、hidden256。不能仅凭基因数决定加层；模型选择只用验证集，测试集不能参与容量选择。

正式CM微调已确认全参数更新。训练优化Gaussian NLL；模型选择与早停监控固定共识输入下的验证方向AUPRC@200（越高越好）。方向分别计算AP再宏平均，条件为drug×dose，先汇总重复再在条件间宏平均；k=200标签沿用固定top/bottom排名定义。Δ RMSE、Δ相关性及低响应/毒性分层作为辅助报告，不能替代已选主指标。该排名标签不等同统计显著DEG，低响应样本仍按已确认原则保留。容量消融的模型选择也须遵循相同验证主指标。

## 微调来源与泛化：待决策

**当前范围覆盖声明（Q15最新修订）：仅Tahoe＋temp_env。下文涉及ScreenSeq的数量与方案是已核查事实和历史备选，不属于本轮训练/评估；本轮取消ScreenSeq外部验证与两来源合并消融。** 其他非24h来源仍全部排除。

本轮只保留24h：CM候选来源为screen_ctox与temp_env，均为iPSC-CM。48h的DToxS primary、DGE、Conv与144h的GSE247541均不训练、不验证、不测试；此前关于其辅助训练/外部验证的建议不再适用于本轮。Tahoe的24h需引用正式论文补证，现obs.time_h=NA不能直接通过24h字段过滤。

冻结注释复算：screen_ctox的42个key与temp_env atc-moa的113个key交集为12，并集143；与temp_env all-smiles的交集也是这12个key。ScreenSeq独有30个key，包含Doxorubicin、Amiodarone、Digoxin、Sunitinib、Dasatinib、Lapatinib、Imatinib等。143是现自动tier并集，不是审定后的最终药物数量。

来源方案已确认：temp_env主实验、ScreenSeq独立来源验证、合并微调消融；对应实验都应包含random-init/Tahoe-init对照并分别报告来源表现。不能把同一药物在不同来源的孔拆到train/test两边后称cold-drug。

ScreenSeq与temp_env现有药理tier共12药：Aspirin、Cisapride、Zidovudine、Diltiazem、Nifedipine、Atenolol、Clozapine、Acetaminophen、Verapamil、Bepridil、Sotalol、Furosemide。独立来源cold-drug验证不能直接让temp_env全量训练后的模型测试这12药；其用途与排除/按fold预测规则需另轮明确。ScreenSeq独有30药可作为候选外部未见药物集，须经过同一资格/QC规则。

Tahoe与CM corpus按parent drug_key重叠31/525。现Tahoe plate13/14验证并非cold-drug，不能用于证明未见药物泛化。

Q3/Q4已确认，按用户选择采用CM未见口径，不排除Tahoe中的同药物；测试药物在所有参与微调的CM来源中保持身份隔离。可另按Tahoe见过/未见分层检查增益，但不改变主评估资格。

下一轮访谈前沿：

Q21—Q24均已按推荐确认。主要研究决策已收拢，用户已批准按[统一spec](../../.scratch/stage2-prnet-formal-transfer/spec.md)拆票；[11张票据及依赖](../../.scratch/stage2-prnet-formal-transfer/ticket-plan.md)已发布。逐化合物资格、样本QC与最终计数仍需事实审计，不能写成已完成。

1. 共识构建权重已按Q17确认：每折训练许可板的合格DMSO行先按板求log表达均值，再跨板等权；唯一对照行及板权重均保存，不按该板药物孔数重复计权。不能使用独立held-out板对照形成共识。
2. 评估仅限temp_env：drug_key五折6:2:2，所有剂量/重复同行，初始化间切分配对。若采用固定共识推断，内部测试也不能读取测试板对照作为模型输入；匹配对照仅用于形成实验真值。共识仅由许可的训练对照构建并冻结，详细控制行池须在切分协议中列明。
3. 测试剂量已由temp_env自然确定为0.2/1/10µM，无需额外选择或单独设计剂量外推测试。资格审核需覆盖460种而非固定113。ScreenSeq外部主结果与共享药诊断已撤出当前范围。

毒性/低响应的保留原则已确认；质量阈值、标签证据与样本权重尚未确认。响应标签仅用于评价/描述时，不得反过来作为测试药物响应输入。

固定共识的可审议实现：

- 方案A（用户Q16已接受）：保留PRnet表达目标，训练输入为log尺度匹配对照c_g，target为处理表达y；validation/test输入固定训练CM共识c*，预测Δ=F(c*,d)-c*。有条件输入偏移，必须在validation即采用部署方式，不能用匹配对照输入的分数替代。对应预测是共识CM背景效应，不承诺供体/批次特异响应。
- 方案B：CM训练及预测均输入c*，目标为实验Δ=y-c_g，模型需允许正负均值输出，区别于原PRnet表达目标/均值ReLU；Tahoe也须重训为相同Δ语义，以其plate×line匹配对照计算标签。不得跨50系平均Tahoe对照。
- 不推荐折中：训练输入c*但target仍为各板原始y，这会把不可观察的板基线差异混入拟合目标。
- c*必须保存训练对照行membership、归一化定义及聚合权重；held-out环境对照只可用于独立计算Δ真值，不进入模型输入。对照池划分与药物折的关系仍须明确。

## 决策状态

已确认：CM未见药物、人源iPSC-CM一般响应、严格24h、仅Tahoe＋temp_env、机制工具纳入、低响应/毒性响应QC合格保留、18,164轴、方向AUPRC@200主指标及Tahoe/CM的checkpoint与早停依据、Δ RMSE辅助、先固定结构再容量消融、全参数微调、原PRnet表达目标与Gaussian NLL训练损失、匹配对照训练/固定共识验证及推断、板内均值后跨板等权的训练共识、3种子及双方相同三学习率搜索预算、首轮接受预先约定实验的正负结果；测试自然覆盖入选temp_env药物的实测三剂量。原ScreenSeq外部验证/合并消融已撤出本轮。用户已批准拆票并完成发布；逐化合物资格、逐样本QC和具体数据gate待执行票据完成，正式训练尚未开始。
