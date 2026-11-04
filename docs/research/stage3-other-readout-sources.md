# Stage 3：其他心肌细胞读数的 500 个独立组合门槛核验

研究日期：2026-11-04。对应 [核验其他心肌细胞读数能否达到 500 个独立药物表型样本](https://github.com/Moke110/arc/issues/23)；边界见 [ADR 0034](../adr/0034-stage3-phenotype-data-boundary.md)。这是有限检索的可行性研究，不是生产数据清单。

**后续核验：** 两个来源的真实逐剂量附表现已取得，实测计数、缺失/零、身份与对照边界见[原始附表样本数核验](stage3-sample-count-audit.md)。下文保留初查时的证据状态；早先统一 0.3% DMSO 的描述也已在后续核验中按 Sirenko Methods 澄清。

## 结论

**本轮没有逐行核实出 ≥500 个可取得、去重后、时间窗合格的单指标组合；不能据此宣称达标，也不能据此断言不存在。** 有两个优先核验对象：

1. **Sirenko 2013 的 24 h 钙瞬变所导出的 beat rate**：原文设计支持 131 种药物 × 4 个浓度 = **524 个名义组合**，时间窗最直接匹配，但缺少补充表逐行完整性核验。
2. **Heliste 2021 的 HL-1 / CellTiter-Glo / 缺氧复氧筛选**：689 种 FDA 药物、24 h 预给药后进入 HR；属于与 MIRI 接近的候选。实际末端读数在开始给药后约 96 h，HR 期间是否持续药物和完整单药数据尚未核实，因此只有条件性可行性。

建议 research resolution：**初步研究已回答，门槛未证实；保留两个候选进入针对性小表核验，读数选择仍开放。** 不把库大小、重复孔、多个成像特征、不同读数的总和当作达标证据。

## 计数与可比较性

本报告的设计组合数是药物 × 非重复条件 × **一个指定 readout**。浓度可以区分条件，但不增加药物种数。重复孔先合并，不能把相同药物/剂量的两孔计为两条。不同物种、细胞来源、HR方案和暴露时长分别保留；跨研究累加之前仍需确定条件等价和药物身份规则。以上来自 ADR 0034 的边界，尚不形成新的纳入决策。

“设计数”只证明实验名义容量。实际可用数还要检查缺失、药物别名/盐型重复、无法读取的信号、对照及效应数值的可得性。beat rate 的停搏值和无法拟合峰形的缺失值也需要区分，不能自动删除或自动标记为 no change。

## 候选证据

| 来源与单指标候选 | 药物数、组合与重复 | 时间与模型 | 当前判定 |
|---|---|---|---|
| Sirenko 2013：beat rate，基于胞内 Ca²⁺ 荧光动力学 | 131 drugs；24 h **4** 浓度，设计 524 组合；duplicate wells 不计额外样本 | human iPSC-CM；24 h，无 HR | 时间匹配，设计门槛匹配；实际表格去重/有效值数未证实 |
| Heliste 2021：CTG ATP 发光活性 | 文称 689 FDA compounds；4 浓度；单 HR 分层设计上限 2,756，不能作为实数；重复孔数未核实 | mouse immortalized atrial HL-1；24 h 预处理 + 48 h hypoxia + 24 h reoxygenation | 可再 HR 定义下的强候选；不是给药后24–48 h直接读数 |
| Sirenko 2017：beat/peak frequency 或 Calcein AM，必须分别计 | 共69 environmental chemicals **and** drugs；不能将69都当药物；6浓度时每个时间点/指标名义上限414；duplicate 不增加组合 | human iCell iPSC-CM；30 min / 24 h | 单研究全化合物设计上限仍不足500；药物子集更小 |
| Burnett 2021：ToxCast Ca²⁺ beat rate | 1,029 chemicals，表1区分125 pharmaceuticals与28 CiPA pharmaceuticals；不是1,029 drugs | 5 human iPSC-CM donors；功能读数90 min | 时间窗不合格，不能用规模抵消 |

表内来源分别为 [Sirenko 2013 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3900303/)、[Heliste 2021 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8371130/)、[Sirenko 2017 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5734940/)、[Burnett 2021 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8762671/)。表中乘法是本报告依据设计的推算，**不是附件行数实测**。

### 1. Sirenko 2013：最直接的 24 h / 524 组合设计

原文 Results → “Analysis of the concentration response”明确区分：30 min 用6浓度，**24 h用4浓度**，均做duplicates。摘要的“up to 6”不能套给24 h。候选优先锁定 beat rate；峰宽、幅度等另列，不能叠加。0.3% DMSO vehicle 为基线；每板30–38个vehicle孔，原文采用相对对照及BMC分析。131药物与四浓度若全部有可用beat-rate值，设计数524；余量仅24，少25个组合就不能达500。[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3900303/)，[出版社页面/摘要](https://www.sciencedirect.com/science/article/abs/pii/S0041008X13004171)

数据入口是原文 Supplemental Table 1（药物身份）及 Tables 2/3（两个时间点的数据线索）；本轮未成功取得这些小表，**未确认 Table 3 是完整逐剂量值还是处理后的模型汇总，更未验证524行均存在**。只有每药一个BMC时，131个BMC不能当524个逐剂量表型值。许可也须单独核实；出版社页面保留版权，PMC可阅读不自动等于数据可任意再分发。[PubMed稳定标识 PMID 24095675](https://pubmed.ncbi.nlm.nih.gov/24095675/)，[出版社版权页](https://www.sciencedirect.com/science/article/abs/pii/S0041008X13004171)

### 2. Heliste 2021：HR 与 CTG 的条件性候选

Drug library screening 方法：689 FDA-approved compounds；40、200、1,000、5,000 nM，0.01% DMSO；HL-1小鼠心房永生化细胞，24 h后转48 h低氧、再复氧24 h。CTG测ATP相关发光活性；不能视为MTT或细胞数。方法未明确报告此筛选的洗脱或补药，**不能推断仅暴露24 h**；需记录预处理24 h、总历时96 h、HR期间药物状态未知。[原文 Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8371130/)

原文按DMSO归一化，并从四剂量选最大HR/常氧相对RLU作为每药代表值。因此选后数据有剂量选择偏差；必须核对是否能取得各剂量及同方案DMSO对照，不能直接把选后值用于固定剂量效应。论文声明数据在文内和补充材料，文章标注CC BY 4.0（第三方内容另核）。**689只是论文宣称的库大小，本轮没有证明689个独立药物都有一个可提取效应值。**[原文 Methods、Data availability、License](https://pmc.ncbi.nlm.nih.gov/articles/PMC8371130/)

小文件入口：[Supplementary Information 2（99.5 KB xlsx）](https://pmc.ncbi.nlm.nih.gov/articles/instance/8371130/bin/41598_2021_96033_MOESM2_ESM.xlsx)，另有 Information 3/4，须确认哪个工作表对应 Dataset S1。此次publisher CDN连接失败，PMC地址返回1,816字节“Preparing to download”HTML而非xlsx；已检查文件类型，**没有解析伪附件或报告虚构行数**。因此重复孔数、身份映射、缺失率、误差列和准确组合数均未核验。HR可再经历是否允许持续药物至96 h，留给时长定义决策，不能默许纳入。

### 3. 2017 小库及较大 ToxCast：为何暂不能补齐

Sirenko 2017覆盖69种环境化合物和药物，human iCell iPSC-CM；30 min与24 h测Ca²⁺ flux，并有高内涵活性/线粒体读数。方法给出0.3–100 µM半对数稀释、duplicates；按6浓度推算单时间点单指标最多414，且仅部分为药物。Supplementary Tables含POD及其他数据文件，未核验逐剂量全表；不能累加多个终点跨过500，也不能与2013直接相加而不查重。[原文 Methods、Supplementary material](https://pmc.ncbi.nlm.nih.gov/articles/PMC5734940/)

Burnett 2021大库确有1,029化合物、5个人源donor，四浓度0.2、1、10、100 µM；Ca²⁺功能读取在90 min。表1的药物类别最多153项也需身份去重，其余多为农药/工业化学品。**对本次24–48 h钙瞬变问题不合格**；没有将其全部高内涵指标自动标记为24 h。本轮未审计官方ToxCast发布包，不将“在ToxCast存在”当作已取得合格样本。[原文 Table 1、Functional cardiotoxicity screening](https://pmc.ncbi.nlm.nih.gov/articles/PMC8762671/)

电生理/收缩的线索如CiPA或多donor数据，应先检查真实暴露时长；本轮未确证一个满足24–48 h和500独立组合的此类附件。131-drug钙瞬变设计无需更换为混合指标即可接近门槛，优先解决其附件核验比扩大模糊候选更有价值。未覆盖所有动物原代、人原代和心肌样细胞系研究，不能将有限搜索的阴性当成这些层没有数据。

## 下一次小表核验的可验收结果

1. 取得2013 Table 1/3，记录版本/许可；输出药物身份去重前后数、每药24 h浓度数、一个指定beat-rate列非缺失组合数、重复结构和对照列。以唯一键计数，明确停止跳动的编码；只拿到BMC则报告131级别的剂量曲线汇总。
2. 取得2021 Dataset S1对应工作表，核对689是否含别名/重复品，输出逐剂量CTG值、HR与常氧各自DMSO值、误差/重复及选后代表值是否可逆恢复。若仅hits或选后分数可得，不宣称完整2,756条。
3. 在时长决策中明确“24–48 h预处理再HR”是否允许药物在HR期间持续存在。记录总给药时长；CTG是ATP代理活性，不能默认映射到MTT/CCK8或采用作者显著性阈值作三分类标签。

本轮只读取论文、元数据和尝试小附件；未载入大型数据，未生产入库，未训练模型。
