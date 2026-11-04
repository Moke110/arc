# Stage 3 正式表型 corpus：已同意设计与 spec 交接

日期：2026-11-04。状态：`grill-with-docs` 访谈已完成；本文件交接已同意设计，
不是已发布 spec、生产数据、已运行验收或模型结果。

依据：[ADR 0038](../adr/0038-dedicated-phenotype-data-and-stage3-corpus.md)、
[ADR 0039](../adr/0039-condition-corpus-and-model-view-separation.md)、
[ADR 0040](../adr/0040-numeric-endpoints-and-source-readout-sequences.md)、
[ADR 0041](../adr/0041-formal-corpus-statistics-export-and-coverage.md)及
[GLOSSARY](../../GLOSSARY.md)。来源事实沿用[原件审计](stage3-sample-count-audit.md)、
[来源纳入决定](stage3-inclusion-decision.md)与[原整理契约](stage3-curation-contract.md)。
新决定优先于旧 wayfinder 文件中的冲突设计：来源归属、步骤名、交换格式、样本键、
MTT/CCK8 phenotype、数值标签生成、0 的含义和 500 验收位置均以本轮决定为准。

## 1. 范围与来源顺序

当前实施范围是实验表型来源正式归档、细粒度条件整理、证据与统计处理、数值 endpoint
派生、corpus 发布及离线验收。后续 model data view 的构建、条件合并/聚合、模型训练、
benchmark 和 Stage 2 预测扰动组生成/匹配均不在范围内。

来源顺序保留：Sirenko 2013 `24hr_beat_count` 为第三类优先来源；A/Linn →
B/Liu icariin → D/Duan SYA/NAC 为 MTT/CCK8 图级试点。各细胞模型、常氧/HR 和
给药顺序分别保存；C/NR 时长冲突保持隔离，D 的 4 h 常氧不冒充目标窗口。
Heliste S1 在持续给药状态、归一化分母等条件补证后才是可重新验收的候选。
不通过扩指标、改时间或混入非药物参比来凑规模；不主动联系作者。

来源可归档、测量可提取、记录可进入 corpus 是三个独立结论。不能判 endpoint 的来源
仍可正式归入 `pheno_data`；这不使其成为带标签 corpus 样本。

## 2. 未来目录与不可变证据

新增第五外部来源类别 `pheno_data/`，一个上游数据产品一个项，未来布局为
`pheno_data/<item_id>/`。遵循 [ADR 0030](../adr/0030-external-data-provenance-and-step-layout.md)：

- `script/` 保存下载及离线验证入口；`source/` 保存未经改写的上游字节、名称和目录。
- 裁图、数字化值、身份映射快照等本地新增材料放 `local/`，逐项登记来源、规则、工具版本和哈希。
- 跟踪 `author_log.json` 与 `manifest.sha256`；前者记录 URL、版本、获取时间、许可事实、验证及 local additions，后者只覆盖 `source/`。
- 镜像、重复出版、同一文件的多个端点不制造新的实验来源；更正版保存版本关系，不覆盖旧原件。

未来消费步骤是 `stage3/01_data_corpus/`，配套 Stage 3 与步骤的 `AGENTS.md`、`README.md`；
步骤沿用 `execution/`、`output/` 结构。交换产物是 **每个 `(condition, phenotype)`
一个 Parquet 文件，加 JSON 清单及 schema**，支持真实 list 字段。
发布包还需包含规则/统计参数、证据关联、排除记录、覆盖与损耗报告、文件哈希和验证结果。
具体文件命名与目录编码由 spec 固定，当前不创建这棵生产目录。

每个 release 不可变；修正产生新 release 并关联旧版。发布记录固定输入哈希、代码 commit、
schema/规则版本、统计软件与参数、输出哈希和运行证据。输出包的校验清单与 source manifest
分开；Stage 3 引用原件路径与哈希，不复制原件制造第二份来源。逐资产保留许可事实；
未知许可不转写为开放许可，受限原件不因本设计自动获得再分发权。

## 3. Condition、phenotype 与样本行

Condition 采用 `condition_*` 字段，至少保留 `condition_cell_type`、
`condition_drug_time`、`condition_hypoxia_time`、`condition_reoxygen_time`；进一步保存
物种、细胞来源/模型、培养背景、氧方案、vehicle、检测时钟及其他影响解释的实验条件。
同名细胞或同名检测不自动证明条件等价；关键未知不得以 `null == null` 当作合并证据。

`condition_drug_time` 在 HR 实验中指 HR 前预处理时长，无 HR 时指直接暴露时长。
`condition_drug_schedule` 另记 HR 期间持续、洗脱、补药或未知；总历时和实际暴露时钟另存。
已知无 HR 时 hypoxia/reoxygen 时间为 0；不知道是否/多久 HR 时是 null。单位规范化后
仍保留原值、原单位和换算证据。不得用总历时替代预处理时长。

MTT 与 CCK8 的 phenotype 均为 `cell_viability`，逐值保留 `assay_method`。
该决定不自动将 CTG、Calcein 成像或细胞计数映射为 `cell_viability`。Beat count 保留测量窗口和原单位，
没有证据不转为 bpm。

一个文件内唯一键为 **`(source, drug, dose, endpoint)`**；加上文件所属 condition/phenotype
才是完整样本身份。`source` 是上游数据集，不是每个附件、板、孔或重复；`drug` 是证据解析的
活性母体，给药盐型等形式单独保留；`dose` 的规范值必须有明确单位及计量基准。
用途资格、母体/别名去重和有效测量要求沿用来源决定，不凭 CAS 或名称唯一判为合格药物。

| 字段/关联 | 最小语义与约束 |
| --- | --- |
| `sample_id`、`source`、`drug`、`dose`、`endpoint` | 稳定标识与规范键；endpoint 仅数值 −1/0/+1，不可 null |
| `condition_*`、`phenotype` | 文件清单与行内容一致；未知条件附状态/理由，不能隐藏成已确认条件 |
| `readout` | 非空数值 list，保存所有合格真实来源测量，不以均值替代 |
| 逐值元数据/关联 | 与 readout 元素稳定对应的 value/observation ID、单位、assay_method、表示类型、原实验/批次、重复号/类型、来源定位和 comparison ID |
| normalized 表示 | 与原测量 ID 关联另存，不作为额外重复；normalized-only 标明表示并保留原数值 |
| comparison 证据 | 处理/对照成员、匹配依据、均值与 FC、p/p_adj、检验与 family、输入来源、endpoint 派生版本 |
| source/evidence/decision 关联 | 原件哈希、sheet/单元格或页/图/面板定位、身份/条件/有效性判定、排除理由及修订历史 |

这些是语义契约；嵌套 list-of-struct 或数值 list 加关联表的具体 schema 由 spec 固定。
若使用并行 list，长度和次序必须一致；不能只留一串数字而失去逐值证据。
一条来源均值仍是一条 summary 测量，即使作者写 n=6 也不得展开六条虚构测量。
SD/SEM/CI 作为统计证据关联，不混入处理 readout list 冒充重复。

优先 raw 作为 readout，作者 normalized 另存；仅有 normalized 时保存实际 supplied 值并标类型，
不制造 raw。raw 与 normalized 是同一测量的不同表示。raw 空白而 normalized=0 时保留两种事实，
不得据此凭空创造 raw 零或独立重复。

**先在原实验比较层计算统计和 endpoint，再按最终键合并。** 同源、同 condition/phenotype、
drug/dose/endpoint 的记录合并为一行并保留全部真实 readout 与比较关联；不同 source 或不同
endpoint 不合并。相同 drug/dose 出现相反 endpoint 时保留多行，冲突及后续聚合归 model view。
不能先按 endpoint 分组来决定检验输入；同一个测量值的多重表示/镜像不重复计入 list。

## 4. 对照、检验与 BH family

FC 定义为 `treated mean / matched-control mean`。分子分母须来自明确、量纲兼容且有效的
对应比较；对照均值必须为正且有限。保留对照的条件、vehicle、时间、归一化语义及来源。
normalized 值只有在分母/尺度明确时才支持 FC；不凭“百分数”自行假定对照为 100，
不使用零分母 pseudocount，不从 Heliste HYP/CTRL 自造 HR 比值。

优先使用作者明确对应比较的**数值 p**，保留比较对象、检验、调整方法与 family 证据。
作者原始 p 进入规定 BH 流程；作者已提供符合所需规则且有证据的 BH-adjusted p 时直接保留，
不得再次 BH。未说明调整方式/比较对象的数值不能自动当作满足规则的 p_adj。
星号或 `ns` 不能替代数值 p，也不反推一个数值。

需要重算且输入与假设适用时：独立组用双侧 Welch t-test；来源明确配对时用双侧 paired t-test。
记录配对关系和独立实验单位；技术孔不能充当独立实验 n。只有汇总值时，必须实际具备适用检验
所需的均值、变异、n 和实验单位等信息；不由误差条或技术重复猜出生物学 n，不制造控制分布。
不足以重算时进入缺 p 策略，而非自动选其他检验。

BH family 固定为 **`source × 完整 condition × phenotype × assay_method × 原实验批次`**。
family 内所有药物/剂量的合法 p **先参与 BH，再筛效应幅度**；非 hits、幅度未达标者也参与，
不得先筛 `FC >= 1.2`/`FC <= 0.8` 再调整。保留 family ID、成员及数量、原始 p、调整结果、
作者/重算路径和缺失范围。不得把已报道显著点当成完整 family，或以少数可见星号补齐 p。
无法证明某调整值适用时不假装已经完成所需 BH；按可获得证据执行缺 p 策略并记录缺口。
原实验批次由来源设计确定，不按药物、剂量或显著性人为拆分。批次未标识时，只有来源证据明确
属于一套共同实验流程才可登记一个 campaign；无法建立所需 family 且无适用作者 BH 结果时，
所需 p_adj 视为不可得。
BH cutoff 是本项目的标注规则；不据此声称再经幅度筛选后的最终标签集合严格满足5%的FDR保证，
同时保留所用方法及实验依赖结构的证据。

## 5. 数值 endpoint 与准入

在身份、条件和测量有效性等基础准入通过后，按下表判定；统计显著性为 `p_adj < 0.05`。
阈值是项目初始操作参数，不声称已完成生物学校准。0 表示联合变化规则未满足，
不等于证明生物学等效；±1 是读数变化方向，不自动等于药物获益/损害。

| FC / 统计证据 | corpus 结果 | 必须保留的解释 |
| --- | --- | --- |
| `FC >= 1.2` 且 `p_adj < 0.05` | +1 | 对照、效应、检验及 family 证据 |
| `FC <= 0.8` 且 `p_adj < 0.05` | −1 | 同上 |
| 达到任一幅度边界，且已知 `p_adj >= 0.05` | 0 | 真实非显著结果，不改写为等效证据 |
| `0.8 < FC < 1.2` | 0 | 幅度已失败；p 若缺失仍标缺失，不声称统计非显著 |
| 幅度达标，但必需 p_adj 无法取得/合法计算 | 排除 corpus | missing/unevaluable significance 及来源缺口 |
| FC 未知、尺度/匹配分母不明，或分母非正/非法 | 排除 corpus | 不伪造 FC；即使存在 p 也不能代替效应证据 |
| 只有星号/ns，且无可复算信息 | 按上面缺 p 分支 | 不将 ns 写成 p=1 或擅自赋 0 |

不能用 NaN 的比较结果落入程序默认 `else -> 0`。源层允许未知、无效或被排除的记录和
null 证据字段；**最终 corpus 不允许未知/null endpoint**。排除不是删原件，逐对象保留全部原因。

## 6. 覆盖、规模与来源现实

不因单个 condition/phenotype 文件不足 500 而过滤它；小文件同样发布。
第三类的 500 规模目标继续报告，交给后续明确选择/合并条件的 model view 验收；
本轮不实现该 view，也不据 corpus 行数宣称目标达成。

覆盖报告同时给出三个口径，并写清 condition/phenotype 范围及分母：

1. 来源测量/值数：区分独立实验、技术重复、汇总和表示；raw+normalized 不算两个实验。
2. corpus 样本行数：唯一 `(condition, phenotype, source, drug, dose, endpoint)`，并按 −1/0/+1 分列。
3. 去除 source 与 endpoint 后的唯一药物–剂量组合数：在声明的条件/phenotype 范围内取并集，
   同时报告药物数，避免跨来源、相反 endpoint 或重复把规模虚增。

另报告来源有值/基础资格/FC可计算/统计可判定/入 corpus/排除的损耗；多原因不相加作互斥分区。
尚未运行的计数为 null+`not_run`，不能写零或 pass。未来 view 报告须再声明自己的合并/计数规则。

源审计显示 Sirenko beat count 有 524 个至少一份 raw 的药物标签–剂量组合、514 个双 raw
齐全组合，但尚未做完整身份/用途/零值与控制验收；附件缺 vehicle 原值/方差。
Heliste 每氧条件 2,756 个有值标签–剂量组合，但缺明确分母、重复/变异信息，且时程/药物持续
状态有缺口。**两张大表均可能因 p/control 不足，使幅度达标者被排除 corpus**。
不能把来源覆盖直接当可标注样本数，也不能保证 500 已达到；原件仍可正常归档。

## 7. 离线验收边界与下一步

| 验收对象/反例 | 预期结果 |
| --- | --- |
| 原件/本地新增/发布包 | source manifest、local 登记、发布哈希各自可复核；原件未改写 |
| 每个 corpus 文件 | condition/phenotype 一致，规范键唯一，readout 非空且逐值证据完整，endpoint 数值且仅 −1/0/+1 |
| 同源同键两次真实实验 | 一行保留全部值和原实验/比较关系；不静默平均 |
| 跨 source 或相反 endpoint | 保持不同样本；相同药物–剂量组合并集计数不随之增加 |
| raw+normalized、summary n=6 | 不增加独立重复；summary 不展开成六个孔 |
| 无 HR 与未知 HR | 前者已知时长 0，后者 null；不互换 |
| MTT/CCK8 与 CTG/Calcein | 前两者共同 phenotype 且方法保留；后两者不自动归入 |
| 技术孔冒充 n、按显著/幅度预筛 BH | 验收失败，统计需在原比较与完整规定 family 上重做 |
| 已有适用 BH p_adj | 不再次调整；保留作者方法/family 证据 |
| FC 恰为 0.8/1.2，p_adj 恰为 0.05 | 幅度达标、显著性失败，endpoint=0；缺 p 时则排除 |
| FC=1.1 且缺 p；FC=1.3 且缺 p | 前者 0 且 p 仍缺失；后者排除，不能都默认 0 |
| control=0、未知分母、技术失败零 | 不造 FC/样本；保留源事实、有效性判定及排除原因 |
| 单文件 20 行或 corpus 总行数 >500 | 小文件保留；总行数不替代后续 model view 的规模验收 |

`to-spec` 接收本设计与 ADR/术语，固定纯实现细节：Parquet 嵌套类型与关联文件、单位/枚举、
稳定 ID/文件名编码、JSON manifest/schema、统计与缺失原因字段、离线 fixture、规则版本、
可重跑入口及发布报告。实现必须同时验收“形成数值 corpus”和“证据不足但完整归档/排除”的
合法输出。设计已定，不借实现选择重开科学规则；随后按 `to-tickets` → `implement-spec` 推进。
