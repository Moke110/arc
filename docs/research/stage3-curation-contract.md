# Stage 3：首批表型整理溯源契约与 spec 交接

> 历史 wayfinder 决策资产。后续用户访谈已形成[正式 corpus 设计](stage3-formal-corpus-design.md)；目录、样本键、phenotype 合并、数值标签、格式和规模门槛等冲突部分以正式设计及 ADR 0038–0041 为准。本文保留当时的取舍与证据。

日期：2026-11-04。对应[确定首批表型数据库的整理溯源契约与 spec 交接条件](https://github.com/Moke110/arc/issues/27)，
上游为[确定 Stage 3 首批表型数据的纳入分层与优先来源](https://github.com/Moke110/arc/issues/24)、
[定义独立药物表型样本与跨研究条件等价及去重规则](https://github.com/Moke110/arc/issues/25)、
[确定表型变化方向的对照及三分类标签证据要求](https://github.com/Moke110/arc/issues/26)。
本契约整合其[纳入决定](stage3-inclusion-decision.md)、[计数决定](stage3-counting-decision.md)、
[方向证据决定](stage3-direction-evidence-decision.md)及[原件审计](stage3-sample-count-audit.md)，不替代原始证据。
约束依据：[ADR 0030](../adr/0030-external-data-provenance-and-step-layout.md)、
[ADR 0034](../adr/0034-stage3-phenotype-data-boundary.md)、
[ADR 0035](../adr/0035-stage3-stratified-phenotype-curation.md)、
[ADR 0036](../adr/0036-stage3-continuous-evidence-before-direction-labels.md)与 [GLOSSARY](../../GLOSSARY.md)。

本契约的载体与布局决定已落档为[ADR 0037](../adr/0037-stage3-phenotype-release-layout.md)。

这是 AFK 代理依据用户授权作出的设计决定，使用 grilling 检查边界、domain-modeling 检查实体关系；
没有虚构用户答复。**本轮规划完成，不等于数据、500 门槛或标签验收完成。**
本文件只交接设计，不发布 spec、不创建生产目录、不下载/导入数据、不运行 Stage 2 或训练。

## 1. 首版范围和两处地图 fog 的结论

首版按 Sirenko 2013 `24hr_beat_count` 优先整理第三类主目标；MTT/CCK8 按
A/Linn → B/Liu icariin → D/Duan SYA/NAC 做图级试点。所有目标图中的合格及不合格剂量均保留审计，
不只取 hits、显著点或最佳剂量。Linn 的 hiPSC-CM/AC16、Liu 的常氧/HR 各自分层。
C/NR 保留时长冲突证据并隔离；D 的 4 h 常氧不计目标。Heliste S1 保留完整剂量系列的扩展清单，
持续给药与归一化分母未明时保持 HR 隔离，常氧约 96 h 不冒充直接 24–48 h。
MTT、CCK8、CTG、beat count 及峰形参数各为独立 readout，不拼成通用 viability。

地图“是否还需许可确认或请求作者提供小论文原始数据”的首版决定为：**不以联系作者为首版前置条件**。
用当前可合法获得的公开材料建立带缺口的 inventory/pilot；未公开 raw、误差或不清楚许可均记录并隔离。
许可事实核对属于逐资产整理，未知许可不等于开放许可，也不授权再分发原图/附件。
数据库不向外打包许可未明确的资产；可发布本地可复核的来源索引、缺口和允许发布的自制整理结果，
不能发布的结果只列元数据及阻断原因。首版不主动请求作者、发邮件或替用户接受额外条款；实际联系须另行明确授权。

地图“三类数据都不足时是否放宽范围”的决定为：**本轮不放宽任何时间、模型、指标或 ≥500 门槛**。
前两类无500要求，小规模或空结果照实交接；第三类若 qualification 后不足500，输出
`insufficient_qualified_combinations`，保留 exploratory 数据，不无限找替代指标凑数。
三类均不足则交接 `insufficient_scope`，附各类损耗及证据缺口。改变范围另起一次决策努力，
不是本地图继续悬置的问题，也不授权自动扩展。数值标签校准归后续独立模型/labeling spec，当前默认 null。

## 2. 载体与未来布局的冻结决定

首版交换格式固定为 **UTF-8 CSV 表 + JSON schema/规则/发布记录**，不引入数据库服务、h5ad 或
模型输入格式。CSV 方言采用 RFC 4180 的字段引用规则、LF 记录分隔；首行为字段名，数值用十进制点，
时间用 UTC ISO 8601；schema 为每列声明类型、
可空性、枚举、单位及外键。表按主键排序；数值不写 NaN/Inf。JSON 用原生 null。
CSV 中 `\N` 专指 null，空字符串仍是字符串；文本列中的真实反斜杠先转义为双反斜杠，
读入先识别整格 `\N` 再还原文本转义。源单元格的空白、`N/A`、公式文本等另存 `source_token`，
不可用 CSV null 规则改写来源语义。列表/多证据使用关联表，不塞逗号数组。

未来外部资产放 **`pert_data/<item_id>/`**；`item_id` 按上游研究产品命名，候选名称为
`Sirenko-2013-CM-Phenotype`、`Heliste-2021-CM-CTG`、`Linn-CM-MTT`（须在源 inventory 核对后登记，示例不是已建项）。
它们是药物扰动实验表型来源，属于既有 `pert_data` 类别；[ADR 0001](../adr/0001-category-taxonomy.md)
的类别容纳新扰动数据，不要求数据必须是转录组。一次下载产品一个项，不能按端点/模型复制原件成多个项；
同一产品多个端点由下游分层。无需第五外部类别，无需挪动 Stage 2 现有项或让现有管线自动消费它们。

每项严格采用 ADR 0030 的 `script/`、`source/`、可选 `local/`、`author_log.json`、`manifest.sha256`：

- `source/` 保留上游原始字节、文件名和目录；不改 XLS 单元格、不重命名或裁图覆盖原件。历史不同版下载各保留上游容器/版本来源，禁止覆盖旧证据。
- 下载/离线校验放项的 `script/`。作者日志记录 URL、上游版本、retrieved_at、媒体类型、许可及证据、文件数/字节/树摘要和验证结果。
- `manifest.sha256` 只列 `source/` 相对路径，按路径排序，树摘要为 manifest 文件字节 SHA256；不把 `local/` 混入 source manifest。
- 裁图、手工数字化 CSV、注册信息快照等本地产物放 `local/`，逐路径登记 `local_additions[]` 的来源、规则/软件版本和 SHA256。生产标准化跨源表归 Stage 3 消费步骤。

未来 Stage 3 使用 `stage3/AGENTS.md`、`stage3/README.md` 路由，首个步骤固定为
`stage3/01_phenotype_database/`，采用 `AGENTS.md`、`README.md`、`execution/`、`output/` 四成员结构；
直接供人查看的 `coverage.csv`、`attrition.csv`、`gate_report.csv` 放步骤根部。
`execution/` 放提取/映射/验证程序、配置、测试和运行日志；`output/releases/<release_id>/` 放
`tables/*.csv`、`schema.json`、`rules.json`、`release.json`、`checksums.sha256` 及机器门禁报告。
此处 output 的 checksums 覆盖发布产物，**不是**外部项的 source manifest。
发布不可变；修正生成新 release，并连到旧 release。Stage 3 只引用项相对路径及文件哈希，不复制外部原件进 output。
本轮不实际创建这棵树；布局已由 ADR 0037 记录，后续 grill-with-docs 校核领域文档一致性，之后 spec 才实施。

## 3. 关系表、粒度与必要字段

以下表属于 schema v1 的最小实体集，所有 `_id` 为发布内稳定、不含语义猜测的字符串主键。
每次修正规范键须保存旧新键映射；同一 release 内主键唯一、外键不得悬空。
规范键的生成规则/字段有序序列记录在 rules 中，不以名称或来源行号充当母体身份。
源身份未解析时有局部 formulation ID，但 `drug_key=null`；不伪造 qualified 键。

| 表／一行粒度 | 主键与关联 | 必备信息（未知可空但须有状态/理由） |
| --- | --- | --- |
| `studies`：一个研究/实验 campaign | `study_id`；同研究多论文用证据关联 | 标题、DOI/PMID/项目号、campaign标识、引用、数据重用/重复出版关系 |
| `assets`：一个不可变文件版本 | `asset_id` → study | 外部项、source/local相对路径、URL、上游版、SHA256、bytes、媒体类型、retrieved_at、license及证据、supersedes_asset_id；镜像同字节仍记URL而不创造实验 |
| `evidence`：一个可定位证据片段 | `evidence_id` → asset | XLS sheet、1-based行列及A1地址/表头；PDF印刷页与文件页、图号/面板/曲线或柱/剂量/图注、bounding box；文本章节/段落、短摘录；定位格式版本 |
| `entity_evidence`：一次字段证据链接 | `(entity_table, entity_id, field_name, evidence_id)` | relation=`supports/contradicts/derived_from`、提取人/代理身份、时间、方法；同一字段允许多个冲突证据 |
| `drugs`：一个经解析活性母体 | `drug_key` | 母体结构/注册标识、已知立体组成、药物用途资格及证据、identity_mapping_rule_version；不把前药/代谢物当别名 |
| `formulations`：一种来源给药形式 | `formulation_id` → study；drug_key可空 | 原名、CAS、厂商/catalog/批号、原结构、盐/水合物/组分比例、实际分子量与浓度计量基准、身份状态；注册查询使用冻结快照及asset引用 |
| `readouts`：一个测量量与检测体系 | `readout_key` | 作者端点名、技术/试剂、量纲/单位、采样窗口、检测范围；MTT/CCK8/CTG/beat count分键 |
| `strata`：一个可比较目标 | `stratum_key` → readout | 纳入类别、物种、模型类型/细胞系或产品、成熟/分化/背景、培养体系、HR/氧方案、持续/洗脱/补药策略、固定背景、vehicle profile、参照/归一化语义、允许剂量/时间轴、跨源等价结论 |
| `conditions`：层内一个精确自变量组合 | `condition_key` → stratum | 原剂量/单位、母体摩尔浓度、给药日程、名义时点及首次给药/检测开始/试剂孵育/读数/HR起止、实际vehicle浓度；逐时钟出处/计算规则，未知不推算 |
| `combinations`：唯一药物×精确条件×指标 | `combination_id` → drug、condition、readout | 唯一约束 `(drug_key,condition_key,readout_key)`；与condition所属层readout一致；qualified/excluded/unresolved、理由及决策版本；不存一个覆盖所有重复的伪均值 |
| `observations`：一个来源实验单位的一项读数 | `observation_id` → study、formulation、condition、readout；combination可空 | role=`treated/vehicle/normal_reference/positive_control/blank`、campaign/供者/分化批/培养批/板/孔、replicate_key与作者重复号、类型technical/biological/independent_study/summary/unknown、汇总n及实验单位、measurement_validity |
| `values`：观察的一种表示/统计量 | `value_id` → observation | representation=`raw/author_normalized/digitized/derived`、statistic=`individual/mean/median/SD/SEM/CI_lower/CI_upper`、有限value或null、unit、source_token、缺失/零/删失状态、归一化分母/公式、extraction_id/derivation_id（适用时） |
| `comparisons`：一项处理与明确参考的对比 | `comparison_id` → treated observation或组合、readout | reference语义、匹配依据、control_validity、contrast_validity、uncertainty_status、作者陈述/p/test及比较对象、连续效应/尺度/公式、描述符号、label=null/−1/0/+1、label_status/reasons、label_rule_version |
| `comparison_members`：一个对比的输入值 | `(comparison_id, value_id, role)` | role=`treated/control/uncertainty`；共享分母、配对关系和权重；作者仅声明归一化基准时无control原值成员，保留reference evidence，不造一个100原始孔 |
| `decisions`：一个对象的一次审核 | `decision_id` →对象表/ID | dimension、status、reason_code、证据、reviewer、时间、rule_version、supersedes_decision_id；保留历史不覆盖排除理由 |
| `conflicts`：一个待核实或已解决冲突 | `conflict_id` →对象及相关evidence/decisions | 类型identity/time/control/value/duplicate/direction/version、影响范围、open/resolved、解决理由/证据/替代版；不按多数票消冲突 |
| `extractions`：一次图值读取及复核 | `extraction_id` → evidence | 软件/版本、线性/对数轴、坐标校准点、曲线/柱与剂量映射、读取值/界限/精度、误差条语义和n、提取人/复核人、差异及解决、校准/裁图local资产哈希 |
| `derivations`：一次派生操作 | `derivation_id` | 输入value/evidence引用（另用关联表）、公式、参数、代码commit、rule_version、输出对象、时间；转换与标签规则版本分开 |
| `derivation_inputs`：一次派生的一项输入 | `(derivation_id, input_index)` → derivation；value_id或evidence_id外键恰一个非null | input_index为派生内稳定序号；role=`operand/reference/uncertainty/rule_evidence`，分别表示运算数、参照输入、不确定性输入、支持规则/参数的证据；多输入各占一行，不用逗号数组 |

原文只给柱均值且n=6时是一条summary observation及多个统计量，不展开六条孔。
comparisons的target_observation_id与target_combination_id必须恰有一个非null；
组合层对比只有存在明确派生规则时才建立，不能自动将跨研究观察合成一个contrast。
raw与normalized属于同一 observation 的不同 value；未公开孔号使用可追溯来源局部replicate_key，
不冒充真实孔。对照也是 observation，有实际条件和role，不贡献药物组合；Sirenko `Control=1`
仍是作者化合物标志，另存原值，不能自动转为vehicle role。
vehicle、normal_reference、blank等无目标药物的观察允许formulation_id为null，combination_id必须为null；
不得为连接处理组而伪造目标drug或formulation。positive_control若确有参比药物可关联其真实formulation，
仍不自动建立目标药物组合。
不满足身份资格但可保留的记录仍进 observations/values，combination关联可空，排除原因不丢。
来源缺关键条件时可建立局部未解析 condition/stratum，但不得靠 `unknown==unknown` 归并或贡献qualified计数。

## 4. 单位、状态、提取与版本约束

浓度规范为 `mol/L`，时间为 `h`，vehicle为明确 `v/v` 或 `w/v` 的 fraction；百分比保留原值并标明分母语义。
质量浓度只有实际给药形式、分子量与计量基准都确定才转换；单位精确换写不制造新剂量。
原单位/数值/有效位数与转换公式始终并存。beat count保留 count及测量窗口；未知窗口不换bpm。
吸光度、RLU和作者相对百分比各自标单位，不能互相替换。SD/SEM/CI必须连到均值、n、独立实验单位和CI水平；
不知n不做SEM→SD，不把数字化误差当生物方差。

状态采用正交维度，而非一个“好/坏”总分：

- identity、用途资格、时间/协议、条件等价、计数资格：`qualified/excluded/unresolved`，逐维保留理由。
- measurement/control/contrast：`valid/invalid/unresolved/not_applicable`；uncertainty：`sufficient/insufficient/unknown/not_applicable`。
- 缺失理由：`source_blank/not_reported/not_readable/not_applicable/parse_failure`；保留source_token。
- 零语义：`observed_zero/source_derived_zero/no_peak_undefined/below_detection/zero_unresolved/not_zero`；技术失败另有validity。
- label_status：`unresolved/label_ready/excluded`。null必须有reason，例如missing_control、unknown_denominator、insufficient_uncertainty、conflicting_evidence、rule_not_frozen。

状态间允许同时有多个缺口。来源证据以 `primary_numeric/author_summary/figure_digitized/author_statement/curator_derived`
标方法和来源层级，不用一个证据等级替代控制、技术和不确定性验收。图值不是天生不合格，raw也不是天生合格。
任何标签非null必须有冻结labeling规则、有效匹配对照及充分不确定性；0需正面等效证据，ns不等于0。
组合计数资格不强求已公开vehicle方差或已定label阈值，但必须符合上游计数/有效测量证据要求；
控制信息不足应准确限定其影响，不能反过来把所有qualified都说成label-ready。

人工图提取先冻结原图资产哈希与定位、轴/组别/剂量映射和可读范围，再提取中心值与误差界限；
无法辨认则保留not_readable。每个待发布图值由另一复核者独立读取校准与值（可为独立代理审核），
不把同次输出自检称独立复核。差异超过双方声明的提取区间时进入conflict，不能静默平均；
重叠也只给合理有效位数和提取区间，不伪装逐孔raw。误差类型/n缺失保持unknown。
不可复核提取可留inventory，不能通过pilot值验收；符号可能受读取误差影响时不赋标签。

文件版本以哈希而非basename识别；author manuscript与publisher版分别登记，不凭同标题声称字节一致。
勘误建立supersedes链接，旧版及其旧decision保留；source manifest不因清洗变化。
一份来源单元格只能出现一次规范提取（保留修订历史），重复出版/镜像链接到同实验，不能新增replicate。
不同研究同组合结果冲突保留全部有效观察；不平均抹平、不挑有利研究，组合标签默认null，
后续labeling spec才决定层级及冲突策略。

## 5. 离线验收与失败出口

未来 verifier 读冻结文件和规则，不请求网络。每次 release 出具逐gate `pass/fail/not_run/not_applicable`、
命令/代码commit、输入哈希、规则版本、发现数及失败对象定位；`not_run`绝不能转pass。
源文件完整性与整体语义门禁分开；一次 `verify_all.sh` 成功不能证明500或标签正确。
G0–G4的pass表示相应验证、判定和隔离规则已正确执行，不要求保留的每条来源记录都qualified。
有明确unresolved状态且正确隔离的记录不会自动使整个pilot结构验收失败；G5仅统计通过资格的集合。

| Gate | 必须检查的内容 | 失败时允许发布什么 |
| --- | --- | --- |
| G0 provenance/rights | source manifest离线校验、asset版本/URL/许可状态、local逐项登记、全部数值和判定回链证据 | 缺文件/定位则相关条目只进缺口清单；许可未知或受限则不对外分发相应资产/受限制结果 |
| G1 structural | schema类型/枚举/null、主键唯一、外键及值到证据关系、readout/层一致、旧新版本链无循环 | fail禁止发布标准化pilot；inventory可列未解析源信息并明确失败 |
| G2 fidelity | 源单元格/图定位、源token与原值一致、raw/normalized同观察不重复；图复核与误差；单位转换/派生可重算 | 错误记录隔离；不能以0填缺失、以作者100伪造对照或按n扩孔 |
| G3 identity/conditions | 用途/母体/形式、摩尔基准、精确时间与时钟、模型/HR/持续给药、vehicle profile与对照语义；逐字段跨源等价 | 未知留unresolved层，不进入qualified并集；不同层不相加 |
| G4 dedup/counts | 重复出版/别名/盐映射、同组合多来源取并集、重复类型；重算各层药物/组合/重复及损耗 | 映射未执行则qualified数为null而非0或raw覆盖数；不声称500 |
| G5 third-category ≥500 | G0–G4通过后，一个固定readout、可比较stratum中的qualified唯一组合数≥500；报告集合与每条资格证据 | 未执行=`not_run`；实际<500=`fail:insufficient_qualified_combinations`，仅exploratory。前两类此gate为not_applicable |
| G6 contrasts/labels | 对照生物/vehicle/时钟匹配、量纲/分母、零与删失、重复/误差、规则冻结、非null标签证据及冲突隔离 | label全null可通过“未伪造标签”的数据库一致性检查；label-ready/training验收仍not_run或fail，不能称training-ready |

inventory可发布条件：来源、范围、资产可得性/许可、已知缺口及状态可复核，G0对已取得资产通过；
未取得资产明确unavailable且不宣称已验真。inventory不是一份有值/qualified数据库。
pilot可发布条件：选定的实际提取记录G0–G2通过，G3/G4对记录逐项判定并保留隔离/排除，
结果表清楚区分有值、qualified与unresolved；**不要求G5通过**，但不得把pilot叫第三类达标数据集。
若图值不可读，发布“试点失败及定位/原因”而非假数据。

coverage/attrition必须逐类别、stratum、readout和来源列出：来源标签数、有值组合、解析母体数、
qualified药物数、qualified唯一组合数、technical/biological/independent_study/summary/unknown重复数、
contrast可计算数、uncertainty充分数、label-ready数及−1/0/+1/null数。每列注明计数实体/分母；
跨源合计qualified用并集而非来源数求和。损耗表按对象列出全部原因，同时提供不重复的最终状态分区，
多原因统计不得误作互斥分区相加。未知计数为null+not_run，不把尚未验收写成零。

至少以下边界成为后续验收夹具（均为规则预期，当前未运行）：

1. Sirenko raw空白/normalized=0保持两种事实；无该重复观测证据，不凭派生零增加有值组合。
2. Sirenko双raw零先标zero_unresolved；未证明真实停搏/有效测量时不贡献qualified，绝不自动赋0。
3. 两个孔、多个供者、独立研究、raw+normalized、盐别名均不直接增加同层同药同条件组合；同药四合格剂量最多四组合、一药物。
4. 名义24 h与23.5 h+30 min时钟并存；Heliste约96 h总历时不改写为已知24 h总用药，未知持续给药不进入已知协议层。
5. HR无药对照与常氧对照分开；归一化分母未知不得造ratio；分母0不加pseudocount。
6. 500有值组合但499 qualified时G5 fail；524覆盖/514双raw齐全不证明Sirenko过门槛；Heliste每氧条件2,756有值亦不证明qualified。
7. 全部label=null且有reason时可有合法pilot；ns/恰100/缺失均不赋0，冲突保留；rule未冻结时label-ready不得为真。

## 6. 交接到固定工作流

本轮可关闭[Stage 3 决策地图：心肌细胞药物表型数据库](https://github.com/Moke110/arc/issues/21)
的规划出口。事实未知仍为数据状态，不是未决定如何处理的设计分支：药物映射后Sirenko是否≥500、
零值技术含义、逐图可读性、未公开控制/误差、Heliste分母/持续给药及逐资产许可仍需实际整理核验。
处理方法已定为证据不足隔离、门禁失败/未运行，不用联系作者或放宽科学边界消除未知。

| 下一步 | 必须输入 | 输出及边界 |
| --- | --- | --- |
| `grill-with-docs` | ADR0034–0037、GLOSSARY、三份上游决定、本契约、来源研究/原件审计、地图结题索引 | 校核契约与已落档布局/载体及领域术语的一致性；以已授权决定为基线，不伪造数据验收；若实际新冲突再记录决策 |
| `to-spec` | 上述已落档设计、精确表字段/枚举/单位与PK/FK、证据定位、版本/规则、离线gates和失败例 | `.scratch/<feature-slug>/spec.md`并发布spec issue；范围为inventory→图级pilot/Sirenko资格审计与可复核release，允许G5失败且标签全null；不夹带训练、正式标签或作者联系 |
| `to-tickets` | 发布spec、验收矩阵/源资产边界 | 按依赖发布tracer-bullet票及native blocked_by：schema/一个真实记录贯通→溯源提取与校验→身份/条件计数→来源扩展与发布审核；每票可离线验收，后续implement-spec才执行 |

后续spec必须同时给“成功数据集”和“证据不足但审计完整”两种输出验收，不能以任一分支失败为由改500。
数据库整理spec不校准最终方向阈值；后续模型/labeling spec须独立冻结效应尺度、最小有意义变化/等效界值、
区间算法、独立实验单位、多重比较、冲突策略及独立校准依据，才可能生成正式标签。
Stage 2预测转录组生成/匹配/评估、扰动组输入表示、模型训练与benchmark均不在本交接内。

## Resolved checklist 与 resolution 建议

- [x] 首批来源/范围沿用纳入决定；来源保留、qualified计数、label-ready三者分开。
- [x] 组合与replicate/value层分开；药物母体/形式、时间/单位、主外键和来源定位已定。
- [x] 缺失与各种零、人工图复核、原件版本/勘误、冲突和派生规则版本有可执行契约。
- [x] 既有四类外部来源兼容；source不动、local逐项声明、manifest只source；Stage 3未来布局/格式已定但未建树。
- [x] inventory/pilot/500/标签分别验收；未运行不通过，标签可全null而不宣称训练可用。
- [x] 首版无需作者联系；未知许可隔离；不足不降门槛/扩范围，两处地图fog已有明确结论。
- [x] 设计到spec/tickets的输入输出明确；未发布spec、未实施，数值标签校准明确归后续独立范围。

建议关闭本决策票：**已决定首批整理溯源契约及后续spec验收边界**。建议地图以“规划路线完成”关闭，
引用四份决定及原件审计，注明“生产入库未执行、qualified ≥500 未验收、正式标签未验收”。
不得把本次结题改写为“已获得500样本数据库”或“已就绪三分类训练集”。
