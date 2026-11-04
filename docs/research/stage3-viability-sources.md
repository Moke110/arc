# Stage 3：24–48 h 给药活性数据来源初步核验

核验日期：2026-11-04。对应 [核验 24–48 h 心肌细胞活性及缺氧复氧药物表型数据源](https://github.com/Moke110/arc/issues/22)，遵循 [ADR 0034](../adr/0034-stage3-phenotype-data-boundary.md)。本轮只核验来源，不下载大型数据、不生产入库、不规定三分类阈值。以下“可整理”首先指可追溯到实验和图表，**不等于已有可批量下载的训练数据库**。

## 初步结论与 research resolution 建议

建议以“**初步来源可行性盘点完成；可启动小规模人工整理试点，批量数据库可行性仍未证实**”记录本轮 resolution。常氧 48 h MTT、24 h CCK8，以及 24 h 药物预处理后缺氧/复氧再测 MTT/CCK8，均找到方法明确的原始论文（下表 A–D）；这回答了精确实验是否存在。它们主要提供论文图级汇总，不能据此宣布已有可扩展数据集，也没有完成跨论文去重。MTT/CCK8 两类没有 500 条门槛；本轮不回答其他读出是否达到 500 条。[A], [B], [C], [D]

优先试点：A 检查人 hiPSC-CM 和 AC16 的 48 h 剂量反应图能否稳定提取；B 检查同一研究内常氧与 H/R CCK8 的条件拆分；D 为动物原代层提供独立验证例。C 可补充，但正文对实验时间有歧义，需保留而非自行修正。以上是整理顺序建议，不是入库授权。

## 精确条件的单篇实验来源

“组合”按药物 × 剂量 × 细胞模型 × 暴露/损伤条件 × 读出计算；不把对照或重复孔增加为独立处理。下面是逐源核验下界，不是可合并训练样本总数。药物标识沿用作者名称，尚未核对 PubChem/InChIKey、盐型或批次。

| 来源 | 模型、药物与顺序 | 对照、读出、重复与数量口径 | 获取与缺口 |
| --- | --- | --- | --- |
| **A：Linn et al., 2024** | 人 MUi019-A hiPSC-CM；另有 AC16。doxorubicin、remdesivir、panduratin A；给药 **48 h** 后 MTT，另加显色 3 h；取 Fig. 3 的细胞毒性实验，勿混入感染后抗病毒实验。[A] | 非给药对照=100%；剂量范围 0–100 μM。**3 种药、6 条药物×模型剂量曲线**（3×2，属本笔记计算）；剂量级组合数未逐点核验。图注分别报告 CC50 的 3 次独立实验、曲线点的 3 个技术重复，不能将两者相乘。图示 mean±SD / mean±SEM。[A] | [全文/PDF][A]，Fig. 3、§2.3、§4.5；CC BY 4.0；声明数据包含于文章。未核验到逐孔机器可读附件；精确剂量网格、vehicle 对照及原始重复需补核。[A] |
| **B：Liu et al., 2021，icariin** | 大鼠 H9C2 细胞系；icariin **2.5、5、10、20 μM，24 h** 常氧 CCK8；另有 **24 h 预处理→4 h 无糖缺氧→24 h 复氧** 后 CCK8。方法的缺氧气体为95% N2/5% CO2，复氧95% air/5% CO2。[B] | 常氧可确认 **1 种药、4 个非零剂量组合**。H/R 剂量效应范围2.5–10 μM，因本轮未逐点核对 Fig. 2B，仅确认至少1个 H/R 组合，不报完整计数。含 H/R/PBS 比较；还有 erastin、Znpp 的机制联合干预，不能当作 icariin 单药重复。[B] | [PMC全文][B]，Methods、Fig. 2/5；CC BY 4.0；图级 mean±SD，图注称细胞实验重复3次，独立培养/孔级结构未完全核验。作者承诺依请求提供 raw data；未验证公开逐孔文件。注意 Fig. 5 图注用“诱导后处理”措辞，与总方法的“先给药后 H/R”有差异，联合干预暂缓纳入。[B] |
| **C：Yuan et al., 2024，NR** | H9c2；nicotinamide riboside（NR；标题写 ribose），0、0.01、0.1、1、5、10 mM。方法明确 **24 h NR预处理→3 h缺氧（1% O2）→换液→3 h复氧（95% air）**。另有无 H/R 的毒性试验。[C] | **1 种主药**，Fig. 2A/B 同一剂量系列，按方法与图注可列 **5常氧+5 H/R=10个非零剂量组合**作为待精确数值提取清单；CCK8，vehicle比较。每组实验重复3次、mean±SD；未说明孔重复与独立培养拆分。EX527机制组合不计入上述10条。[C] | [全文、PDF/XML入口][C]；CC BY；原始数据需向通讯作者请求。正文§Results另有“在H/R培养24 h”措辞，与方法3+3 h不一致；**10条是候选清单，不是已验收样本**。药品名称/盐型也需消歧。[C] |
| **D：Duan et al., 2013，SYA** | 新生大鼠原代心肌细胞；safflor yellow A（SYA）40、60、80 nmol/L、N-acetylcysteine（NAC）200 μmol/L，均 **24 h预处理→3 h无氧→6 h复氧**，Fig. 2B MTT。[D] | **2种药、4个非零单药×剂量×A/R×MTT组合**；对照为正常与A/R模型。Fig. 2A的常氧毒性实验仅4 h，不能混为24 h常氧数据。本轮未核验重复类型/数量及误差定义。[D] | [PMC全文/图入口][D]；本轮已核方法、结果及图注文本，未核验原始逐孔下载、专门数据许可及补充文件。适合补足动物原代层的人工整理试点。[D] |

本轮严格来源计数只报告各行，不求跨源总和：A 的6是曲线级，B/C/D是剂量级，不能相加。C还存在时长冲突，故任何“已可用观察总数”均未成立。对 B 的常氧4条和 D 的H/R4条也仅确认实验定义，尚未转录图值或验收误差条。

## 数据平台与规模线索：尚未找到精确端点的批量入口

**DToxS 官方入口**称所有数据免费，涉及人原代心肌样细胞及健康供者的 iPSC-CM，但当前主页没有可供本轮直接枚举的 MTT/CCK8 表型文件清单。“免费”也不替代具体文件的再使用许可核验。[E]

DToxS 关联的 van Hasselt et al. 2020 论文确有 **48 h** 剂量反应；两个成人来源心肌样细胞模型、8个非零剂量（5 nM至100 μM）、vehicle、4复孔。实际读出为 **Hoechst核计数、MitoTracker Red**，不是MTT/CCK8；不能作为本票的精确命中。论文 Fig. 1 与 Supplementary Fig. 1、Source Data 是值得后续追的表型入口，GEO/RNAseq记录不能代替该表型。未核验具体表型药物数、文件行数或原始孔级粒度；**不能把全文26种KI的研究药物数乘以8剂量/2模型推算本端点规模**。[F]

因此，DToxS 是“平台级邻近线索”，A–D 是“精确端点的小论文来源”。本轮没有证实常氧或H/R的 MTT/CCK8 可从一个官方库批量导出。此为有限检索结果，不是数据库不存在的断言。

## 排除与待核验线索

| 线索 | 判定 |
| --- | --- |
| Stachydrine/SIRT1-Nrf2，2023 | H9c2常氧 **24 h**，0、10、20、50、100、200 μM，CCK8；可追5个非零常氧剂量。其H/R实验仅 **2 h预处理**、4 h缺氧/12 h复氧，排除出精确24–48 h预处理H/R；至少3次独立实验、mean±SD。尚未追逐孔文件。[G] |
| Acacetin/autophagy，2021 | 文中有12、24、36 h常氧CCK8时间探索，但后续H/R采用12 h预处理；暂作时间过滤反例，未计算剂量×时间笛卡尔积，未作为已验收数据。[H] |
| Quercetin/SIRT1/TMBIM6，2021 | 原文称“human cardiomyocytes”、24 h给药、5剂量、3 h缺氧/2 h复氧、CCK8，但未明确原代/hiPSC/细胞系身份；复氧写95% O2而非air。暂列待核验，不把它当作已确认的人原代层。[I] |

## 搜索范围、检索词与局限

本轮用 web 搜索追原始论文、PMC、出版社PDF/图注及 DToxS/LINCS 官方页；跨 human iPSC、human primary、H9c2/AC16、neonatal rat primary 检索，不限定已有转录组。实际检索词包括：

- `DToxS cardiomyocyte MTT 48 hours viability data`；`"DToxS" cardiomyocyte viability`。
- `cardiomyocytes hypoxia reoxygenation "24 h" "CCK-8" "pretreatment" drug`。
- `"cardiomyocytes" "MTT" "48 h" "screening" compounds`；`"cardiomyocytes" "CCK-8" "24 h" "iPSC" drug`。
- `"cardiomyocytes" "MTT" "24 h" "neonatal" "pretreatment" "reoxygenation"`。
- 候选全名 + `supplementary` / `Data Availability Statement` / `concentrations`；DToxS关联论文全名。

这是限时初查，不是系统综述：未系统遍历PubMed/Europe PMC全部命中，未联系作者、未做化学结构标准化、未数字化图、未审计所有勘误/撤稿与许可。PMC部分 open 请求触发验证码，Nature HTML曾跳认证、MDPI HTML曾429；本轮以检索到的原文文本、可访问出版社PDF与图注交叉核对，未将访问失败写成缺失数据。人原代精确MTT/CCK8仍未确认；DToxS原代样模型只属邻近端点。

下一轮最有信息量的工作：核对 A 的逐剂量坐标/数据附件，向 B/C 的可请求原始数据入口确认可获取性（联系作者需另行授权），解决 C 的时长冲突，再决定“图级汇总是否足以入库”与可比较条件。独立培养、技术孔、供者层级应分字段保存；不显著不能直接赋0类，且H/R保护应对照H/R vehicle，不能默认对照常氧。[ADR 0034](../adr/0034-stage3-phenotype-data-boundary.md)

## 主来源索引

- [A：Linn et al., 2024，出版社PDF](https://www.mdpi.com/1422-0067/25/3/1427/pdf?version=1706090030)，doi:10.3390/ijms25031427，Fig. 3、§4.5、Data Availability；[PMC副本](https://pmc.ncbi.nlm.nih.gov/articles/PMC10855687/)。
- [B：Liu et al., 2021，PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8564343/)，doi:10.1002/2211-5463.13276，Methods、Fig. 2/5、Data Availability；[出版社](https://febs.onlinelibrary.wiley.com/doi/10.1002/2211-5463.13276)。
- [C：Yuan et al., 2024，出版社全文](https://www.spandidos-publications.com/10.3892/etm.2024.12475)，doi:10.3892/etm.2024.12475，Cell culture and H/R treatments、Fig. 2、Statistical analysis、Availability of data。
- [D：Duan et al., 2013，PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4002788/)，题名 *Safflor yellow A protects neonatal rat cardiomyocytes against anoxia/reoxygenation injury in vitro*，Experimental groups、Fig. 2。
- [E：DToxS 实验室官方数据主页](https://iyengarlab.org/dtoxs/)；[NIH LINCS 官方中心介绍](https://lincsproject.org/LINCS/centers/data-and-signature-generating-centers/dtoxs)。
- [F：van Hasselt et al., 2020，出版社PDF](https://www.nature.com/articles/s41467-020-18396-7.pdf)，doi:10.1038/s41467-020-18396-7，Dose–response experiments、Fig. 1、Data availability。
- [G：Zhu et al., 2023，出版社全文](https://link.springer.com/article/10.1186/s13019-023-02363-6)，doi:10.1186/s13019-023-02363-6，Fig. 1。
- [H：Liu et al., 2021，PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8263275/)，doi:10.1155/2021/9979843，*Acacetin Protects Myocardial Cells against Hypoxia-Reoxygenation Injury through Activation of Autophagy*，Methods、Fig. 1。
- [I：Chang et al., 2021，PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8024107/)，doi:10.1155/2021/5529913，§2.1–2.2、Fig. 1。

[A]: https://www.mdpi.com/1422-0067/25/3/1427/pdf?version=1706090030
[B]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8564343/
[C]: https://www.spandidos-publications.com/10.3892/etm.2024.12475
[D]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4002788/
[E]: https://iyengarlab.org/dtoxs/
[F]: https://www.nature.com/articles/s41467-020-18396-7.pdf
[G]: https://link.springer.com/article/10.1186/s13019-023-02363-6
[H]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8263275/
[I]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8024107/
