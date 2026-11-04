# Stage 3：原始表型附表的样本数核验

核验日期：2026-11-04（Asia/Shanghai）。对应任务
[取得并逐行核验 500 样本候选的原始表型附表](https://github.com/Moke110/arc/issues/28)，
上游为[其他读数来源研究](stage3-other-readout-sources.md)，范围遵循
[ADR 0034](../adr/0034-stage3-phenotype-data-boundary.md)。

## 结论

**两个候选的真实逐剂量附表均已取得，数值覆盖已实测；此前的“只有设计容量、附件未取得”障碍已解除。**
这完成了本票的获取和计数任务，不代表已完成药物资格、化学身份归并、实验条件或三分类标签的验收。

| 来源与单一指标 | 文件内身份口径 | 实测组合数，不计重复为额外样本 | 尚未解决的纳入边界 |
| --- | --- | --- | --- |
| Sirenko 2013，24 h iPSC-CM beat count | 131 个不同 CAS 字符串，131 个不同名称 | 4 剂量；524 组合至少一份 raw 值，514 组合两份 raw 齐全 | 库含非治疗化合物；盐型/母体归并未做；零与技术失败需区分；vehicle 原始值/方差未附 |
| Heliste 2021，HL-1 CTG relative viability | 689 个不同药物标签，尚无结构身份键 | 每个固定剂量、每个氧条件均 689 个有限数值；HR 2,756，常氧 2,756 | 24 h 预处理后 HR，总历时约 96 h；持续给药状态、归一化分母及重复误差未明确 |

表内计数来自下列原始 XLS/XLSX，而非简单将论文库大小乘以剂量数。
HR 与常氧不相加验收一个条件；七个心肌指标不相加验收一个指标。
药物标签、CAS 条目、化学母体、合格治疗药物也不视为同一个计数口径。

来源：[Sirenko Table 3 原件][sirenko-table]、[Heliste Dataset S1 原件][heliste-s1]。
原文分别为 [Sirenko 2013][sirenko-paper]、[Heliste 2021][heliste-paper]。

## 获取路径、版本和完整性

旧 PMC OA API 已退役，404 不能作为论文无附件的证据；官方说明见
[PMC Cloud 文档](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/)。Heliste 的正式附件位于
[新版 public bucket](https://pmc-oa-opendata.s3.amazonaws.com/README.txt)。
Sirenko 的公开 author-manuscript 前缀只列 XML/TXT/JSON，没有附件；改由 Elsevier 正式附件 CDN 取得 XLS。
没有绕过登录或验证码。此前 PMC 的“Preparing to download”HTML 没有被当作 Excel 解析。

| 文件与用途 | bytes | SHA-256 |
| --- | ---: | --- |
| [Sirenko publisher mmc3，24 h Table 3][sirenko-table] | 231424 | `103b1601944ace5ebe11c269884a3634ab23b17ef2d7446131051a6e5c77b496` |
| [Sirenko mmc4，AUC 汇总][sirenko-auc] | 27136 | `dab37016296b2254163f1a07d12e58c026ad617132a77fac212e77f3208084f6` |
| [Sirenko mmc5，BMC/ToxPi 汇总][sirenko-bmc] | 52736 | `7246df19653c5737d2a9b56ec144fbf4618d23ffa04b45883bafc5c652e4ad09` |
| [Heliste Information 2 / Dataset S1，药物逐剂量值][heliste-s1] | 101938 | `e3bce1a6b017d1d71a2033e89986b462a294fdf576bfeff16f9b75fe60c73aac` |
| [Heliste Information 3 / Dataset S2，shRNA 数据][heliste-s2] | 2374209 | `7b26d85da51e2b3832d66921653c1c866adf772a37d34aa2aeb33909637e692e` |
| [Heliste Information 4 / Dataset S3，hits 名称][heliste-s3] | 10205 | `181119aaec3d8221317fc150a34ae62d8b436a60d4ddc547e33a0f1387beb2e3` |

Heliste 使用 PMC 版本 `PMC8371130.1`，本地附件 MD5 与[官方 JSON][heliste-json]一致。
论文为 CC BY 4.0，第三方材料例外按原文保留。Sirenko publisher 原件是 OLE/BIFF XLS；其内容与正文 Table 3 的
24 h 描述吻合，但未宣称与拿不到的 NIHMS 附件逐字节同版。PMC author-manuscript
[元数据][sirenko-json]标为 TDM / 非 Open Access，原文保留 Elsevier 版权；没有将可下载推断为可任意再分发。

这里只保留来源、哈希和核验方法；原件及临时展开表存于 `/tmp/arc-stage3-sirenko-audit/`、
`/tmp/arc-stage3-heliste-audit/`，不作为持久库或后续训练入口。正式入库仍需来源布局与溯源设计。

## Sirenko：24 h 逐剂量表

### 实表结构与身份

[Table 3][sirenko-table]有七张表，每张 132 行（含表头）×22 列，131 个化合物，身份顺序一致。
A:E 为 `Plate_Location`、`Catalog_number`、`CAS_number`、`Chemical_Name`、`Control`；
F:M 为 0.3、3、10、30 µM 各两个重复的 raw 值；N 空白；O:V 为相同剂量/重复对应的 normalized 值。
原始 beat count 列未标成 beats/min，核验保留原名，没有换算成 bpm。

CAS、名称和 plate location 各有 131 个唯一值，无表内精确重复。CAS 字符串唯一不证明化学母体唯一。
Catalog 仅 23 个实际编号，其余 108 行为 `N/A`，不能用此列去重。
名称含 sucrose、xylitol、citric acid 等，不能未经资格整理就把 131 个条目称为 131 种合格治疗药物。
完整 Supplemental Table 1 未取得；获取尝试有界停止，现有 Table 3 的身份列足以完成本次表内计数。

`Control` 列有 107 个 0、24 个 1；这与论文的 107 心毒/24 非心毒库分类相对应，
“1 为非心毒类”仍是对应推断，未取得 Table 1 图例。**这些是化合物行，不能当作 DMSO 孔或变化方向标签。**

### 每个指标单独计数

数值零按“有记录”计入，空白不计；没有自动把零判成停搏、无变化或技术失败。

| sheet | 至少一份 raw | 两份 raw 齐全 | 至少一份 normalized | 两份 normalized 齐全 | raw 空白单元格 | 两份 raw 都零 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `24hr_beat_count` | 524 | 514 | 524 | 516 | 10 | 71 |
| `24hr_peak_amplitude` | 522 | 522 | 522 | 522 | 4 | 71 |
| `24hr_peak_width` | 524 | 523 | 524 | 523 | 1 | 70 |
| `24hr_peak_decay` | 522 | 522 | 522 | 522 | 4 | 71 |
| `24hr_spacing` | 524 | 524 | 523 | 523 | 0 | 76 |
| `24hr_peak_rise` | 522 | 522 | 522 | 522 | 4 | 71 |
| `24hr_viable_cell_number` | 523 | 523 | 523 | 523 | 2 | 2 |

最后一项为 Calcein AM 阳性细胞成像计数，不是 MTT/CCK8；它是另一个可供后续选择的单指标。
当没有峰时，峰形参数的零可能只是计算回退，不能按“数值存在”保证生理有效。

### 缺失、零与对照

Beat count 共 1,048 个设计重复单元格，其中 1,038 个为 raw 数字。十个组合只保留一份原始重复。
Ascorbic acid、Biotin 在 0.3 µM 的 r1，原始 F120/F121 空白，但 normalized O120/O121 被填为零；
这两处不能作为观测到的停搏。其余八处 raw 空白在 normalized 中仍空白。

71 个组合两份 raw 都为零；453 个组合至少有一个正 raw 值。453 不是自动的可用样本数：
原文明确描述药物可完全阻断跳动，真实零可能正是表型信号。表内没有逐项 arrest/missing/technical-failure 标识，
所以“至少一份数值”“双重复齐全”“至少一个正值”须分别报告，不能互相替代。

正文 Methods：人 iCell iPSC-CM，24 h 实为给药 23.5 h 后加入检测试剂，再过 30 min 读数，无 HR。
正文 Results 与数据表均为四剂量、双重复；不把 30 min 的六剂量套到本表。
Methods 还明确 DMSO 通常为 0.1%，30 µM 孔为 0.3%，并非全剂量统一 0.3%。
正文称按板 vehicle 归一化，每板有 30–38 vehicle 孔估计控制方差；
本表没有这些独立 vehicle 孔、完整均值/SD，不能据此重建所有控制分布或三分类阈值。
以上依据 [Sirenko 原文 Methods / Results / Fig. 2][sirenko-paper]。

## Heliste：完整药物表，而非仅有最佳剂量

[Dataset S1][heliste-s1]只有 `Sheet1`，693×9；第 4 行为字段，第 5–693 行为 689 个药物标签。
字段为 `Drug`，以及 `CTRL_40nM`、`CTRL_200nM`、`CTRL_1000nM`、`CTRL_5000nM`、
`HYP_40nM`、`HYP_200nM`、`HYP_1000nM`、`HYP_5000nM`。

- 精确标签与去首尾/合并空白/忽略大小写后均为 689 唯一标签；无相同八数值向量的重复行。
- 八列各 689 个有限数值；无空白、非数值、无穷、零、负值或公式单元格。
  因此每个固定剂量、每个氧条件已有 689 个有值标签；HR 与常氧各 2,756 组合。
- 无 CAS、结构键或 library ID；盐型/商品名/别名归并未做。689 是标签数，不是已标准化的化学母体数。
- 标题虽称 raw viability data，表值为 relative viability；没有原始 RLU、板/孔/重复号、SD/SE/CI/n，
  也没有 DMSO/vehicle 行。不能确定每值是单孔、均值还是其他汇总。
- Dataset S2 为 shRNA 表，不计入药物样本；S3 只有 48 个 hits 名称，不代替完整 S1。

[正文 Methods → Drug library screening][heliste-paper]：HL-1 小鼠心房细胞，药物与细胞同时铺板，
24 h 后进行 48 h 缺氧，再 24 h 复氧；平行常氧相同时间。总历时约 96 h。
没有明确写筛选中是否洗脱/换液/补药，记录为“24 h 预处理，HR 期间药物状态未知”，不默认为只给药 24 h。
常氧表也不能当成给药 24–48 h 直接读数。

正文以四剂量最佳效果选每药代表值筛 hits，但 S1 保留所有剂量，已解除“可能只有选后值”的疑问。
另一个疑问仍在：Methods 仅称以 DMSO-treated control cells 归一化，没有明确 HYP/CTRL 各用本氧条件对照，
还是共享常氧对照。Fig. 2b/c 的 HR 对照尺度约 0.78/0.8，而 S1 HYP 多在 100 附近；
缺乏对照原值和明确转换，不能自行用 HYP/CTRL 恢复原始 HR/常氧比值、作者筛选分数或保护阈值。
CTG 是 ATP 相关活性读数，不自动映射为 MTT/CCK8 或直接细胞数。

## 复算方法与独立检查

子代理分别保留 `audit.py`、JSON、逐行 CSV 与来源哈希于上述临时目录；主代理另外独立复算两表，
使用标准库 ZIP/XML 检查 S1，再以 `xlrd` 检查 Table 3，结果一致。
本机已有 `xlrd 2.0.2`、`openpyxl 3.1.5`，未安装或修改环境。
下例从上述原件重新读取关键计数，不执行药物归并、缺失插补或阈值标签化：

```python
from pathlib import Path
import math
import xlrd
import openpyxl

root = Path('/tmp')
book = xlrd.open_workbook(
    root / 'arc-stage3-sirenko-audit/publisher-mmc3.xls', on_demand=True)
for sheet in book.sheets():
    counts = [0, 0]  # 至少一个有限 raw，两个有限 raw
    for row in range(1, sheet.nrows):
        for dose in range(4):
            cells = [sheet.cell(row, 5 + 2 * dose + rep) for rep in range(2)]
            n = sum(c.ctype == xlrd.XL_CELL_NUMBER and math.isfinite(c.value)
                    for c in cells)
            counts[0] += n >= 1
            counts[1] += n == 2
    print(sheet.name, counts)
book.release_resources()

path = root / 'arc-stage3-heliste-audit/41598_2021_96033_MOESM2_ESM.xlsx'
book = openpyxl.load_workbook(path, read_only=True, data_only=True)
rows = list(book.active.values)
data = rows[4:]
labels = {' '.join(r[0].split()).casefold() for r in data}
print('labels', len(data), len(labels))
for column in range(1, 9):
    valid = sum(isinstance(r[column], (int, float))
                and not isinstance(r[column], bool)
                and math.isfinite(r[column]) for r in data)
    print(rows[3][column], valid)
book.close()
```

仅处理小型附件与元数据，无训练、大矩阵或生产入库。机器预检 `free -g` available 为 42 GiB。

## 后续决策的输入

下一项 [确定 Stage 3 首批表型数据的纳入分层与优先来源](https://github.com/Moke110/arc/issues/24)
可基于真实表来讨论：严格 24 h、数据粒度较好的 iPSC-CM 单指标，还是 HR 相关但时长/对照有缺口的 HL-1 数据，
以及两者是否分层保留。没有在本票选择终点或放宽门槛。

独立样本与跨源去重、零/缺失处理、方向标签阈值及归一化证据分别由已有后续决策票承接。
需要联系作者澄清时先取得用户指示；本轮未联系作者。

[sirenko-paper]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3900303/
[sirenko-table]: https://ars.els-cdn.com/content/image/1-s2.0-S0041008X13004171-mmc3.xls
[sirenko-auc]: https://ars.els-cdn.com/content/image/1-s2.0-S0041008X13004171-mmc4.xls
[sirenko-bmc]: https://ars.els-cdn.com/content/image/1-s2.0-S0041008X13004171-mmc5.xls
[sirenko-json]: https://pmc-oa-opendata.s3.amazonaws.com/PMC3900303.1/PMC3900303.1.json
[heliste-paper]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8371130/
[heliste-s1]: https://pmc-oa-opendata.s3.amazonaws.com/PMC8371130.1/41598_2021_96033_MOESM2_ESM.xlsx
[heliste-s2]: https://pmc-oa-opendata.s3.amazonaws.com/PMC8371130.1/41598_2021_96033_MOESM3_ESM.xlsx
[heliste-s3]: https://pmc-oa-opendata.s3.amazonaws.com/PMC8371130.1/41598_2021_96033_MOESM4_ESM.xlsx
[heliste-json]: https://pmc-oa-opendata.s3.amazonaws.com/PMC8371130.1/PMC8371130.1.json
