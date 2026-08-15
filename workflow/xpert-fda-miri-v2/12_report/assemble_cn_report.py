#!/usr/bin/env python
"""Generate 12_report/output/report_cn.md — Chinese translation of the drug report.

Changes vs report.md:
  * Full Chinese translation.
  * First (summary) table drops the NMPA column and adds an AMI/MIRI column
    formatted "n/m": n = drug AMI-improvement evidence, m = drug MIRI-improvement
    evidence, with 2=clinical improvement, 1=preclinical improvement,
    0=contradictory reports/ineffective, -1=preclinical worsening, -2=clinical worsening.

Covers the 12_report final-drug set (10_filter survivors with integrated-score
percentile > 0.5), read via run_drug_agents.load_final_drugs().
"""
import json
import re
from datetime import date
from pathlib import Path

import pandas as pd

from run_drug_agents import load_final_drugs

WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = WORKFLOW_ROOT / "12_report"
CACHE_DIR = REPORT_DIR / "cache"
OUT = REPORT_DIR / "output" / "report_cn.md"

# (n=AMI evidence, m=MIRI evidence) — scored from the A1 evidence (see report notes)
SCORES = {
    "MAPROTILINE": (0, 0), "ALFUZOSIN": (0, 0), "TRIFLUPROMAZINE": (0, 0),
    "BROMPHENIRAMINE": (0, 0), "MIANSERIN": (0, 1), "LEVALLORPHAN": (0, 0),
    "DROPERIDOL": (0, 1), "LABETALOL": (0, 1), "PINDOLOL": (0, 1),
    "NALTREXONE": (0, -1), "LOXAPINE": (0, 0), "PROMAZINE": (0, 0),
    "DILEVALOL": (0, 0), "CHLORPHENIRAMINE": (0, 1), "MIRTAZAPINE": (0, 1),
    "PHENOXYBENZAMINE": (0, 0), "CYCLIZINE": (-2, 0), "CHLORPROMAZINE": (0, 0),
    "AMBRISENTAN": (0, 0), "BENFLUOREX": (0, 0), "LORATADINE": (2, 0),
    "ANTAZOLINE": (0, 0), "HALOPERIDOL": (0, 0), "CYPROHEPTADINE": (0, 1),
    "TRIMIPRAMINE": (-1, 0), "TRAZODONE": (0, 1), "TERBUTALINE": (0, 1),
    "DESIPRAMINE": (0, -1), "PROPRANOLOL": (2, 1), "PROPAFENONE": (-1, 0),
    "PROCAINAMIDE": (0, 0), "TAMSULOSIN": (0, 1), "PHENTOLAMINE": (0, 0),
    "ATOMOXETINE": (0, 0), "DYCLONINE": (0, 0), "DULOXETINE": (0, 0),
    "SOTALOL": (-2, 0), "SULCONAZOLE": (0, 0), "PIOGLITAZONE": (2, 1),
    "CABERGOLINE": (0, 1), "RIMONABANT": (0, 1), "ALPIDEM": (0, 0),
    "BUFLOMEDIL": (0, 1),
}

# Chinese content per drug, authored from the A1 evidence.
# cn: 分类 | cn_scope: 临床用途 | cn_us: 美国FDA | cn_clin: MIRI/AMI临床 | cn_pre: 临床前
CN = {
"MAPROTILINE": dict(cn="四环类抗抑郁药（TeCA）；选择性去甲肾上腺素再摄取抑制剂（NRI）；ATC N06AA21；靶点 SLC6A2；CYP2D6 底物",
  cn_scope="抗抑郁药，用于重性抑郁障碍/心境恶劣；多用于心脏病人群情绪障碍、慢性疼痛辅助",
  cn_us="1980 年获 FDA 批准（抗抑郁）；美国市场已停售，OpenFDA 无现行标签",
  cn_clin="仅 AMI 后抑郁情绪管理证据：126 例急性心梗患者中 25mg bid 改善抑郁、无 ECG 心脏毒性；GPRD 病例对照 AMI 风险趋势 OR 0.59（不显著）；无直接改善 AMI/MIRI 的证据",
  cn_pre="犬心梗后模型显示传导减慢轻于阿米替林、致心律失常风险低（安全性好），但无直接心肌 I/R 保护作用；无 MIRI 动物研究"),
"ALFUZOSIN": dict(cn="喹唑啉类 α1-肾上腺素受体拮抗剂（尿选择性）；ATC G04CA01",
  cn_scope="良性前列腺增生/下尿路症状",
  cn_us="2003 年获批治疗 BPH",
  cn_clin="无 AMI/MIRI 临床证据；仅 6523 例 BPH 患者心血管耐受性研究（合并缺血性心脏病者耐受良好）",
  cn_pre="无直接 MIRI/AMI 动物研究；TxGNN 预测为阴性；同类喹唑啉 α1 拮抗剂在离体心脏中效果不一"),
"TRIFLUPROMAZINE": dict(cn="吩噻嗪类典型抗精神病药（脂肪侧链）；ATC N05AA05；多巴胺 D2 拮抗剂",
  cn_scope="抗精神病（精神病性障碍）、止吐",
  cn_us="1957 年获批（历史性）；现美国无现行标签",
  cn_clin="无 AMI/MIRI 临床证据；知识图谱将其列为冠心病/冠脉血栓禁忌",
  cn_pre="无阳性心脏保护证据；TxGNN 对缺血再灌注/心梗预测为阴性"),
"BROMPHENIRAMINE": dict(cn="第一代烷基胺类 H1 受体拮抗剂（抗组胺药）；ATC R06AB01",
  cn_scope="OTC 抗组胺药：过敏性鼻炎、感冒症状",
  cn_us="1957 年获批，OTC 销售",
  cn_clin="无直接 AMI/MIRI 临床证据；心脏肥大细胞综述仅将氯苯那敏/西替利嗪列为潜在 H1 保护药",
  cn_pre="间接：菲尼拉敏类（同类）保护肠系膜 I/R；相关 H1 抗组胺药（氯苯那敏）减轻离体大鼠心 I/R 心律失常；TxGNN 阴性"),
"MIANSERIN": dict(cn="四环类抗抑郁药（NaSSA）；α2 拮抗、5-HT2A/2C/3 拮抗；ATC N06AX03",
  cn_scope="抗抑郁药（重性抑郁），主要用于欧洲/日本/亚洲",
  cn_us="未获美国 FDA 批准",
  cn_clin="无 MIANSERIN 用于 MIRI/AMI 的临床试验（心梗后抑郁试验用的是米氮平）",
  cn_pre="大鼠直接阳性：离体工作心脏（1-10μM）抗再灌注心律失常并改善再灌注功能恢复；冠状动脉结扎大鼠降低死亡率、室颤持续时间与再灌注心律失常（奎尼丁样心脏抑制作用）"),
"LEVALLORPHAN": dict(cn="阿片受体拮抗剂（吗啡烷衍生物，激动/拮抗）；靶点 OPRM1/OPRK1/OPRD1",
  cn_scope="麻醉药拮抗剂、呼吸兴奋剂",
  cn_us="未检索到美国 FDA 现行记录（ChEMBL 记录 1982 年首次批准）",
  cn_clin="无 AMI/MIRI 临床试验或心梗后应用记录",
  cn_pre="仅间接：脊髓/脑缺血模型有保护（κ-阿片介导）；新生大鼠心肌细胞增强收缩幅度；TxGNN 对 I/R、心梗预测为阴性"),
"DROPERIDOL": dict(cn="丁酰苯类抗精神病药；D2 样多巴胺拮抗剂；ATC N05AD08",
  cn_scope="止吐（围手术期恶心呕吐）、神经安定镇痛",
  cn_us="1970 年获批（注射剂）；现行标签为手术相关恶心呕吐",
  cn_clin="矛盾：8 例急性前壁透壁心梗神经安定镇痛（芬太尼+氟哌利多）疼痛缓解、ST 段抬高总和下降 57%；但 112 例不稳定心绞痛 RCT 显示该混合液显著恶化结局（住院不稳定 OR 5.93，12 个月 AMI 发生 OR 3.57）",
  cn_pre="直接阳性：犬/猪冠脉血栓模型中（0.2-0.8mg/kg）消除反复冠脉血流下降、抑制血小板聚集；急性冠脉闭塞模型提高侧支循环、抗心律失常、改善缺血心肌血供与收缩力"),
"LABETALOL": dict(cn="α、β 肾上腺素受体双重拮抗剂（降压药）；ATC C07AG01",
  cn_scope="高血压（含妊娠期高血压/先兆子痫）",
  cn_us="1984 年获批治疗高血压",
  cn_clin="166 例疑似急性心梗 RCT：降压降心率但 CKMB 释放显著增加、R 波评分/射血分数无改善，作者结论为不推荐常规用于正常血压疑似梗死患者（即无 AMI 改善证据）",
  cn_pre="兔冠脉闭塞/再灌注模型（1mg/kg）恢复被脑死亡去甲肾上腺素风暴消除的缺血预适应心脏保护（梗死面积约 50%→13.9%）；TxGNN 阴性"),
"PINDOLOL": dict(cn="非选择性 β 受体阻滞剂（具有内在拟交感活性 ISA）；ATC C07AA03；II 类抗心律失常药",
  cn_scope="高血压",
  cn_us="1982 年获批（美国在售）治疗高血压",
  cn_clin="521 例心梗后患者有抗心律失常/抗心绞痛作用；AMI 中降低心脏指数；但文献指出具 ISA 的 β 阻滞剂（如吲哚洛尔）降低心肌缺血发病/死亡获益差、二级预防获益少",
  cn_pre="大鼠离体缺血-再灌注心脏（50μM）加速机械功能与能量状态恢复、抑制游离脂肪酸蓄积；（20μM）保护膜磷脂、降低脂质过氧化与 CK 释放"),
"NALTREXONE": dict(cn="阿片受体纯拮抗剂（μ/κ/δ）；ATC N07BB04",
  cn_scope="酒精依赖、阿片拮抗",
  cn_us="1984 年获批（酒精依赖）",
  cn_clin="无 AMI/MIRI 临床证据（ClinicalTrials.gov 零试验；心血管数据主要来自纳曲酮-安非他酮组合的 CV 安全性试验）",
  cn_pre="在心肌 I/R 模型中消除阿片介导及预适应/适应介导的心脏保护而非自身保护：大鼠冠脉闭塞/再灌注研究（5mg/kg）阻断阿片激动剂与慢性低氧适应的缩梗死效应 → 临床前恶化"),
"LOXAPINE": dict(cn="二苯并氧氮卓类典型抗精神病药；D2 拮抗；ATC N05AH01",
  cn_scope="精神分裂症（口服/吸入）",
  cn_us="1975 年获批口服治疗精神分裂症；2012 年吸入剂（Adasuve）获批",
  cn_clin="无 AMI/MIRI 临床证据；QT 研究中性；知识图谱将心脏病列为禁忌；过量可致室速",
  cn_pre="无直接心肌 I/R 证据；TxGNN 预测为阴性；唯一 I/R 研究在视网膜（轻度改善 B 波恢复）"),
"PROMAZINE": dict(cn="吩噻嗪类典型抗精神病药；D2 拮抗；ATC N05AA03",
  cn_scope="抗精神病、镇静",
  cn_us="（历史性精神安定药）",
  cn_clin="无 AMI/MIRI 临床试验；仅过量中毒病例报告出现一过性心肌缺血样 ECG 改变（毒性而非获益）",
  cn_pre="TxGNN 对心梗/心肌缺血/I-R 无评分；未见心脏保护模型；大鼠肺 I/R 中作为 ecto-ATPase 底物作用（无心脏保护）"),
"DILEVALOL": dict(cn="非选择性 β 受体阻滞剂（拉贝洛尔的 R,R 异构体，具血管扩张/ISA）",
  cn_scope="降压（未广泛上市）",
  cn_us="未检索到现行美国记录（拉贝洛尔母体在售）",
  cn_clin="无 DILEVALOL 特异性临床证据（仅检索到母体拉贝洛尔数据）",
  cn_pre="阴性：兔离体工作心（25min 缺血/30min 再灌注）功能恢复无改善；猪 85-90% LAD 闭塞模型 300μg/kg 无抗缺血作用；TxGNN 未入前 100"),
"CHLORPHENIRAMINE": dict(cn="第一代烷基胺类 H1 受体拮抗剂；ATC R06AB04",
  cn_scope="OTC 抗组胺药：过敏性鼻炎、感冒症状",
  cn_us="1957 年获批，OTC 销售",
  cn_clin="无直接临床试验；但 H1-组胺轴与急性冠脉综合征相关（防止组胺诱导组织因子表达，与血栓形成相关；Kounis 综合征过敏心梗中抗组胺药的使用）",
  cn_pre="大鼠离体 Langendorff 心脏直接阳性（10μM）：减轻 I/R 心功能不全、提高左室收缩压、完全抑制再灌注心律失常、增加再灌注冠脉血流；综述提出氯苯那敏为心脏保护药（预防缺血再灌注中肥大细胞脱颗粒）"),
"MIRTAZAPINE": dict(cn="去甲肾上腺素能及特异性 5-羟色胺能抗抑郁药（NaSSA）；ATC N06AX11",
  cn_scope="重性抑郁障碍",
  cn_us="获批抗抑郁药",
  cn_clin="心梗后抑郁：MIND-IT 子研究显示米氮平/西酞普兰治疗无反应者心脏事件增加（HR 2.66-2.92）；贝叶斯荟萃分析示抗抑郁药降低 ACS 后重复住院（OR 0.62）与复发心梗（OR 0.45）——属抑郁治疗间接效应，非直接 AMI 改善",
  cn_pre="小鼠急性心梗→慢性心衰模型验证心脏保护（1 周早期与 4 周慢性期均有效，GRB14/RAB3A 参与）；大鼠肾 I/R 减轻氧化应激（MDA/MPO↓、GSH/GST 保存）；TxGNN 预测阴性"),
"PHENOXYBENZAMINE": dict(cn="非选择性长效 α 肾上腺素受体拮抗剂；ATC C04AX",
  cn_scope="嗜铬细胞瘤血压控制、雷诺现象",
  cn_us="获批（嗜铬细胞瘤）",
  cn_clin="无专门临床试验；仅嗜铬细胞瘤合并急性心梗的病例报告用其控制血压",
  cn_pre="兔冠脉闭塞/再灌注预适应模型（4mg/kg）既不阻断预适应缩梗死效应、自身也不改变梗死面积（无直接保护）；火鸡扩张性心肌病模型中保护作用弱于 β 阻滞剂/钙拮抗剂"),
"CYCLIZINE": dict(cn="哌嗪类 H1 受体拮抗剂/止吐药；ATC R06AE03",
  cn_scope="恶心呕吐（术后、晕动病、眩晕）、过敏",
  cn_us="1954 年获批（Marezine，现停药）",
  cn_clin="AMI/MIRI 中为有害：Lancet 1988 示静注赛克利嗪升高体/肺循环与心室充盈压、抵消二醋吗啡扩静脉作用，警示勿用于心梗/左心衰；Emerg Med J 综述建议急性冠脉事件中避免使用；列为心梗禁忌",
  cn_pre="无直接 MIRI 保护证据（'cyclizine AND reperfusion' 0 结果）；早期动物药理仅描述抗心律失常/局麻活性，无 I/R 心脏保护；TxGNN 阴性"),
"CHLORPROMAZINE": dict(cn="吩噻嗪类典型抗精神病药；D2 拮抗；ATC N05AA01",
  cn_scope="抗精神病、止吐、镇静",
  cn_us="1954 年获批（首个抗精神病药）",
  cn_clin="无心脏 AMI/MIRI 临床证据；唯一 I/R 临床工作为脑卒中（氯丙嗪+异丙嗪冬眠疗法，CHILL-S 试验）",
  cn_pre="多为机制/工具性：作为钙调蛋白抑制剂用于增强心肌保存；抑制缺血心肌细胞对促凋亡心脏成纤维细胞外泌体的内吞（减轻凋亡、心梗表型更小，Br J Pharmacol）；非直接心脏保护"),
"AMBRISENTAN": dict(cn="选择性内皮素 A 型受体（ETA）拮抗剂；ATC C02KX02",
  cn_scope="肺动脉高压",
  cn_us="获批治疗肺动脉高压",
  cn_clin="研究性：ENDORA-PCI 试验（NSTEACS 行 PCI 前使用以减轻围手术期微血管损伤/围手术期心梗，进行中，无结果）",
  cn_pre="支持但间接：幼年大鼠肠 I/R 模型减轻心脏酶学升高、氧化应激（MDA）、NF-κB、caspase-3 并改善心肌纤维组织学；大鼠肾 I/R 保护；无直接心脏 I/R 研究"),
"BENFLUOREX": dict(cn="苯乙胺类食欲抑制剂/降糖辅助用药；ATC A10BX06",
  cn_scope="（曾用于糖尿病辅助，已撤市）",
  cn_us="未检索到现行记录（安全性问题撤市）",
  cn_clin="无 AMI/MIRI 临床证据；安全性信号为瓣膜性心脏病、肺动脉高压、室颤/心脏骤停",
  cn_pre="无心肌 I/R 或心脏保护研究；TxGNN 未入前列"),
"LORATADINE": dict(cn="第二代非镇静 H1 受体拮抗剂；ATC R06AX13",
  cn_scope="过敏性鼻炎、荨麻疹",
  cn_us="1993 年获批，OTC",
  cn_clin="小型 RCT 临床阳性：20 例急性下壁心梗稳定期患者加用氯雷他定 10mg/d×7 天显著改善运动负荷试验缺血参数（最大 ST 压低、受累导联数、总 ST 评分下降）；无更大结局试验",
  cn_pre="无氯雷他定特异心肌 I/R 动物研究；机制性：H1 抗组胺药/肥大细胞稳定剂可能通过防止肥大细胞脱颗粒限制 I/R 心肌损伤（大鼠模型）；地氯雷他定经 Nrf2-GCL-GSH 减轻大鼠肾 I/R；TxGNN 阴性"),
"ANTAZOLINE": dict(cn="第一代 H1 抗组胺药（咪唑啉类）；具局麻与 Vaughan-Williams I 类抗心律失常活性",
  cn_scope="抗组胺、眼科局麻（历史性）",
  cn_us="历史性用药",
  cn_clin="无直接 AMI/MIRI 试验；仅急诊房颤复律回顾性研究（合并冠心病患者安全有效，既往心梗不影响转复成功率）",
  cn_pre="无直接心肌 I/R 或心脏保护证据；多通道抗心律失常特性（抑制 hERG/hKv）、脑缺血 NMDA 阻断神经保护；TxGNN 对所有缺血/心梗疾病为阴性"),
"HALOPERIDOL": dict(cn="丁酰苯类典型抗精神病药；D2 拮抗；ATC N05AD01",
  cn_scope="精神分裂症、急性激越、谵妄",
  cn_us="获批抗精神病药",
  cn_clin="猪 AMI 研究：冠脉闭塞前给药不加重 AMI 心律失常（低剂量室颤 3/11 vs 对照 7/11，不显著），高剂量延长 QT；医院 QTc-DDI 研究示其（如与胺碘酮联用）为急性心梗患者药物性 QTc 延长常见原因；观察队列示抗精神病药使用与心血管风险+32%相关",
  cn_pre="阳性研究主要针对衍生物 F2（N-正丁基氟哌啶醇碘化物）：小鼠心肌 I/R 模型（0.4mg/kg 再灌注前）经 SIRT1-Nrf2 及 AMPK/FoxO1 通路减轻损伤；氟哌啶醇本身抑制铁死亡（机制提示与 I/R 相关），但无自身直接心脏 I/R 数据"),
"CYPROHEPTADINE": dict(cn="第一代抗组胺药；H1 拮抗、抗 5-HT（5-HT2）；ATC R06AX02",
  cn_scope="过敏、偏头痛预防、食欲刺激",
  cn_us="获批（过敏/荨麻疹）",
  cn_clin="无完成的心梗 RCT；进行中 CYPRO-MR 2 期试验（NCT05469165）测试 4mg 口服减轻心梗后缺血性二尖瓣反流（经 5-HT2B 阻断）；羊转译模型示心梗后降低 5-HT 与二尖瓣反流",
  cn_pre="直接阳性：PMID 7976381（1994）赛庚啶 2.5/5μM 保护离体 Langendorff 大鼠心脏免于低流量 I/R 损伤（LDH/CK 释放↓、SOD/GSH-Px 保存、MDA↓、心功能改善、心律失常/停搏减少）；另一中国研究示对钙反常心肌损伤保护"),
"TRIMIPRAMINE": dict(cn="三环类抗抑郁药（TCA）；非选择性单胺再摄取抑制剂；ATC N06AA06",
  cn_scope="重性抑郁障碍",
  cn_us="获批抗抑郁药",
  cn_clin="FDA 标签明确禁忌：心肌梗死急性恢复期禁用；GPRD 病例对照（8688 例首次 AMI）无显著 AMI 风险关联（OR 0.92）；轻度心脏病抑郁患者中有效且心脏不良反应轻微",
  cn_pre="无直接心肌 I/R 研究；TxGNN 对心肌缺血/I-R/心梗预测未入前列（部分精神科适应证入列）"),
"TRAZODONE": dict(cn="5-羟色胺拮抗与再摄取抑制剂（SARI）；ATC N06AX05",
  cn_scope="抑郁障碍（含伴失眠）",
  cn_us="获批抗抑郁药",
  cn_clin="矛盾：巢式病例对照（24155 例 AMI）示当前使用曲唑酮 AMI 风险下降（AOR 0.76，剂量依赖）；但日本 JADER 数据库示其危及生命室性心律失常 ROR 最高（17.1）、QT 延长/TdP（aROR 12.8）",
  cn_pre="小鼠全脑 I/R 减轻氧化损伤与线粒体酶功能障碍（NO-cGMP 通路）；机制相关：5-HT2A/5-HT2 拮抗在离体大鼠心及大鼠/豚鼠心脏 I/R 模型具保护作用（cinanserin/ketanserin/LY53857 降低 LDH、改善收缩恢复；sarpogrelate 缩梗死、保 ATP）"),
"TERBUTALINE": dict(cn="选择性 β2-肾上腺素受体激动剂（支气管扩张药）；ATC R03AC03/R03CC03",
  cn_scope="哮喘/COPD 支气管痉挛",
  cn_us="获批（吸入/口服支气管扩张剂）",
  cn_clin="矛盾：6 例心源性休克 AMI 患者输注改善血压/心脏指数/肺楔压/尿量；美托洛尔+特布他林对合并 COPD 的 AMI 患者血流动力学有利；但 FDA 黑框警告将心肌缺血列为不良反应",
  cn_pre="β2-AR 激活可增加心脏对 I/R 的耐受；选择性 β2 激动剂 arformoterol 在冠脉闭塞动物模型缩梗死；再灌注时 β2-AR 刺激经脾脏 IL-10 抗炎通路限制再灌注损伤（小鼠 LAD 结扎）"),
"DESIPRAMINE": dict(cn="三环类抗抑郁药（TCA）；非选择性单胺再摄取抑制剂；ATC N06AA01",
  cn_scope="重性抑郁障碍",
  cn_us="1964 年获批（抗抑郁）",
  cn_clin="无直接临床试验；台湾 NHIRD 病例交叉研究探讨抗抑郁药与 AMI 住院；药遗传学研究发现治疗可正常化血液黏度（可能心血管有利）；BDNF Val66Met 小鼠中减少动脉粥样硬化—均为间接",
  cn_pre="有害：Wistar 大鼠再灌注期 0.8mg/kg 使梗死面积 32.0%→46.1%，并消除缺血后处理保护；Kardiologiia 研究示再灌注期去甲肾上腺素再摄取阻断加重心肌损伤并消除预适应保护 → 临床前恶化"),
"PROPRANOLOL": dict(cn="非选择性 β 肾上腺素受体拮抗剂（β 阻滞剂）；ATC C07AA05；II 类抗心律失常",
  cn_scope="高血压、心绞痛、心梗后、心律失常、偏头痛",
  cn_us="获批：用于降低心梗急性期后心血管死亡率（FDA 标签）；AHA/ACC 指南推荐",
  cn_clin="强临床证据：BHAT 试验示心梗后普萘洛尔全因死亡下降约 26%；FDA 标签适应证、指南一级推荐（β 阻滞剂降低心梗后再梗死与死亡）",
  cn_pre="系统综述示普萘洛尔在冠脉闭塞动物模型缩梗死、可预防再灌注心脏损伤；兔 I/R 以等效 β 阻断剂量减少凋亡心肌细胞约 39%；离体兔心保存再灌注功能恢复；S-普萘洛尔经 RACK1 下调保护 H9C2 细胞免于 I/R 凋亡"),
"PROPAFENONE": dict(cn="Ic 类抗心律失常药（钠通道阻滞）；ATC C01BC03；弱 β 阻滞活性",
  cn_scope="室上性/室性心律失常",
  cn_us="获批抗心律失常药",
  cn_clin="Ic 类抗心律失常荟萃分析讨论 CAST（心梗后伴室早+射血分数降低者死亡率增加）；普帕酮在 CAST 型心梗后人群仍属禁忌（潜在临床危害）",
  cn_pre="犬冠脉闭塞/再灌注研究中 diprafenone（强效类似物）减少急性闭塞期室颤但未减少再灌注心律失常、部分动物室速增强；直接普帕酮特异心脏保护证据缺乏"),
"PROCAINAMIDE": dict(cn="IA 类抗心律失常药；钠通道阻滞；ATC C01BA02",
  cn_scope="室性/室上性心律失常",
  cn_us="获批抗心律失常药",
  cn_clin="不治疗心梗本身；AMI 中仅用于抑制缺血并发的室速（10mg/kg 终止单形性室速优于利多卡因 12/15 vs 3/14）；FDA 标签警示按 CAST 不要抑制无症状非致命性心梗后心律失常",
  cn_pre="冠脉闭塞-再灌注心律失常模型（大鼠/犬/猪）中 20mg/kg 预防大鼠、犬（非猪）室颤；无直接缩梗死/减轻再灌注损伤证据；TxGNN 未入前 20"),
"TAMSULOSIN": dict(cn="选择性 α1A/α1D 肾上腺素受体拮抗剂；ATC G04CA02",
  cn_scope="良性前列腺增生",
  cn_us="1997 年获批治疗 BPH",
  cn_clin="无 AMI/MIRI 适应证或试验；仅 BPH+2 型糖尿病队列：5α-还原酶抑制剂相对坦索罗辛心梗风险更高（HR 1.20-1.33）；无意外心血管安全性信号",
  cn_pre="直接阳性：大鼠异丙肾上腺素诱导心梗模型口服 0.8mg/kg×21 天显著降低心脏/体重指数与 CK-MB、抗氧化损伤（GSH/SOD↑、脂质过氧化物↓）、降低 TGF-β/p-Smad2/3、增强 ILK 表达（经 ILK/TGF-β/Smad 通路保护）"),
"PHENTOLAMINE": dict(cn="非选择性 α 肾上腺素受体拮抗剂（咪唑啉类）；外周血管扩张/降压",
  cn_scope="嗜铬细胞瘤高血压危象、可乐定相关高血压",
  cn_us="获批但 FDA 标签禁忌用于心肌梗死/冠心病",
  cn_clin="历史性用于急性心梗左心衰（12 例改善心指数、降低肺楔压）；可卡因相关 AMI 中 α 阻滞；但 FDA 标签禁忌用于心梗/冠心病，警示低血压后可致心梗/脑血管事件",
  cn_pre="无支持性保护证据：兔 I/R 非选择性 α 阻滞不消除预适应缩梗死；羊 CPB 模型酚妥拉明消除 CPB 诱导的预适应肌保护；犬冠脉闭塞-再灌注仅减少单个室早；TxGNN 阴性"),
"ATOMOXETINE": dict(cn="选择性去甲肾上腺素再摄取抑制剂（NRI）；ATC N06BA09",
  cn_scope="注意缺陷多动障碍（ADHD）",
  cn_us="2002 年获批治疗 ADHD",
  cn_clin="无 AMI/MIRI 治疗试验；青年女性病例报告阿托莫西汀诱发心梗；大型队列无 MI 风险升高（HR 0.56；JAMA 2011 无严重心血管事件）；无 MI 药警戒信号",
  cn_pre="无直接心脏 I/R 研究；TxGNN 对 AMI/I-R 未入前 20；仅脑缺血模型有神经保护（沙鼠 CA1 神经元、新生鼠缺氧缺血）"),
"DYCLONINE": dict(cn="局部麻醉药（哌啶丙苯酮类）；ATC N01BX02/R02AD04",
  cn_scope="表面麻醉（黏膜）",
  cn_us="OTC/处方局麻药",
  cn_clin="阴性：PubMed/ClinicalTrials.gov 检索无任何 MIRI/AMI 临床试验、心梗后应用或流行病学证据",
  cn_pre="阴性：无心肌 I/R 保护数据；TxGNN 对 I-R/心梗/心肌缺血未入前 200；机制上抑制 ALDH2（其激活通常具心脏保护），不支持 MIRI/AMI 获益"),
"DULOXETINE": dict(cn="5-羟色胺-去甲肾上腺素再摄取抑制剂（SNRI）；ATC N06AX21",
  cn_scope="抑郁障碍、广泛性焦虑、神经病理性疼痛",
  cn_us="2004 年获批",
  cn_clin="无 MIRI/AMI 治疗作用：Medicare 队列（34009 度洛西汀 vs 233060 加巴喷丁）AMI/卒中/死亡相当（HR 0.98）；SNRI 对心梗后抑郁治疗存在证据缺口；轻度升心率/舒张压",
  cn_pre="抗血小板/抗血栓特性（浓度依赖性抑制人/鼠血小板聚集、延长小鼠闭塞时间）；相关 SNRI 去甲文拉法辛在大鼠 LAD 闭塞-再灌注中显著缩梗死；TxGNN 预测阴性"),
"SOTALOL": dict(cn="非选择性 β 阻滞剂+III 类抗心律失常药（钾通道/hERG 阻滞）；ATC C07AA07/C01BD01",
  cn_scope="室性/室上性心律失常",
  cn_us="获批抗心律失常药",
  cn_clin="临床恶化：SWORD 试验（心梗后 EF≤40% 患者 d-索他洛尔）因死亡率增加提前终止（5.0% vs 3.1%，RR 1.65，主要为心律失常性死亡）；急性/心梗后及心源性休克禁用；VANISH2 中作为一线抗心律失常但导管消融更优",
  cn_pre="大鼠冠脉闭塞/再灌注心律失常模型（5mg/kg+10mg/kg/h）显著减少室性心律失常（7/20 vs 16/20）、缩短 VT 时长；但离体鼠心 I/R 中 β 阻滞阻断温度预适应心脏保护（β-肾上腺素能信号参与预适应）"),
"SULCONAZOLE": dict(cn="咪唑类抗真菌药（麦角甾醇合成抑制剂）",
  cn_scope="局部真菌感染",
  cn_us="获批局部抗真菌药",
  cn_clin="无任何 AMI/MIRI 临床证据（外用抗真菌药，无心血管适应证）",
  cn_pre="TxGNN 对心梗/心肌缺血/I-R/STEMI/ACS 均未入列；无心肌 I/R 动物研究；有抗炎活性但与心脏保护无关"),
"PIOGLITAZONE": dict(cn="噻唑烷二酮类胰岛素增敏剂；PPAR-γ 激动剂；ATC A10BG03",
  cn_scope="2 型糖尿病",
  cn_us="1999 年获批（2 型糖尿病）",
  cn_clin="临床改善：IRIS 试验（3876 例既往卒中/TIA 伴胰岛素抵抗）降低复发卒中/心梗（HR 0.70）与 ACS（HR 0.72）；PROactive 试验降低致死/非致死心梗等大血管事件；综述示其在一级/二级预防中一致降低心梗风险",
  cn_pre="大鼠 I/R 预适应（3mg/kg IV）显著缩小梗死、减少凋亡与 CK-MB/cTnI/LDH 释放、降低炎症因子（IL-1β/IL-6/TNF-α/NF-κB/Bax），升高 Bcl-2/HSP90（HSP90 依赖、PPAR-γ 介导）；PPAR-γ 纳米药改善 AMI 后心脏修复"),
"CABERGOLINE": dict(cn="麦角碱衍生物；长效多巴胺 D2 受体激动剂（泌乳素抑制剂）；ATC G02CB03/N04BC06",
  cn_scope="高泌乳素血症、帕金森病",
  cn_us="1996 年获批（高泌乳素血症）",
  cn_clin="无心梗/MIRI 临床试验；病例报告青年女性服用卡麦角林（泌乳素瘤）发生多支自发性冠脉夹层（SCAD）呈现 ACS；与瓣膜反流相关而非心脏保护",
  cn_pre="小鼠网络分析研究示早期心梗→慢性心衰转变期潜在心脏保护；大鼠心脏 D2 激动（溴隐亭）增强缺血后处理保护（PKC-ε 转位、缩梗死/凋亡），为卡麦角林 D2 激动提供机制依据；TxGNN 阴性"),
"RIMONABANT": dict(cn="选择性大麻素 CB1 受体拮抗剂/反向激动剂；减肥药；ATC A08AX01",
  cn_scope="（曾用于肥胖，已撤市）",
  cn_us="未获批美国（CRESCENDO 因神经精神不良事件提前终止后撤回申请）",
  cn_clin="CRESCENDO 试验（18695 例血管高危患者）示利莫那班 20mg 对主要血管事件无获益（HR 0.97，NS），因神经精神事件/自杀提前终止；无人类心梗后适应证",
  cn_pre="阳性：C57BL/6J 小鼠急性心梗 10mg/kg 缩小梗死（29.6% vs 49.8% 对照），CB1-/- 小鼠效应消失；大鼠冠脉永久结扎 10mg/kg/d 改善收缩/舒张功能与重构、降低纤维化/TGF-β1/MMP-9/BNP"),
"ALPIDEM": dict(cn="咪唑并吡啶类非苯二氮卓镇静催眠/抗焦虑药；GABA-A ω1(BZ1) 部分激动剂",
  cn_scope="（曾用于广泛性焦虑，法国短暂上市后因肝毒性撤市）",
  cn_us="未获批美国；欧洲撤市",
  cn_clin="无 AMI/MIRI 临床证据（仅心血管术前情境性焦虑的间接使用，无心结局数据）",
  cn_pre="无心脏缺血再灌注动物或机制研究；TxGNN 药物知识图谱中无 alpidem 条目，无法生成重定位预测"),
"BUFLOMEDIL": dict(cn="外周血管扩张药；α1 拮抗+钙通道阻滞；ATC C04AX20",
  cn_scope="外周动脉闭塞性疾病（间歇性跛行）",
  cn_us="获批外周血管疾病（历史性）",
  cn_clin="无直接 MIRI/AMI 临床试验；冠心病综述将口服血液流变学药物（丁咯地尔/己酮可可碱）列为顽固性心绞痛辅助、可能改善 AMI 溶栓后缺血区再灌注；TxGNN 对 ACS 排第 130/200（低优先级）",
  cn_pre="犬 LAD 闭塞 15min/再灌注 1h 模型：预处理 10mg/kg 提高射血分数（89% vs 69%）、缩小功能失调区；仓鼠骨骼肌 I/R 减轻毛细血管灌注衰竭、白细胞黏附与微血管大分子渗漏"),
}

DETAIL_LABELS = {
    "ENGLISH_NAME": "英文名", "CHINESE_NAME": "中文名",
    "US_FDA_APPROVAL": "美国 FDA 审批", "CHINA_NMPA_APPROVAL": "中国 NMPA 审批",
    "CLASSIFICATION": "药理分类", "CLINICAL_SCOPE": "临床应用范围",
    "MIRI_AMI_CLINICAL": "MIRI/AMI 临床关联", "MIRI_AMI_PRECLINICAL": "MIRI/AMI 临床前关联",
}

SCORE_MEANING = {
    2: "临床改善", 1: "临床前改善", 0: "矛盾报道/无效", -1: "临床前恶化", -2: "临床恶化",
}

CHINESE_OVERRIDE = {
    "LEVALLORPHAN": "左啡诺 (烯丙左吗喃)",
    "LOXAPINE": "洛沙平",
    "CYCLIZINE": "赛克利嗪 (环利嗪)",
    "DESIPRAMINE": "地昔帕明 (去甲丙咪嗪)",
    "ALPIDEM": "阿吡坦",
}


def parse_solution(text):
    if not text:
        return {}
    out = {}
    for k in DETAIL_LABELS:
        m = re.search(rf"(?m)^{re.escape(k)}\s*:\s*(.*)$", text)
        out[k] = m.group(1).strip() if m else ""
    return out


def main():
    ranked = load_final_drugs()
    ranked["integrated"] = ranked["integrated"].round(4)
    n_drugs = len(ranked)

    rows = []
    quick = []
    details = []
    for i, (_, r) in enumerate(ranked.iterrows(), 1):
        drug = r["drug_name"]
        d = json.load(open(CACHE_DIR / f"{drug}.json"))
        sol = parse_solution(d.get("solution") or "")
        n, m = SCORES[drug]
        cn = CN[drug]
        zh = sol["CHINESE_NAME"]
        if (not zh or zh.lower() == "not found") and drug in CHINESE_OVERRIDE:
            zh = CHINESE_OVERRIDE[drug] + " ★"
        rows.append((i, drug, sol["ENGLISH_NAME"], zh, cn["cn_us"],
                     cn["cn"], cn["cn_scope"], f"{n}/{m}"))
        quick.append((drug, zh, cn["cn_clin"], cn["cn_pre"]))
        details.append((i, drug, sol, cn))

    # ---------- summary table ----------
    head = ("| 排名 | 药物 | 英文名 | 中文名 | 美国 FDA | 药理分类 | 临床用途 | "
            "AMI/MIRI |\n| ---: | --- | --- | --- | --- | --- | --- | --- |")
    lines = [head]
    for (i, drug, en, zh, us, cls, scope, score) in rows:
        def c(x, k=80):
            s = re.sub(r"\s+", " ", str(x)).strip().replace("|", "/")
            return s[:k] + ("…" if len(s) > k else "")
        lines.append(f"| {i} | **{drug}** | {c(en)} | {c(zh)} | {c(us)} | {c(cls,60)} | {c(scope,60)} | {score} |")
    summary = "\n".join(lines)

    # ---------- MIRI/AMI quick table ----------
    quick_lines = ["| 药物 | 中文名 | 临床证据 | 临床前证据 |", "| --- | --- | --- | --- |"]
    for (drug, zh, clin, pre) in quick:
        quick_lines.append(f"| {drug} | {zh} | {clin} | {pre} |")
    quick_tbl = "\n".join(quick_lines)

    # ---------- details ----------
    detail_blocks = []
    for (i, drug, sol, cn) in details:
        block = [f"### {i}. {drug}  —  AMI/MIRI 评分：{SCORES[drug][0]}/{SCORES[drug][1]}"
                 f"（{SCORE_MEANING[SCORES[drug][0]]} / {SCORE_MEANING[SCORES[drug][1]]}）"]
        # use authored Chinese where present, else translate-inline from raw solution
        for key, label in DETAIL_LABELS.items():
            v = sol.get(key, "").strip()
            if not v:
                block.append(f"- **{label}**：not found")
            elif key in ("MIRI_AMI_CLINICAL", "MIRI_AMI_PRECLINICAL"):
                authored = cn["cn_clin"] if key == "MIRI_AMI_CLINICAL" else cn["cn_pre"]
                block.append(f"- **{label}**：{authored}")
            elif key == "CLASSIFICATION":
                block.append(f"- **{label}**：{cn['cn']}")
            elif key == "CLINICAL_SCOPE":
                block.append(f"- **{label}**：{cn['cn_scope']}")
            elif key == "US_FDA_APPROVAL":
                block.append(f"- **{label}**：{cn['cn_us']}")
            else:
                block.append(f"- **{label}**：{v}")
        detail_blocks.append("\n".join(block))

    parts = [
        f"# xpert-fda-miri-v2 MIRI 药物重定位 — {n_drugs} 候选药调研报告（中文版）",
        "",
        "> 数据来源：`11_integrate/output/ranked_drugs.csv`（step 11 整合评分，integrated "
        "percentile > 0.5），每种药物由 biomni A1 agent 独立调研"
        "（KEGG/ChEMBL/openFDA/OpenTargets/PubMed/TxGNN），"
        f"{date.today().isoformat()} 生成。中文版为原报告的翻译，"
        "并据各药 MIRI/AMI 证据重新评分。",
        "",
        "## AMI/MIRI 评分说明",
        "",
        "每个药物的 AMI/MIRI 列为 `n/m`：",
        "- **n = 该药改善 AMI（急性心肌梗死）的证据强度**；**m = 该药改善 MIRI（心肌缺血再灌注损伤）的证据强度**；",
        "- 分值定义：**2** = 临床改善（有人类临床试验/结局改善证据）；**1** = 临床前改善（动物模型/机制证据）；"
        "**0** = 矛盾报道/无效（证据冲突或未见获益）；**-1** = 临床前恶化（动物模型显示有害）；"
        "**-2** = 临床恶化（临床证据显示有害，如禁忌/增加死亡率）。",
        "- 评分依据见每药详情中的临床/临床前证据；分值为主观综合判断，供候选排序参考。",
        "",
        "## 汇总表",
        "",
        summary,
        "",
        "## MIRI/AMI 关联线索速览",
        "",
        quick_tbl,
        "",
        "## 单药详情",
        "",
        "\n\n".join(detail_blocks),
        "",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fp:
        fp.write("\n".join(parts))
    print(f"Wrote {OUT} ({len(rows)} drugs)")
    print("scores:", {d: s for d, s in SCORES.items()})


if __name__ == "__main__":
    main()
