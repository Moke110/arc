# Non-ergot D2/D3 dopamine receptor agonists (pramipexole, ropinirole) and coronary heart disease — evidence-based safety assessment
# 非麦角类 D2/D3 多巴胺受体激动剂（普拉克索、罗匹尼罗）与冠心病——循证安全性评估

Prepared for the MIRI drug-repositioning screen (step-12 drug survey; candidates include **pramipexole 普拉克索** and **ropinirole 罗匹尼罗**).
Date: 2026-10-14. Citations are inline (journal + PMID, or FDA label / openFDA / FAERS / ChEMBL endpoint).
Data sources queried: PubMed/Europe PMC, DailyMed & openFDA drug-label API, openFDA FAERS (data cut `last_updated = 2026-10-13`), ChEMBL, ClinicalTrials.gov, Biomni tool suite.
为 MIRI 药物重定位筛选（step-12 药物调研；候选药含**普拉克索**与**罗匹尼罗**）准备。日期：2026-10-14。所有断言的引用均就近标注（期刊 + PMID，或 FDA 说明书 / openFDA / FAERS / ChEMBL）。检索源：PubMed/Europe PMC、DailyMed 与 openFDA 说明书 API、openFDA FAERS（数据截至 2026-10-13）、ChEMBL、ClinicalTrials.gov、Biomni 工具套件。

---

## 0. Scope note / 范围说明

**EN.** This report grades the cardiac/coronary safety of the two **non-ergot** oral dopamine agonists **pramipexole** and **ropinirole**, and separates three different questions that are easily conflated:
(i) **coronary/atherosclerotic risk** (does the drug cause MI/CHD?),
(ii) **haemodynamic/functional risk** (hypotension, syncope, falls, arrhythmia, heart failure),
(iii) **mechanistic cardiotherapeutic plausibility** (could the drug *protect* the ischaemic-reperfused myocardium?).
The single most important structural point: pramipexole and ropinirole **lack the ergoline (ergot) scaffold**, and this — not the dopamine agonism itself — is what separates them from bromocriptine/cabergoline/pergolide, whose valvulopathy and vasoconstriction are driven by **serotonin 5-HT2B** agonism (Clin Pharmacol Ther 2026, PMID:42209433; Mov Disord 2004, PMID:15197703).

**中文.** 本报告对两种**非麦角类**口服多巴胺激动剂——**普拉克索**与**罗匹尼罗**——的心脏/冠脉安全性分级，并明确区分三件容易被混淆的事：
(i)  **冠脉/动脉粥样硬化风险**（药物是否致心肌梗死/冠心病？）；
(ii) **血流动力学/功能风险**（低血压、晕厥、跌倒、心律失常、心力衰竭）；
(iii) **机制层面的心脏治疗可能性**（该药能否*保护*缺血再灌注心肌？）。
最关键的结构学要点：普拉克索与罗匹尼罗**不含麦角（ergoline）骨架**；正是这一点——而非多巴胺激动本身——把它们与溴隐亭/卡麦角林/培高利特区分开，后者的瓣膜病与血管收缩由 **5-HT2B 血清素受体**激动驱动（Clin Pharmacol Ther 2026, PMID:42209433；Mov Disord 2004, PMID:15197703）。

**Regulatory status / 监管状态.**
- **FDA.** Pramipexole (Mirapex; first approval 1997) and ropinirole (Requip; first approval 1997) are FDA-approved for Parkinson's disease and for restless legs syndrome (RLS); both are also marketed as extended-release forms (Mirapex ER, Requip XL) and as numerous generics (ChEMBL `first_approval` = 1997, CHEMBL301265 / CHEMBL589; archived brand SPLs on DailyMed: Mirapex setid `3736d691-…`; Requip XL setid `2deb758c-…`; current generic labels retrieved via openFDA label API).
- **NMPA (中国).** Both 普拉克索 (pramipexole) and 罗匹尼罗 (ropinirole) are approved and marketed in China for PD; pramipexole also for RLS. 两者均在中国获批上市用于帕金森病，普拉克索另有不宁腿综合征（RLS）适应证。

**Grading legend / 分级图例:** ✅ acceptable 可接受 · ⚠️ caution / monitoring 慎用/需监测 · ⛔ avoid 避免/禁忌 · 🎯 indicated 有适应证。

---

## 1. Pharmacology / mechanism of concern 药理学 / 机制层面的担忧

### 1.1 Receptor profile: pramipexole vs ropinirole 受体谱比较

**EN.** Both are **non-ergot, full dopamine D2-family agonists** with essentially **no D1-like (D1/D5) activity**.

- **Pramipexole** is the more **D3-preferring** agent. Aggregated Ki values in ChEMBL (assays heterogeneous, so ranges not point values): **D3 ≈ 0.7–2.8 nM**, **D2 ≈ 21–140 nM**, **D4 ≈ 2–138 nM**, **D1 > 20,000 nM** (ChEMBL ligand CHEMBL301265). It is a classic D3 > D2 >> D1 ligand (ChEMBL, CHEMBL301265); the D3-preference is standard in the literature (e.g. Expert Opin Drug Metab Toxicol 2011, PMID:21892895), and early advanced-PD trials already noted "the gastrointestinal and cardiovascular tolerability of the compound was satisfactory" (Neurology 1997, PMID:9222185).
- **Ropinirole** is a **less D3-selective, more "balanced" D2/D3** agonist: ChEMBL Ki **D3 ≈ 2.9–61 nM**, **D2 ≈ 3.7–2674 nM** (assay-dependent), **D4 ≈ 7.8 nM**, **D1A ≈ 32,500 nM**, with only weak **5-HT1A** affinity (~1.7 µM) (ChEMBL ligand CHEMBL589).
- Because both are D2-family full agonists, they **share the class haemodynamic signature** (see §2): D2-mediated inhibition of sympathetic outflow and stimulation of renal/vasodilatory dopamine receptors reduce blood pressure ("Dopamine and vascular dynamics control", Curr Neurovasc Res 2011, PMID:21722093).

**中文.** 两者均为**非麦角类、D2 家族完全激动剂**，且**几乎无 D1 样（D1/D5）活性**。**普拉克索**更偏向 **D3**：ChEMBL 汇总 Ki（各实验体系不一，故为区间）**D3 ≈ 0.7–2.8 nM**、**D2 ≈ 21–140 nM**、**D4 ≈ 2–138 nM**、**D1 > 20,000 nM**（配体 CHEMBL301265）。**罗匹尼罗**的 D3 选择性较弱、D2/D3 更"均衡"：**D3 ≈ 2.9–61 nM**、**D2 ≈ 3.7–2674 nM**、**D4 ≈ 7.8 nM**、**D1A ≈ 32,500 nM**，另仅弱亲和 **5-HT1A**（~1.7 µM）（配体 CHEMBL589）。由于均为 D2 家族完全激动剂，二者**共享类效应血流动力学特征**（见 §2）：D2 介导的交感传出抑制 + 肾/舒血管多巴胺受体刺激使血压下降（Curr Neurovasc Res 2011, PMID:21722093）。

### 1.2 Dopamine receptors in the coronary vasculature and myocardium 冠脉与心肌中的多巴胺受体

**EN.** Dopamine receptors are **present and functional in the coronary circulation**:
- **D1-like receptors mediate coronary vasodilation** — demonstrated directly in the anaesthetised dog (J Pharmacol Exp Ther 1989, PMID:2563287), and dopamine/dopamine agonists dilate both **renal and coronary** beds (Br J Pharmacol 1982, PMID:6982091).
- In the **ischaemic** heart, dopamine's coronary action is complex and **only partly dopaminergic**: in ischaemic rabbit myocardium, dopaminergic blockade modified but did not abolish dopamine's effect on coronary blood flow, with adrenergic components involved (J Cardiovasc Pharmacol 1985, PMID:2410711).
- A novel dopamine agonist produced combined **systemic, pulmonary and coronary** haemodynamic effects in awake pigs (Br J Pharmacol 1997, PMID:9134223).
- Peripherally, dopamine is an **autocrine/paracrine regulator** of vascular tone and of **renal sodium handling**, which underlies the class's BP and oedema effects (Ital J Anat Embryol 2002, PMID:12437142).
- **In the heart, D2-family (especially D3) receptors are expressed** and modulate cardiomyocyte and fibroblast biology (Front Cardiovasc Med 2021, PMID:34708087; BBRC 2020, PMID:32819568).

**Net:** dopamine agonism in the coronary bed is predominantly **vasodilatory (D1-like)**; the D2/D3 components are more relevant to **systemic BP, myocardial signalling and fibroblast/remodelling biology** than to direct coronary constriction.

**中文.** 多巴胺受体在**冠脉循环中客观存在且有功能**：**D1 样受体介导冠脉扩张**——在麻醉犬中直接证实（J Pharmacol Exp Ther 1989, PMID:2563287），多巴胺/多巴胺激动剂同时扩张**肾与冠脉**床（Br J Pharmacol 1982, PMID:6982091）。在**缺血**心肌中，多巴胺的冠脉效应复杂且**仅部分为多巴胺能**：缺血兔心肌中阻断多巴胺受体只部分改变、并未消除多巴胺对冠脉血流的作用，提示有肾上腺素能成分参与（J Cardiovasc Pharmacol 1985, PMID:2410711）。新型多巴胺激动剂在清醒猪中产生**全身、肺、冠脉**综合血流动力学效应（Br J Pharmacol 1997, PMID:9134223）。外周多巴胺是血管张力与**肾脏钠处理**的**自分泌/旁分泌调节因子**，这构成该类药降压与水肿效应的基础（Ital J Anat Embryol 2002, PMID:12437142）。**心脏中 D2 家族（尤其 D3）受体有表达**，调控心肌细胞与成纤维细胞生物学（Front Cardiovasc Med 2021, PMID:34708087；BBRC 2020, PMID:32819568）。**结论**：冠脉床上的多巴胺激动以**扩血管（D1 样）**为主；D2/D3 成分更关乎**全身血压、心肌信号与成纤维/重构生物学**，而非直接的冠脉收缩。

### 1.3 The ergot contrast — why the scaffold matters 麦角类对照——为何骨架是关键

**EN.** The cardiac toxicity of the **ergot** dopamine agonists (bromocriptine, cabergoline, pergolide, dihydroergocriptine) is **not** a dopamine-receptor phenomenon; it is a **5-HT2B-receptor phenomenon**:
- Pergolide/cabergoline cause a **fibrotic valvulopathy** clinically and echocardiographically nearly identical to fenfluramine/carcinoid valve disease, whose pathogenesis is **5-HT2B-mediated fibrogenesis in valve fibroblasts** (Mov Disord 2004, PMID:15197703). Comprehensive narrative review of the 5-HT2B literature confirms the mechanistic link to valvular heart disease (Clin Pharmacol Ther 2026, PMID:42209433).
- **The pharmacology is separable from the scaffold:** the ergot **lisuride**, which is a 5-HT2B **antagonist**, is *not* associated with cardiac valvulopathy (Clin Neuropharmacol 2006, PMID:16614540); the **bulky N6 substituent of cabergoline** confers its 5-HT2A/5-HT2B agonism (J Pharmacol Exp Ther 2011, PMID:21518772); and a 5-HT2B antagonist (cyproheptadine) **prevents pergolide-induced valvulopathy in rats** (Am J Physiol Heart Circ Physiol 2009, PMID:19346455).
- Epidemiologically, in a systematic review of observational studies, **cardiac valve regurgitation** risk was **concentrated in the ergot agents** (pergolide/cabergoline incidence-rate ratios 2.0–7.1) while **pramipexole** specifically carried a **higher HF risk (IRR 1.40–1.81)** and cabergoline an HF risk of 1.30–2.39 (CNS Drugs 2015, PMID:26585874); both are **dose/duration dependent** (Ontario administrative data, Can J Neurol Sci 2008, PMID:18574930; Neurology case-control 2006, PMID:17030757).
- Pharmacovigilance confirms the split: in **VigiBase** (WHO), heart-failure disproportionality was significant for **pergolide/cabergoline** but **absent for pramipexole, ropinirole, rotigotine and apomorphine** (Eur J Clin Pharmacol 2017, PMID:27796464). In **FAERS**, "ergot-DAs exhibited higher signal intensities for cardiac disorders" while non-ergot DAs were dominated by psychiatric signals (BMC Pharmacol Toxicol 2025, PMID:40065414).

**Critical caveat, straight from the labels:** the *pramipexole* label states that ergot-derived fibrotic complications "are believed to be related to the **ergoline structure** of these compounds, **whether other, non-ergot-derived dopamine agonists can cause them is unknown**," and that "cases of possible fibrotic complications … have been reported in the postmarketing experience with pramipexole … a contribution of pramipexole cannot be completely ruled out" (openFDA, pramipexole dihydrochloride tablets, eff. 2026-09-09; §5). The *ropinirole* label carries the identical structure of statement but names ropinirole explicitly (openFDA, ropinirole tablets, eff. 2026-09-04).

**Honest gap:** I did **not** find a retrieved primary source that reports *measured* 5-HT2B binding/agonist potency for pramipexole or ropinirole (Europe PMC title/abstract searches for `pramipexole`/`ropinirole` × `5-HT2B` returned nothing on point). The class separation therefore rests on (i) the **absence of the ergoline scaffold**, (ii) the **label's own attribution** of valvulopathy to the ergoline structure, and (iii) the **differential valvulopathy/HF epidemiology** above — not on a published 5-HT2B Ki table for these two drugs.

**中文.** **麦角类**多巴胺激动剂（溴隐亭、卡麦角林、培高利特、双氢麦角隐亭）的心脏毒性**并非**多巴胺受体现象，而是 **5-HT2B 受体现象**：培高利特/卡麦角林所致**纤维性瓣膜病**在临床与超声上与芬氟拉明/类癌瓣膜病几乎相同，后者发病机制为**瓣膜成纤维细胞中 5-HT2B 介导的纤维化**（Mov Disord 2004, PMID:15197703）；5-HT2B 文献的叙事性综述确认了这一机制联系（Clin Pharmacol Ther 2026, PMID:42209433）。**药理学可与骨架分离**：麦角类 **lisuride** 是 5-HT2B **拮抗剂**，*不*伴心脏瓣膜病（Clin Neuropharmacol 2006, PMID:16614540）；**卡麦角林庞大的 N6 取代基**赋予其 5-HT2A/5-HT2B 激动活性（J Pharmacol Exp Ther 2011, PMID:21518772）；5-HT2B 拮抗剂**赛庚啶**可在大鼠中**预防培高利特诱导的瓣膜病**（Am J Physiol Heart Circ Physiol 2009, PMID:19346455）。流行病学上，瓣膜反流风险**集中于麦角类**（培高利特/卡麦角林 IRR 2.0–7.1），且**呈剂量/疗程依赖**（Can J Neurol Sci 2008, PMID:18574930；Neurology 2006, PMID:17030757）。药物警戒亦证实此种分化：**VigiBase**（WHO）中心衰的不均衡性信号在**培高利特/卡麦角林**显著，而在**普拉克索、罗匹尼罗、罗替高汀、阿扑吗啡**中**缺失**（Eur J Clin Pharmacol 2017, PMID:27796464）；**FAERS**中"麦角类 DAs 的心脏疾病信号强度更高"，非麦角类以精神类信号为主（BMC Pharmacol Toxicol 2025, PMID:40065414）。**说明书自身的关键告诫**：*普拉克索*说明书指出麦角类纤维化并发症"被认为与这些化合物的**麦角碱结构**相关，**其他非麦角类多巴胺激动剂是否会引起尚不清楚**"，并称"在普拉克索上市后经验中曾报告可能的纤维化并发症……不能完全排除普拉克索的贡献"；*罗匹尼罗*说明书结构相同但明确点名罗匹尼罗（openFDA，见 §5）。**诚实缺口**：未检索到对普拉克索或罗匹尼罗**实测 5-HT2B 结合/激动效价**的一手文献；类别的区分依据是（i）**无麦角碱骨架**、（ii）**说明书自身**把瓣膜病归因于麦角碱结构、（iii）上述**瓣膜病/心衰流行病学差异**，而非这两个药已发表的 5-HT2B Ki 数据。

---

## 2. Hypotension and orthostatic hypotension 低血压与体位性低血压

**EN.** This is the **principal, best-established cardiovascular liability** of both drugs, and it is **labelled** in both.

**Pramipexole label (§5.2, verbatim):** *"Dopamine agonists, in clinical studies and clinical experience, appear to impair the systemic regulation of blood pressure, with resulting orthostatic hypotension, especially during dose escalation. Parkinson's disease patients, in addition, appear to have an impaired capacity to respond to an orthostatic challenge. For these reasons, both Parkinson's disease patients and RLS patients being treated with dopaminergic agonists ordinarily require careful monitoring for signs and symptoms of orthostatic hypotension, especially during dose escalation, and should be informed of this risk. In clinical trials of pramipexole, however, and despite clear orthostatic effects in normal volunteers, the reported incidence of clinically significant orthostatic hypotension was not greater among those assigned to pramipexole tablets than among those assigned to placebo. This result, especially with the higher doses used in Parkinson's disease, is clearly unexpected in light of the previous experience with the risks of dopamine agonist therapy. While this finding could reflect a unique property of pramipexole, it might also be explained by the conditions of the study and the nature of the population enrolled in the clinical trials. Patients were very carefully titrated, and **patients with active cardiovascular disease or significant orthostatic hypotension at baseline were excluded**. Also, clinical trials in patients with RLS did not incorporate orthostatic challenges with intensive blood pressure monitoring done in close temporal proximity to dosing."* (openFDA, pramipexole dihydrochloride tablets, eff. 2026-09-09.)

**Ropinirole label (§5.2 Syncope, verbatim):** *"Syncope, sometimes associated with bradycardia, was observed in association with treatment with ropinirole in both patients with Parkinson's disease and patients with RLS. In controlled clinical trials in patients with Parkinson's disease, syncope was observed more frequently in patients receiving ropinirole than in patients receiving placebo (early Parkinson's disease without levodopa [L-dopa]: ropinirole 12%, placebo 1%; advanced Parkinson's disease: ropinirole 3%, placebo 2%). Syncope was reported in 1% of patients treated with ropinirole for RLS in 12 week, placebo-controlled clinical trials compared with 0.2% of patients treated with placebo."*
**Ropinirole label (§5.3 Hypotension/Orthostatic Hypotension, verbatim):** *"… In 12 week, placebo-controlled trials of patients with RLS, the adverse event orthostatic hypotension was reported by 4 of 496 patients (0.8%) treated with ropinirole compared with 2 of 500 patients (0.4%) receiving placebo. **In 2 Phase 2 trials in patients with RLS, 14 of 55 patients (25%) receiving ropinirole experienced an adverse event of hypotension or orthostatic hypotension compared with none of the 27 patients receiving placebo.**"* (openFDA, ropinirole tablets, eff. 2026-09-04.)

**Supporting clinical evidence.**
- **Ropinirole** adverse-reaction table (early PD, ropinirole vs placebo): **syncope 12% vs 1%**, hypertension 5% vs 3%, **hypotension 2% vs 0%**, **orthostatic symptoms 6% vs 5%**, chest pain 4% vs 2%, leg oedema 7% vs 1%, dependent oedema 6% vs 3% (openFDA ropinirole label §6.1). A dedicated case report describes **ropinirole-associated orthostatic hypotension driving a prescribing cascade in an elderly man** (Cureus 2021, PMID:34268037).
- **Pramipexole** in RLS: a 3-month head-up-tilt/HRV study found pramipexole **lowered systolic BP and significantly increased tilt-induced BP variability**, i.e. measurable autonomic cardiovascular effects (Sleep Med 2015, PMID:26298781).
- **Class effect, cross-checked with another non-ergot agonist:** piribedil caused **symptomatic hypotension with bradycardia (BP 70/45 mmHg, HR 47)** at 75–100 mg, which resolved on switching to pramipexole (BMC Neurol 2018, PMID:30591018) — evidence that BP dysregulation is a **non-ergot DA class effect**, not scaffold-specific.
- In PD generally, dopamine agonists are a recognised cause of orthostatic hypotension (Drugs Aging 2019, PMID:30623310); differential effects of anti-parkinsonian regimens on cardiovascular reflex dysfunction have been documented since the pre-agonist era (Acta Neurol Scand 2004, PMID:14653849).

**Relevance for CHD patients / 对 CHD 患者的意义.** Orthostatic hypotension causes **falls, syncope and (in susceptible patients) transient hypotension that can compromise coronary perfusion** against a fixed stenosis; syncope and falls are the two commonest cardiovascular AEs in PD trials (J Parkinsons Dis 2020, PMID:32116264), and falls are a major PD morbidity (Int J Environ Res Public Health 2019, PMID:31234571). Critically, **the pramipexole development programme explicitly excluded patients with active cardiovascular disease** (label §5.2 above), so the reassuring "no excess orthostatic hypotension vs placebo" **cannot be extrapolated to a CHD cohort**.

**中文.** 这是两药**最主要、证据最充分的**心血管风险，且**两药说明书均已标注**。**普拉克索说明书 §5.2（原文）**指出：多巴胺激动剂"在临床研究与临床实践中似乎损害血压的全身调节，导致体位性低血压，尤其在剂量递增期"；PD 患者本身对体位应激反应受损；故 PD 与 RLS 患者"通常需要密切监测体位性低血压的体征与症状，尤其在剂量递增期"；但"在普拉克索临床试验中，尽管在正常志愿者中体位效应明确，被分入普拉克索组者临床显著体位性低血压的报告发生率并不高于安慰剂组"；说明书自己给出解释——"患者经过非常细致的滴定，**且基线存在活动性心血管疾病或显著体位性低血压的患者被排除**"。**罗匹尼罗说明书 §5.2（晕厥，原文）**：晕厥（有时伴心动过缓）在 PD 与 RLS 患者中均有报告；PD 对照试验中罗匹尼罗组晕厥多于安慰剂（早期 PD 不用左旋多巴：**12% vs 1%**；晚期 PD：3% vs 2%）；RLS 12 周试验中 1% vs 0.2%。**§5.3（低血压/体位性低血压，原文）**：RLS 12 周试验中体位性低血压 4/496（0.8%）vs 2/500（0.4%）；"**在 2 项 RLS 2 期试验中，55 例罗匹尼罗患者有 14 例（25%）出现低血压或体位性低血压不良事件，而 27 例安慰剂患者中为 0**"。**佐证**：罗匹尼罗不良反应表中晕厥 12% vs 1%、低血压 2% vs 0%、腿部水肿 7% vs 1%（openFDA §6.1）；个案报告示**罗匹尼罗相关体位性低血压导致老年男性用药级联**（Cureus 2021, PMID:34268037）。普拉克索在 RLS 中可**降低收缩压并显著增加直立倾斜时的血压波动**（Sleep Med 2015, PMID:26298781）。**类效应交叉验证**：另一非麦角激动剂吡贝地尔可致**症状性低血压伴心动过缓（70/45 mmHg，HR 47）**，换用普拉克索后缓解（BMC Neurol 2018, PMID:30591018）——提示血压调节障碍是**非麦角类 DAs 的类效应**，而非骨架特异。**对 CHD 的意义**：体位性低血压导致**跌倒、晕厥**，在易感者可因短暂低血压在固定狭窄下**损害冠脉灌注**；晕厥与跌倒是 PD 试验中最常见的心血管不良事件（J Parkinsons Dis 2020, PMID:32116264；Int J Environ Res Public Health 2019, PMID:31234571）。关键：**普拉克索研发项目明确排除了活动性心血管疾病患者**，故"与安慰剂相比体位性低血压无增加"这一结论**不能外推至 CHD 人群**。

---

## 3. Clinical cardiovascular safety data 临床心血管安全性数据

### 3.1 Heart failure — the most contested signal 心衰——争议最大的信号

**EN.** The evidence is **internally inconsistent**, and the inconsistency itself is the finding:

| Study | Design | Drug-specific result |
|---|---|---|
| UK GPRD nested case-control (Pharmacoepidemiol Drug Saf 2012, PMID:22109939) | 26,814 antiparkinsonian users; 783 incident HF | **Any DA RR 1.58 (1.26–1.96)**; **pramipexole RR 1.86 (1.21–2.85)**; cabergoline 2.07; **ropinirole 1.23 (0.85–1.97, ns)**; pergolide 1.42 (ns) |
| Taiwan NHIRD nested case-control (Drugs Aging 2013, PMID:23881697) | 27,135 users; 1,707 HF cases | Ergot-DAs OR 1.46 (1.00–2.12); **non-ergot OR 1.24 (0.84–1.82, ns)**; pramipexole 1.40 (ns); ropinirole 1.22 (ns) |
| Cerner Health Facts matched case-control (Mov Disord Clin Pract 2016, PMID:30363519) | 14,122 PD inpatients | **HF the only significant outcome; pramipexole AOR 1.28 (1.07–1.53)**, higher as monotherapy (1.50) |
| VigiBase disproportionality (Eur J Clin Pharmacol 2017, PMID:27796464) | 154 HF ICSRs | **Signal for pergolide/cabergoline (ergots); no signal for pramipexole, ropinirole, rotigotine, apomorphine** |
| Meta-analysis of RCTs (J Clin Med Res 2016, PMID:27222673) | 27 RCTs | **Non-ergot DAs did NOT increase incident HF (pooled RR 0.95, 0.30–2.90)**; no excess mortality (RR 0.62, ns) or CV events (RR 1.07, ns) |
| SR/MA of 44 studies (36 RCT + 8 NRS) (Brain Sci 2024, PMID:39199470) | GRADE-assessed | "Non-ergot DA use **may** be associated with adverse cardiovascular outcomes, including heart failure"; **certainty low/very low**; only RCT found no difference ropinirole vs bromocriptine (OR 0.39, 0.07–2.04) |
| Speckle-tracking echocardiography (Anatol J Cardiol 2018, PMID:30297579) | 91 PD patients | **No subclinical LV systolic dysfunction** with levodopa+pramipexole or levodopa+ropinirole (GLS −18.5/−18.5/−18.9%) |
| FAERS comparative (BMC Pharmacol Toxicol 2025, PMID:40065414) | 19.7 M reports | **Ergot-DAs higher cardiac signal intensities** than non-ergot agents |
| **This report's own FAERS analysis** (§6.2) | 20.7 M reports | **Cardiac failure ROR ≈ 1.5 for both**; see §6.2 |

Published narrative reviews reach opposite conclusions depending on which stratum they emphasise: one states the FDA "linked pramipexole use with an increased risk of heart failure (HF)" and recommends HF monitoring in patients with cardiovascular disease or oedema (Consult Pharm 2015, PMID:25760663); a second notes cabergoline's HF effect is explained by valve fibrosis but the basis for pramipexole's is unclear (Expert Opin Drug Saf 2014, PMID:24547918).

**Interpretation / 解读.** (a) The **RCT-level evidence does not support an HF signal** for non-ergot DAs (PMID:27222673). (b) The **observational HF signal is concentrated on pramipexole**, is **not confirmed for ropinirole** in the two largest case-control studies, and is **confounded by a plausible non-cardiac mechanism**: dopamine agonists cause **peripheral oedema / leg oedema** (label; Arch Neurol 2007, PMID:17420306; Postgrad Med J 2016, PMID:26941267; Intern Med J 2021, PMID:33890374) — and the Consult Pharm review itself proposes that **oedema, not direct cardiotoxicity, may be the mechanism** driving the HF association (PMID:25760663). (c) A major analytic confounder is that **PD is itself associated with substantial background cardiovascular events**: a meta-analysis of 236 PD RCTs found pooled all-CV adverse events of **3.33% (95% CI 2.14–4.70%) even in the placebo arms**, with hypertension and orthostatic hypotension the commonest (J Parkinsons Dis 2020, PMID:32116264).

**中文.** 证据**自相矛盾**，而矛盾本身就是发现：队列/病例对照研究（英国 GPRD，Pharmacoepidemiol Drug Saf 2012, PMID:22109939）显示任意 DA **RR 1.58**、**普拉克索 RR 1.86**（显著）、卡麦角林 2.07，而**罗匹尼罗 1.23（不显著）**；台湾 NHIRD（Drugs Aging 2013, PMID:23881697）非麦角类 OR 1.24（不显著）；Cerner 数据库（Mov Disord Clin Pract 2016, PMID:30363519）**仅心衰达显著，普拉克索 AOR 1.28**；而 **VigiBase 的不均衡性分析中普拉克索/罗匹尼罗均无信号**，仅麦角类有信号（Eur J Clin Pharmacol 2017, PMID:27796464）。**RCT 荟萃分析（27 项）显示非麦角类不增加新发心衰（RR 0.95）**，死亡与心血管事件亦无增加（J Clin Med Res 2016, PMID:27222673）；44 项研究的 SR/MA 认为"可能相关"但**证据确定性低/极低**（Brain Sci 2024, PMID:39199470）；斑点追踪超声未见亚临床左室收缩功能障碍（Anatol J Cardiol 2018, PMID:30297579）。**解读**：(a) **RCT 层面不支持**非麦角类的心衰信号；(b) 观察性心衰信号**集中于普拉克索**，在两项最大病例对照研究中**罗匹尼罗未被证实**，且存在**合理的非心脏混杂机制**——多巴胺激动剂致**外周/腿部水肿**（说明书；Arch Neurol 2007, PMID:17420306；Postgrad Med J 2016, PMID:26941267；Intern Med J 2021, PMID:33890374），有综述即认为**水肿而非直接心脏毒性**才是心衰关联的机制（Consult Pharm 2015, PMID:25760663）；(c) 重要混杂：**PD 本身背景心血管事件率高**——236 项 PD RCT 的 meta 分析显示**安慰剂组**全因心血管不良事件即达 **3.33%**（J Parkinsons Dis 2020, PMID:32116264）。

### 3.2 Myocardial infarction and coronary endpoints 心肌梗死与冠脉终点

**EN.** **No evidence found** of an increased **MI** risk for either drug in any controlled human study retrieved. Specifically:
- In the Cerner matched case-control analysis (Mov Disord Clin Pract 2016, PMID:30363519), **acute myocardial infarction was among the outcomes tested and was *not* significantly associated** with non-ergot DA use — HF was the only positive outcome.
- No RCT or cohort retrieved reported MI as a drug-attributable signal for pramipexole or ropinirole.
- In **RLS**, treatment with (predominantly **non-ergot**) dopamine agonists was associated with **LOWER** future CVD risk than untreated RLS — adjusted HR for CVD 1.26 (1.20–1.32) in all RLS vs non-RLS, but **significantly lower in treated vs untreated RLS**; **notably, "except for ergot-dopamine use"** (J Am Heart Assoc 2021, PMID:33550813). This is the single most reassuring human dataset for the non-ergot class and directly addresses MI/angina/stroke/AF/HF as a **composite**.
- RLS itself carries elevated CVD prevalence (JAMA 2026 review: 29.6% with CAD, stroke or HF; PMID:41563785), and cohort studies of RLS and incident CVD are mixed (no association in Women's Health Study/Physicians' Health Study, BMJ Open 2012, PMID:22447047; positive in Kaiser Permanente, Sleep 2015, PMID:26083613; positive signal limited to severe/long-duration RLS, Sleep 2017 systematic review, PMID:28329227) — i.e. **disease-related, not necessarily drug-related**, risk.

**中文.** **未找到**任何受控人体研究提示这两个药增加**心梗**风险。具体：Cerner 病例对照分析（Mov Disord Clin Pract 2016, PMID:30363519）**把急性心梗列为检验结局且未达显著**，唯一阳性结局是心衰。未检索到以 MI 为药物归因信号的 RCT 或队列。在 **RLS** 中，以**非麦角类**为主的 DA 治疗与**更低的**未来 CVD 风险相关——全部 RLS vs 非 RLS 的调整 HR 为 1.26（1.20–1.32），但**接受治疗者显著低于未治疗者**，**唯"麦角类多巴胺"使用者例外**（J Am Heart Assoc 2021, PMID:33550813）；这是非麦角类最有力的人体安全性数据，且终点涵盖 MI/心绞痛/卒中/房颤/心衰的**复合终点**。RLS 本身 CVD 患病率高（JAMA 2026 综述：29.6% 合并 CAD、卒中或心衰；PMID:41563785），而 RLS 与 CVD 的队列研究结果不一（BMJ Open 2012, PMID:22447047 无关联；Sleep 2015, PMID:26083613 有关联；Sleep 2017 系统综述 PMID:28329227 提示仅重度/长病程 RLS 有关联）——即风险**与疾病相关，而非必然与药物相关**。

### 3.3 Hypertension / blood pressure 高血压/血压

**EN.** Both labels list hypertension among adverse reactions; **neither labels hypertension as a warning-level class effect**. Ropinirole: **"Hypertension 5% vs 3%" (placebo)** in early PD, and in a fixed-dose advanced-PD ER study hypertension was "reported as an adverse reaction in 3% of patients on ropinirole extended-release tablets, compared with 1% of patients on placebo" (openFDA ropinirole label §5/§6); the same label notes ropinirole ER had **no significant effect on average heart rate** (openFDA, ropinirole ER). Pramipexole's older full postmarketing listing includes **hypertension, tachycardia and chest pain** among "Cardiovascular disorders" (openFDA, PRAMIPEXOLE label, eff. 2016-10-13). Net: **BP effects are bidirectional and modest at therapeutic doses** — hypotension/orthostatic hypotension is the dominant labelled liability, with a smaller hypertension signal.

**中文.** 两药说明书均将高血压列入不良反应，但**均未把高血压列为"警示"级类效应**。罗匹尼罗：早期 PD **高血压 5% vs 安慰剂 3%**；晚期 PD 固定剂量 ER 研究中"罗匹尼罗缓释片组高血压 3% vs 安慰剂 1%"（openFDA 罗匹尼罗说明书 §5/§6）；同说明书指出罗匹尼罗缓释片**对平均心率无显著影响**。普拉克索旧版完整上市后清单在"心血管系统疾病"中列出**高血压、心动过速、胸痛**（openFDA，PRAMIPEXOLE 说明书，eff. 2016-10-13）。总体：**治疗剂量下血压效应为双向且幅度有限**——低血压/体位性低血压是主导的标示风险，另有较小的升压信号。

### 3.4 RLS-specific evidence 不宁腿综合征特异证据

**EN.** (a) **Augmentation:** dopamine agonists cause **augmentation** (iatrogenic worsening of RLS), annual incidence **7–10%** with ropinirole/pramipexole/rotigotine, which is precisely why DAs are **no longer first-line** for RLS (JAMA 2026 review, PMID:41563785; Mov Disord 2018 evidence-based review, PMID:29756335). Augmentation/early-morning rebound is a labelled warning for ropinirole (§5.9, verbatim in §5 below). (b) **Impulse-control disorders (ICD):** a labelled warning for both agents (§5.3 pramipexole / §5.6 ropinirole); in FAERS, **pramipexole showed the most significant ICD signal of the six DAs** (BMC Pharmacol Toxicol 2025, PMID:40065414), and ICDs differ by DA (Clin Park Relat Disord 2024, PMID:38292814). (c) **CV signals in RLS:** the class CV effect in RLS is the **high trial-level incidence of hypotension/orthostatic hypotension** seen in the label's Phase-2 RLS data (25% vs 0%, §2) plus measurable tilt-BP changes (Sleep Med 2015, PMID:26298781); by contrast, **treated RLS had *lower* incident CVD than untreated RLS** (J Am Heart Assoc 2021, PMID:33550813). Chinese RLS RCTs with pramipexole reported it "effective and well-tolerated" (Sleep Med 2012, PMID:22137119; Sleep Med 2015, PMID:25533540).

**中文.** (a) **加重（augmentation）**：多巴胺激动剂致 RLS 症状出现医源性恶化，**年发生率 7–10%**，这正是 DAs **不再是一线**的原因（JAMA 2026 综述，PMID:41563785；Mov Disord 2018 循证综述，PMID:29756335）；加重/清晨反跳是罗匹尼罗的标示警示（§5.9）。 (b) **冲动控制障碍（ICD）**：两药均标示（§5.3 普拉克索 / §5.6 罗匹尼罗）；FAERS 中**普拉克索在六种 DA 中 ICD 信号最强**（BMC Pharmacol Toxicol 2025, PMID:40065414）。 (c) **RLS 中的心血管信号**：主要来自说明书 2 期 RLS 数据中**低血压/体位性低血压的试验级高发生率（25% vs 0%，§2）**及可测的倾斜血压变化（Sleep Med 2015, PMID:26298781）；与之相对，**接受治疗的 RLS 未来 CVD 风险低于未治疗者**（J Am Heart Assoc 2021, PMID:33550813）。中国 RLS 普拉克索 RCT 报告"有效且耐受良好"（Sleep Med 2012, PMID:22137119；Sleep Med 2015, PMID:25533540）。

---

## 4. Pleiotropic / mechanistic cardiac effects 多效性 / 机制性心脏效应

**EN.** This is where the **MIRI repositioning hypothesis has genuine (preclinical) support** — and it is important to label each claim as *observed* vs *speculative*.

### 4.1 Observed experimental evidence — cardioprotection in ischaemia-reperfusion 实验观察——缺血再灌注心肌保护

1. **Pramipexole, direct MIRI model (highest relevance).** In an **in-vivo mouse myocardial I/R** model plus an **H9c2 hypoxia/reoxygenation** model, **pramipexole pre-treatment significantly reduced myocardial infarct size**, decreased serum CK and LDH, **reduced ROS**, **restored mitochondrial membrane potential (ΔΨm)**, increased **survival of H9c2 cells**, and increased **autophagic flux**, acting at least partly via an **AMPK-dependent** pathway (Biochem Biophys Res Commun 2016, PMID:27063800). *This is the single most directly relevant paper to a MIRI drug-repositioning screen.*
2. **D2 receptors and ischaemic post-conditioning.** Dopamine D2 receptors contribute to **ischaemic post-conditioning cardioprotection via autophagy in isolated rat hearts** (Int J Cardiol 2016, PMID:26599746); D2 activation mediates post-conditioning cardioprotection **via PKC-ε translocation** (Mol Cell Biochem 2013, PMID:23636617) and **post-conditioning-attenuated cardiomyocyte apoptosis** (Exp Cell Res 2014, PMID:24412422); **preconditioning cardioprotection in hyperlipidaemic rat hearts involves D2 signalling** (BMC Cardiovasc Disord 2015, PMID:26216101). Note: these four use **bromocriptine** (an ergot) as the D2 agonist, so they establish *receptor* biology, not the non-ergot drugs' efficacy.
3. **Rotigotine (non-ergot D2/D3 agonist) in a mouse myocardial-ischaemia model.** Rotigotine and its 14-day extended-release microsphere formulation **reduced troponin I and CK-MB and alleviated ischaemia-induced cardiomyocyte histopathology**; the effect was blocked by the dopamine-receptor blocker **chlorpromazine**, implicating dopamine receptors; protection persisted 14 days with the depot formulation (Eur J Pharm Sci 2019, PMID:30959104). *This is a second, independent non-ergot-DA cardioprotection dataset.*
4. **DRD2 in MIRI.** In a mouse MIRI model, **DRD2 acted as an anti-inflammatory factor suppressing NF-κB-dependent inflammation and preventing MIRI** (via the RUNX1/miR-142-3p/DRD2 axis) (Exp Ther Med 2022, PMID:36561980).
5. **Cardiac-fibrosis / remodelling biology.** The D3 agonist **pramipexole abolished morphine-induced interstitial fibrosis and cardiomyocyte hypertrophy** in mice (BBRC 2020, PMID:32819568). However — a **directionally opposite** finding: in mouse LV cardiac fibroblasts, **pramipexole (D3 agonist) *enhanced* fibroblast migration after scratch injury**, and D3 loss-of-function attenuated proliferation/migration and raised Col3a1 (Front Cardiovasc Med 2021, PMID:34708087). So **D3 modulation of cardiac fibroblasts is bidirectional/context-dependent** — a caution against assuming pure anti-fibrotic benefit.
6. **Autophagy is a plausible shared mechanism.** **DRD3 (but not DRD2) activates autophagy through MTORC1 inhibition** via AMPK–RPTOR–ULK1 (Autophagy 2020, PMID:31538542); D3-mediated autophagy is neuroprotective (Cells 2025, PMID:40358175); pramipexole inhibits astrocytic NLRP3 via Drd3-dependent autophagy (Acta Pharmacol Sin 2023, PMID:35896696). Since the MIRI-protective effect of pramipexole was **autophagy-dependent (AMPK)** (PMID:27063800), a **D3→AMPK→autophagy** axis is a coherent, testable mechanistic bridge for MIRI — but this bridge has **not** been demonstrated directly in myocardium for pramipexole.
7. **D2/D3 neuroprotection in ischaemia** (supportive class-level biology): D2/D3 activation participates in **ischaemic preconditioning-mediated brain protection** (Brain Res 2020, PMID:32919985); the D2/D3 agonist **piribedil** protects against global cerebral ischaemia in rats (Neurosci Lett 2018, PMID:30107200); the D2/D3 agonist **pramipexole** protected against post-ischaemic nigrostriatal degeneration as early as 1996 (Brain Res 1996, PMID:9117424).

### 4.2 Antioxidant / mitochondrial (MPTP-model) mechanism 抗氧化 / 线粒体（MPTP 模型）机制

**EN.** Pramipexole has a **well-documented, receptor-independent antioxidant/mitochondrial action**:
- **Pramipexole inhibits lipid peroxidation and reduces MPP⁺-induced substantia-nigra injury** (Neurosci Lett 2000, PMID:10704769) and protects against dopamine/levodopa-induced cytotoxicity (Life Sci 1999, PMID:10227583); it protects against **MPTP toxicity in non-human primates** (J Neurochem 2006, PMID:16464239).
- **Cellular/mitochondrial basis:** pramipexole prevents ischaemic cell death via **mitochondrial pathways** (Dis Model Mech 2019, PMID:31235613); it preserves **mitochondrial membrane potential** after cerebral I/R (J Stroke Cereb Dis 2023, PMID:37105127); it suppresses **ferroptosis** after cerebral I/R via the **Nrf2/SLC7A11/GPX4** axis (J Stroke Cereb Dis 2025, PMID:39490461; Sci Rep 2026, PMID:41741591); overview of mitochondrial bioenergetics and dopamine agonists (Mol Neurobiol 2024, PMID:38468113).
- **Enantiomer separation supports the antioxidant-only story:** **dexpramipexole** is the R(+) enantiomer with minimal dopamine-agonist activity; in the **EXHALE randomised trial** in eosinophilic asthma it was safe with **no cardiac safety signal** (J Allergy Clin Immunol 2023, PMID:37277072), yet it shares pramipexole's antinociceptive/mitochondrial pharmacology (Neurosci Lett 2023, PMID:37558176) — evidence that the mitochondrial/antioxidant effects are **separable from D2/D3 agonism**.
- **Systemic anti-inflammatory / antioxidant effects across organs:** pramipexole reduces inflammation in experimental inflammation models (Immunopharmacol Immunotoxicol 2017, PMID:28162057), attenuates rotenone-induced **hepatic** injury (Int Immunopharmacol 2026, PMID:42391877), and halts leflunomide **nephrotoxicity** (Eur J Pharmacol 2025, PMID:41109403); neuroprotection via ferroptosis inhibition also reported (Sci Rep 2026, PMID:41741591).

### 4.3 Human translational caveat — dopamine agonism is NOT an established human cardioprotectant 人体转化警告——多巴胺激动剂并非已证实的人体心脏保护剂

**EN.** The D1-agonist **fenoldopam** has been tested for organ protection in **human cardiac surgery**: a meta-analysis of 7 RCTs (1,107 patients) found a **reduction in AKI incidence** (RR 0.42, 0.26–0.69) (Shock 2019, PMID:30601331), but the **Cochrane 2024 review** ("Fenoldopam for preventing and treating acute kidney injury") does not establish a mortality/outcome benefit (Cochrane Database Syst Rev 2024, PMID:39607014), and nephroprotective claims in cardiac surgery remain disputed (J Clin Med 2024, PMID:39407923; P R Health Sci J 2025, PMID:40163743). Likewise, **"renal-dose" dopamine in critically ill patients is NOT protective** (Interact Cardiovasc Thorac Surg 2004, PMID:17670192; Curr Surg 2006, PMID:16757377). **Conclusion:** dopaminergic agonism has repeated human *organ-protection* failures — so the preclinical MIRI data in §4.1 must be treated as **hypothesis-generating, not clinical precedent**.

### 4.4 Mechanistic speculation, stated as such 明确的机制推测

**EN.** *Speculative (not directly demonstrated):* that pramipexole's **AMPK/autophagy** myocardial protection (PMID:27063800) is the same D3→AMPK→autophagy axis characterised in non-cardiac cells (PMID:31538542, PMID:40358175); that **D3-preferring** selectivity confers an advantage over ropinirole for MIRI (no head-to-head myocardial comparison exists); that antioxidant/mitochondrial actions are dose-achievable at clinically tolerable doses without hypotension. *Contradictory:* the pro-migratory effect of pramipexole on cardiac fibroblasts (PMID:34708087) sits against its anti-fibrotic effect (PMID:32819568). **No** retrieved study has tested pramipexole or ropinirole in a **human** ischaemia-reperfusion or post-MI setting.

**中文.** 这是 **MIRI 重定位假设唯一有真正（临床前）支持**的地方，务必逐条区分*实测*与*推测*。**实测**：1）**普拉克索直接 MIRI 模型**——小鼠在体心肌 I/R + H9c2 缺氧/复氧模型中，普拉克索预处理**显著缩小梗死面积**、降低 CK/LDH、**减少 ROS**、**恢复线粒体膜电位**、提高细胞存活、增加**自噬流**，机制至少部分经 **AMPK**（Biochem Biophys Res Commun 2016, PMID:27063800）——这是与 MIRI 筛选最直接相关的一篇。2）**D2 受体与缺血后处理**：D2 受体经**自噬**参与离体大鼠心脏缺血后处理保护（Int J Cardiol 2016, PMID:26599746），经 **PKC-ε 转位**（Mol Cell Biochem 2013, PMID:23636617），并**减轻后处理抗心肌细胞凋亡**（Exp Cell Res 2014, PMID:24412422）；高脂大鼠心脏预适应保护涉及 D2 信号（BMC Cardiovasc Disord 2015, PMID:26216101）——注意这四项以**溴隐亭（麦角类）**为 D2 激动剂，界定的是*受体*生物学而非非麦角药的疗效。3）**罗替高汀（非麦角类 D2/D3）**在小鼠心肌缺血模型中**降低肌钙蛋白 I 与 CK-MB、减轻心肌细胞病理**，被多巴胺受体阻断剂**氯丙嗪**阻断，缓释微球保护可持续 14 天（Eur J Pharm Sci 2019, PMID:30959104）——第二个独立的非麦角药心脏保护数据集。4）**DRD2**在小鼠 MIRI 中作为**抗炎因子抑制 NF-κB 依赖炎症**（Exp Ther Med 2022, PMID:36561980）。5）**心脏纤维化/重构**：D3 激动剂**普拉克索消除吗啡诱导的心肌间质纤维化与心肌细胞肥大**（BBRC 2020, PMID:32819568）；但**方向相反**的结果——在小鼠左室成纤维细胞中，**普拉克索（D3 激动）反而增强划痕后迁移**，而 D3 功能缺失减弱增殖/迁移（Front Cardiovasc Med 2021, PMID:34708087）：**D3 对心脏成纤维细胞的调控是双向/依赖语境的**，不可假定为纯粹抗纤维化获益。6）**自噬是合理的共享机制**：**DRD3（而非 DRD2）经 AMPK–MTORC1–ULK1 激活自噬**（Autophagy 2020, PMID:31538542），D3 介导的自噬具神经保护（Cells 2025, PMID:40358175）；普拉克索经 Drd3 依赖自噬抑制星形胶质 NLRP3（Acta Pharmacol Sin 2023, PMID:35896696）。**抗氧化/线粒体（MPTP 模型）**：普拉克索有**不依赖受体的抗氧化/线粒体作用**——抑制脂质过氧化、减轻 MPP⁺ 所致黑质损伤（Neurosci Lett 2000, PMID:10704769）、抗多巴胺/左旋多巴细胞毒（Life Sci 1999, PMID:10227583）、在非人灵长类抗 MPTP 毒性（J Neurochem 2006, PMID:16464239）；普拉克索经**线粒体通路**防止缺血性细胞死亡（Dis Model Mech 2019, PMID:31235613）、维持脑 I/R 后**线粒体膜电位**（J Stroke Cereb Dis 2023, PMID:37105127）、经 **Nrf2/SLC7A11/GPX4** 抑制**铁死亡**（J Stroke Cereb Dis 2025, PMID:39490461；Sci Rep 2026, PMID:41741591）。**对映体分离支持"抗氧化独立于 DA 激动"**：**右旋普拉克索（dexpramipexole）**DA 激动活性极弱，**EXHALE 随机试验**中安全且**无心脏安全性信号**（J Allergy Clin Immunol 2023, PMID:37277072），却共享普拉克索的镇痛/线粒体药理（Neurosci Lett 2023, PMID:37558176）。**转化警告**：D1 激动剂**非诺多泮**在人心脏手术中虽有 AKI 减少的 meta 分析结果（Shock 2019, PMID:30601331），但 **Cochrane 2024** 未确立临床结局获益（PMID:39607014），且"肾脏剂量"多巴胺在危重症中**并无保护作用**（PMID:17670192；PMID:16757377）——故 §4.1 的临床前 MIRI 数据应视为**产生假设、而非临床先例**。**明确的机制推测**：普拉克索的**AMPK/自噬**心肌保护（PMID:27063800）与在非心脏细胞中界定的 **D3→AMPK→自噬**轴（PMID:31538542、PMID:40358175）是同一通路（*未经心肌直接证实*）；**D3 高选择性**相对罗匹尼罗对 MIRI 有优势（**无**心肌头对头比较）；抗氧化/线粒体效应可在临床耐受剂量下不伴低血压地实现（*未证实*）。**无**任何检索到的研究在**人体**缺血再灌注或梗死后情境中测试过普拉克索或罗匹尼罗。

---

## 5. Label warnings (verbatim) 说明书警示（原文）

### 5.1 Pramipexole (Mirapex / Mirapex ER / generics) 普拉克索

**EN.** Source: openFDA drug-label API, "Pramipexole Dihydrochloride" tablets, effective 2026-09-09 (37 pramipexole labels inspected, incl. ER forms); brand Mirapex SPL archived on DailyMed (setid `3736d691-e8b9-4fdd-aeca-9e96cf907c51`).

- **§4 Contraindications (verbatim):** *"None."*
- **§5.1 Falling Asleep During Activities of Daily Living and Somnolence** — "Sudden onset of sleep may occur without warning."
- **§5.2 Symptomatic Orthostatic Hypotension** — full text quoted in §2 above, including the explicit statement that **patients with active cardiovascular disease were excluded** from the trials and that **careful monitoring during dose escalation** is required.
- **§5.3 Impulse Control / Compulsive Behaviors (verbatim):** *"Case reports and the results of a cross-sectional study suggest that patients can experience intense urges to gamble, increased sexual urges, intense urges to spend money uncontrollably, binge eating, and/or other intense urges and the inability to control these urges while taking one or more of the medications, including pramipexole dihydrochloride, that increase central dopaminergic tone. In some cases, although not all, these urges were reported to have stopped when the dose was reduced or the medication was discontinued."*
- **§5.4 Hallucinations and Psychotic-like Behavior** (risk increases with age); **§5.5 Dyskinesia**; **§5.6 Postural Deformity**; **§5.8 Rhabdomyolysis**; **§5.9 Retinal Pathology**; **§5.11 Withdrawal Symptoms**.
- **§5.10 Events Reported with Dopaminergic Therapy — Fibrotic Complications (verbatim):** *"Cases of retroperitoneal fibrosis, pulmonary infiltrates, pleural effusion, pleural thickening, pericarditis, and cardiac valvulopathy have been reported in patients treated with ergot-derived dopaminergic agents. … Although these adverse events are believed to be related to the ergoline structure of these compounds, whether other, non-ergot-derived dopamine agonists can cause them is unknown. Cases of possible fibrotic complications, including peritoneal fibrosis, pleural fibrosis, and pulmonary fibrosis have been reported in the postmarketing experience with pramipexole dihydrochloride tablets. While the evidence is not sufficient to establish a causal relationship … a contribution of pramipexole dihydrochloride tablets cannot be completely ruled out."*
- **Pharmacovigilance / postmarketing adverse reactions:** current label lists **"Cardiac Disorders: cardiac failure"** and **"Nervous System Disorders: syncope"** among postmarketing events (openFDA). The **2016 free-text listing** (PRAMIPEXOLE, eff. 2016-10-13) additionally enumerated, under "Cardiac disorders": *"… atrioventricular block second degree, bradycardia, bundle branch block, cardiac arrest, cardiac failure, cardiac failure congestive, cardiomegaly, coronary artery occlusion, cyanosis, extrasystoles, left ventricular failure, myocardial infarction, nodal arrhythmia, sinus arrhythmia, sinus bradycardia, sinus tachycardia, supraventricular extrasystoles, supraventricular tachycardia, tachycardia, ventricular fibrillation, ventricular extrasystoles, ventricular hypertrophy"* — plus hypertension and chest pain elsewhere in the Cardiovascular list.
- **§7 Drug Interactions (verbatim):** *"Dopamine antagonists: May diminish the effectiveness of pramipexole (7.1)."* — the **only** drug-interaction entry. **No MAOI interaction or contraindication is listed** (the label's clinical-trial text notes that **MAO-B inhibitors** were permitted concomitant therapy).
- **No boxed warning.**

**中文.** 来源：openFDA 说明书 API，"Pramipexole Dihydrochloride"片，生效 2026-09-09；品牌 Mirapex 的 SPL 存于 DailyMed。**§4 禁忌（原文）："None."（无）**。**§5.1** 日间活动突睡可用；**§5.2 症状性体位性低血压**（全文见 §2，含"**排除了活动性心血管疾病及基线显著体位性低血压患者**"与剂量递增期需密切监测）；**§5.3 冲动控制/强迫行为**（原文见上）；**§5.4** 幻觉/精神病样行为（随年龄增加）；**§5.5** 异动症；**§5.6** 姿势畸形；**§5.8** 横纹肌溶解；**§5.9** 视网膜病变；**§5.11** 停药症状。**§5.10 纤维化并发症**：麦角类致腹膜后纤维化、胸腔积液、心包炎、**心脏瓣膜病**，这些"被认为与**麦角碱结构**相关，**其他非麦角类多巴胺激动剂是否会引起尚不清楚**"，且普拉克索上市后有"可能的纤维化并发症"报告、"**不能完全排除**普拉克索的贡献"。**上市后不良反应**：现行说明书列出 **"心脏疾病：心力衰竭"** 与 **"神经系统：晕厥"**；**2016 年自由文本清单**另列出**房室传导阻滞、心动过缓、心脏骤停、心力衰竭、充血性心衰、心脏肥大、冠脉闭塞、发绀、早搏、左室衰竭、心肌梗死、窦性心律不齐/心动过缓/心动过速、室上性心动过速、心动过速、室颤、室性早搏、心室肥厚**，以及高血压、胸痛。**§7 药物相互作用**仅一条：**多巴胺拮抗剂**（可减弱疗效）；**未列出任何 MAOI 相互作用或禁忌**（说明书临床试验部分反而说明**MAO-B 抑制剂**为允许的合并用药）。**无黑框警告**。

### 5.2 Ropinirole (Requip / Requip XL / generics) 罗匹尼罗

**EN.** Source: openFDA drug-label API, ropinirole tablets, effective 2026-09-04 (45 ropinirole labels inspected, incl. ER); Requip XL SPL archived on DailyMed (setid `2deb758c-659e-4a87-…`).

- **§4 Contraindications (verbatim):** *"Ropinirole tablets are contraindicated in patients known to have a hypersensitivity/allergic reaction (including urticaria, angioedema, rash, pruritus) to ropinirole or to any of the excipients."* → **no cardiac contraindication.**
- **§5 Warnings headline list (verbatim):** *"• Sudden onset of sleep and somnolence may occur (5.1) • **Syncope may occur (5.2)** • **Hypotension, including orthostatic hypotension may occur (5.3)** • May cause hallucinations and psychotic-like behaviors (5.4) • May cause or exacerbate dyskinesia (5.5) • May cause problems with impulse control or compulsive behaviors (5.6)."*
- **§5.2 Syncope** and **§5.3 Hypotension/Orthostatic Hypotension** — verbatim quantitative text quoted in §2 above (early-PD syncope 12% vs 1%; advanced PD 3% vs 2%; RLS 1% vs 0.2%; **RLS Phase-2 hypotension/orthostatic hypotension 25% vs 0%**).
- **§5.4 Hallucinations/psychotic-like behavior; §5.5 Dyskinesia; §5.6 Impulse control/compulsive behaviors (verbatim in §2/§3.4); §5.7 Withdrawal-emergent hyperpyrexia and confusion; §5.8 Withdrawal Symptoms; §5.9 Augmentation and early-morning rebound in RLS (verbatim §2); §5.10 Fibrotic complications; §5.11 Retinal pathology.**
- **§5.10 Fibrotic Complications (verbatim):** *"Cases of retroperitoneal fibrosis, pulmonary infiltrates, pleural effusion, pleural thickening, pericarditis, and cardiac valvulopathy have been reported in some patients treated with ergot-derived dopaminergic agents. … Although these adverse reactions are believed to be related to the ergoline structure of these compounds, whether other, non-ergot-derived dopamine agonists **such as ropinirole** can cause them is unknown. Cases of possible fibrotic complications, including pleural effusion, pleural fibrosis, interstitial lung disease, and cardiac valvulopathy have been reported in the development program and postmarketing experience for ropinirole. While the evidence is not sufficient to establish a causal relationship …"*
- **§6.1 Adverse reactions (PD, ropinirole vs placebo):** syncope 12% vs 1%; hypertension 5% vs 3%; hypotension 2% vs 0%; orthostatic symptoms 6% vs 5%; chest pain 4% vs 2%; leg oedema 7% vs 1%; dependent oedema 6% vs 3%. Ropinirole ER: **"no significant effect … on average heart rate"**; in a fixed-dose advanced-PD study, hypertension 3% vs 1% (openFDA).
- **§7 Drug Interactions (verbatim):** *"Inhibitors or inducers of CYP1A2: May alter the clearance of ropinirole …; Hormone replacement therapy (HRT): …; Dopamine antagonists (e.g., neuroleptics, metoclopramide): May reduce efficacy of ropinirole."* **No MAOI interaction or contraindication is listed.**
- **No boxed warning.**

**中文.** 来源：openFDA 说明书 API，罗匹尼罗片，生效 2026-09-04；Requip XL 的 SPL 存于 DailyMed。**§4 禁忌（原文）**：已知对罗匹尼罗或辅料过敏（荨麻疹、血管性水肿、皮疹、瘙痒）——**无心脏禁忌**。**§5 警示**：突睡/嗜睡；**可发生晕厥**；**可发生低血压（含体位性低血压）**；幻觉/精神病样行为；可致或加重异动症；冲动控制/强迫行为。**§5.2/§5.3 定量原文见 §2**（早期 PD 晕厥 12% vs 1%；晚期 3% vs 2%；RLS 1% vs 0.2%；**RLS 2 期低血压/体位性低血压 25% vs 0%**）。**§5.10 纤维化并发症**：麦角碱结构相关，**"其他非麦角类多巴胺激动剂如罗匹尼罗是否会引起尚不清楚"**，但罗匹尼罗研发与上市后经验中曾报告胸腔积液、胸膜纤维化、间质性肺病与**心脏瓣膜病**，"证据不足以确立因果关系"。**§6.1 不良反应（PD，罗匹尼罗 vs 安慰剂）**：晕厥 12% vs 1%、高血压 5% vs 3%、低血压 2% vs 0%、体位性症状 6% vs 5%、胸痛 4% vs 2%、腿部水肿 7% vs 1%。缓释剂**对平均心率无显著影响**。**§7 相互作用**：CYP1A2 抑制剂/诱导剂、激素替代治疗、多巴胺拮抗剂；**未列出 MAOI 相互作用或禁忌**。**无黑框警告。**

### 5.3 MAOI interaction — the honest answer MAOI 相互作用——诚实的回答

**EN.** **Neither the pramipexole nor the ropinirole US label lists an MAOI interaction or contraindication** (openFDA, both labels, §7 above). For pramipexole only dopamine antagonists appear; for ropinirole only CYP1A2/HRT/dopamine antagonists. The pramipexole label's trial description even records that **MAO-B inhibitors were permitted** co-medication. Therefore the "MAOI interaction" sometimes attributed to dopaminergic drugs is **not a labelled feature of these two agents**; the classic MAOI-hypertensive ("cheese") concern belongs to **MAO-A inhibition** and to older ergot agents, and in PD practice the relevant co-class is **MAO-B inhibitors** (selegiline/rasagiline), whose own labels carry the tyramine/interaction cautions.

**中文.** **普拉克索与罗匹尼罗的美国说明书均未列入 MAOI 相互作用或禁忌**。普拉克索仅列多巴胺拮抗剂；罗匹尼罗仅列 CYP1A2/HRT/多巴胺拮抗剂。普拉克索说明书甚至记载**MAO-B 抑制剂为允许合并用药**。因此常被归于"多巴胺能药物"的"MAOI 相互作用"**并非这两药的标示特征**；经典的 MAOI 高血压（"奶酪"）风险属于 **MAO-A 抑制**及老一代麦角类；在 PD 实践中相关共类为 **MAO-B 抑制剂**（司来吉兰/雷沙吉兰），其自身说明书含酪胺/相互作用警示。

---

## 6. Reported MI / chest pain / arrhythmia signals 已报告的心梗/胸痛/心律失常信号

### 6.1 Case reports and label pharmacovigilance 病例报告与说明书药物警戒

**EN.**
- **Reversible heart failure attributed to pramipexole** — dedicated case report (Mov Disord Clin Pract 2014, PMID:30363854); a prescribing alert also summarised "Pramipexole: heart failure" (Prescrire Int 2013, PMID:24171217).
- **Peripheral/leg oedema with pramipexole and other dopaminergic drugs** — multiple reports and a risk-factor study (Postgrad Med J 2016, PMID:26941267; Rev Med Chil 2017, PMID:29652961; Intern Med J 2021, PMID:33890374; Arch Neurol 2007, PMID:17420306). Oedema is relevant because it is the mechanism most often invoked to explain the HF association (Consult Pharm 2015, PMID:25760663).
- **Syncope / hypotension with ropinirole** (Cureus 2021, PMID:34268037); **hypotension + bradycardia with another non-ergot DA, piribedil** (BMC Neurol 2018, PMID:30591018); a **fatal ropinirole** toxicology case has been published (Forensic Toxicol 2022, PMID:36454487).
- **Label-level MI/arrhythmia terms:** pramipexole's historical postmarketing list includes **myocardial infarction, coronary artery occlusion, ventricular fibrillation, ventricular extrasystoles, tachycardia, cardiac arrest, left ventricular failure** (openFDA, eff. 2016-10-13); the current label lists **cardiac failure** and **syncope** (openFDA). Ropinirole's label lists **chest pain 4% vs 2%** and **leg/dependent oedema** (openFDA §6.1).
- **No dedicated case report or label term for acute MI attributable to ropinirole** was retrieved.

**中文.** 有**普拉克索致可逆性心衰**的专门病例报告（Mov Disord Clin Pract 2014, PMID:30363854）及处方警示摘要（Prescrire Int 2013, PMID:24171217）；**普拉克索及多巴胺能药物致外周/腿部水肿**的多个报告与危险因素研究（Postgrad Med J 2016, PMID:26941267；Rev Med Chil 2017, PMID:29652961；Intern Med J 2021, PMID:33890374；Arch Neurol 2007, PMID:17420306）——水肿常被用来解释心衰关联（Consult Pharm 2015, PMID:25760663）；**罗匹尼罗致晕厥/低血压**（Cureus 2021, PMID:34268037）；**另一非麦角药吡贝地尔致低血压+心动过缓**（BMC Neurol 2018, PMID:30591018）；**罗匹尼罗致死毒理案例**（Forensic Toxicol 2022, PMID:36454487）。**说明书层面**：普拉克索旧版上市后清单含**心肌梗死、冠脉闭塞、室颤、室性早搏、心动过速、心脏骤停、左室衰竭**（openFDA，eff. 2016-10-13），现行版含**心力衰竭**与**晕厥**；罗匹尼罗说明书含**胸痛 4% vs 2%**与腿部水肿。**未检索到**归因于罗匹尼罗的急性 MI 专门病例报告或说明书条目。

### 6.2 FAERS disproportionality analysis (computed for this report) FAERS 不均衡性分析（本报告自行计算）

**EN.** Using the **openFDA FAERS** endpoint (`/drug/event.json`, data cut listed `last_updated = 2026-10-13`), I computed crude reporting odds ratios (ROR, 95% CI by Woolf's method) for the drugs vs the whole database. **Baseline:** total FAERS reports **20,692,305**; pramipexole-mentioning reports **21,182**; ropinirole-mentioning reports **21,217**.

| Reaction (MedDRA PT) | Pramipexole n | ROR (95% CI) | Ropinirole n | ROR (95% CI) |
|---|---:|---:|---:|---:|
| **Myocardial infarction** | 176 | **0.85 (0.73–0.98)** | 177 | **0.85 (0.73–0.99)** |
| Cardiac failure | 251 | **1.46 (1.29–1.66)** | 258 | **1.50 (1.33–1.70)** |
| Cardiac failure congestive | 106 | **1.29 (1.07–1.57)** | 165 | **2.02 (1.73–2.35)** |
| Cardiac arrest | 70 | 0.83 (0.66–1.05) | 95 | 1.13 (0.92–1.38) |
| Angina pectoris | 26 | 0.89 (0.61–1.31) | 27 | 0.92 (0.63–1.35) |
| Coronary artery disease | 37 | 1.19 (0.86–1.65) | 29 | 0.93 (0.65–1.34) |
| **Hypotension** | 636 | **2.94 (2.72–3.18)** | 452 | **2.07 (1.88–2.27)** |
| **Orthostatic hypotension** | 159 | **9.07 (7.75–10.61)** | 94 | **5.34 (4.36–6.54)** |
| **Syncope** | 171 | **1.68 (1.45–1.96)** | 198 | **1.95 (1.69–2.24)** |
| **Atrial fibrillation** | 144 | **1.46 (1.24–1.72)** | 209 | **2.13 (1.86–2.44)** |
| Ventricular fibrillation | 14 | 1.24 (0.74–2.10) | 10 | 0.89 (0.48–1.65) |
| Palpitations | 123 | 1.04 (0.87–1.24) | 110 | 0.93 (0.77–1.12) |
| Tachycardia | 169 | **1.29 (1.11–1.50)** | 156 | **1.19 (1.01–1.39)** |
| Oedema peripheral | 505 | **4.16 (3.81–4.54)** | 313 | **2.55 (2.28–2.85)** |
| Condition aggravated | 1113 | 3.79 (3.57–4.03) | 708 | 2.36 (2.19–2.54) |

**Reading of this table / 对表格的解读.** (1) **No MI signal** (ROR < 1 for both). (2) **No angina or coronary-artery-disease signal.** (3) **Robust signals for hypotension, orthostatic hypotension, syncope and peripheral oedema** — i.e. the *labelled* pharmacodynamic effects, which is the expected internal validation. (4) **HF/CHF signals** (ROR 1.3–2.0), consistent with the observational literature but **confounded by notoriety bias** (the 2012–2013 pramipexole-HF publicity), by PD comorbidity, and by oedema reporting; (5) an **atrial-fibrillation signal (1.5–2.1)** that is not in the labels and deserves follow-up but may reflect the elderly RLS/PD population.

**Mandatory caveats / 必须的告诫.** FAERS **cannot establish incidence or causality**. Disproportionality is never corrected for age, indication, comedication or reporting bias; a high ROR can arise purely because a drug's user population is elderly and comorbid, and the RORs here for *class-effect* terms (oedema, hypotension) are consistent with simple risk-factor confounding. This is a **signal-detection aid**, not a risk estimate.

**中文.** 用 **openFDA FAERS** 端点自行计算粗略报告比值比（ROR，Woolf 法 95% CI）。**基线**：FAERS 总报告 **20,692,305**；提及普拉克索 **21,182**、罗匹尼罗 **21,217**。**解读**：(1) **无心梗信号**（两者 ROR < 1）；(2) **无心绞痛/CAD 信号**；(3) **低血压、体位性低血压、晕厥、外周水肿信号稳健**——即说明书已载的药效学效应，构成预期的内部验证；(4) **心衰/充血性心衰信号**（ROR 1.3–2.0），与观察性文献一致，但受**关注偏倚**（2012–2013 年普拉克索-心衰的宣传）、PD 共病与水肿上报的混杂；(5) **房颤信号（1.5–2.1）**未见于说明书，值得随访，但可能反映老年 RLS/PD 人群。**必须的告诫**：FAERS **不能确立发生率或因果**；不均衡性从不校正年龄、适应证、合并用药或报告偏倚；**类效应**术语（水肿、低血压）的 ROR 与单纯危险因素混杂一致。这是**信号检测辅助**，不是风险估计。

---

## 7. Bottom line — graded recommendation table 结论——分级推荐表

**EN — direct answers to the seven questions.**
1. Both are non-ergot D2-family agonists (pramipexole D3-preferring; ropinirole more balanced), with negligible D1 activity; the coronary bed's dopamine receptor is predominantly **vasodilatory (D1)**; the ergot scaffold (→5-HT2B agonism) is the valvulopathy/vasoconstriction culprit, and non-ergot agents **lack it**.
2. **Orthostatic hypotension/syncope is the dominant, labelled CV liability**, and the pramipexole trials **excluded active cardiovascular disease**, so reassurance does not extend to CHD patients.
3. **Controlled/measured data do not show an MI or mortality excess**; HF signals exist at the observational level, are **pramipexole-predominant**, **not replicated for ropinirole** in the largest studies, and are **low-certainty**; RCT meta-analysis is negative; treated RLS had **lower** incident CVD than untreated RLS.
4. There is **real preclinical cardioprotection** (pramipexole mouse MIRI via AMPK/autophagy; rotigotine mouse myocardial ischaemia; D2/D3 post-conditioning biology) plus receptor-independent antioxidant/ferroptosis/mitochondrial actions — **but no human IRI/MI evidence, and dopaminergic organ-protection has a history of human failure (fenoldopam, renal-dose dopamine)**.
5. **Contraindications:** pramipexole = **none**; ropinirole = **hypersensitivity only**. Labelled CV warnings = **orthostatic hypotension, syncope, hypotension** (ropinirole), **falling asleep/somnolence, ICDs, dyskinesia, hallucinations** (both). **No MAOI interaction** is labelled for either.
6. **No MI signal** in controlled data or in the FAERS disproportionality computed here; label-level MI/arrhythmia entries exist only in pramipexole's old free-text postmarketing list; signals are for **hypotension/syncope/oedema/HF**, with a novel **AF** signal worth watching.
7. **Grade:** acceptable-with-monitoring in stable CHD; caution near-mi/in unstable CHD; mechanistic MIRI plausibility is **+1 (preclinical only)**.

**中文——七个问题直答。** 1）两者均为非麦角 D2 家族激动剂（普拉克索偏 D3，罗匹尼罗较均衡），D1 活性可忽略；冠脉多巴胺受体以**扩血管（D1）**为主；麦角碱骨架（→5-HT2B 激动）是瓣膜病/血管收缩的元凶，非麦角类**不含该骨架**。2）**体位性低血压/晕厥是主导的标示心血管风险**，且普拉克索试验**排除了活动性心血管疾病患者**，故不能外推至 CHD。3）**受控/实测数据未见 MI 或死亡率超标**；心衰信号仅存在于观察性层面，**以普拉克索为主**、在最大研究中**未被罗匹尼罗重复**、**证据确定性低**；RCT 荟萃分析为阴性；接受治疗的 RLS 未来 CVD 反而更低。4）存在**真实的临床前心脏保护**（普拉克索经 AMPK/自噬的小鼠 MIRI；罗替高汀小鼠心肌缺血；D2/D3 缺血后处理生物学）及不依赖受体的抗氧化/铁死亡/线粒体作用——**但无人 IRI/MI 证据，且多巴胺能器官保护有人体失败史（非诺多泮、肾剂量多巴胺）**。5）**禁忌**：普拉克索**无**；罗匹尼罗**仅过敏**。标示心血管警示：**体位性低血压、晕厥、低血压**（罗匹尼罗）、**突睡/嗜睡、ICD、异动症、幻觉**（两者）。**两药均未标示 MAOI 相互作用**。6）受控数据及本人 FAERS 不均衡性分析**均无心梗信号**；MI/心律失常条目仅见于普拉克索旧版上市后清单；信号集中在**低血压/晕厥/水肿/心衰**，另有值得关注的**房颤**新信号。7）**分级**：稳定 CHD 属"可接受但需监测"；不稳定 CHD/近 MI 属"慎用"；MIRI 机制可能性为 **+1（仅临床前）**。

Legend / 图例: ✅ acceptable 可接受 · ⚠️ caution / monitoring 慎用·需监测 · ⛔ avoid 避免 · 🎯 indicated 有适应证 · ❓ evidence insufficient 证据不足

| Clinical scenario / 临床情境 | Pramipexole 普拉克索 | Ropinirole 罗匹尼罗 | Basis / 依据 |
|---|:---:|:---:|---|
| **Stable CHD, no heart failure** 稳定 CHD，无心衰 | ⚠️ caution — lowest effective dose, titrate slowly, monitor BP/orthostatics | ⚠️ caution — same; syncope is the standout labelled risk | Labels §5.2/§5.3 (openFDA); PMID:22109939; PMID:23881697; PMID:33550813 (treated RLS had lower CVD) |
| **Stable CHD with LV dysfunction or HF** 稳定 CHD 合并左室功能不全/心衰 | ⚠️⚠️ higher caution — pramipexole is the agent carrying the observational HF signal | ⚠️ caution — HF signal weaker/absent for ropinirole in largest studies; oedema monitoring anyway | PMID:22109939; PMID:30363519; PMID:27796464; PMID:27222673; PMID:25760663 |
| **Recent MI (< 4–12 weeks) / unstable angina** 近期 MI（<4–12 周）/不稳定心绞痛 | ⛔ avoid unless a compelling indication and haemodynamically stable | ⛔ avoid | Trials **excluded active CVD** (label §5.2); hypotension/syncope risk in a fixed-stenosis heart; **no outcome data** |
| **Acute peri-PCI / MIRI window** 急性围 PCI / MIRI 窗口期 | ❓ no human evidence — do not use therapeutically outside a trial | ❓ no human evidence | No human IRI data; class precedent (fenoldopam, renal-dose dopamine) is negative (PMID:30601331; PMID:39607014; PMID:17670192) |
| **Coronary vasoconstriction / valvulopathy risk** 冠脉收缩/瓣膜病风险 | ✅ low — lacks the ergot scaffold; no valvulopathy signal (VigiBase/FAERS) | ✅ low — same | PMID:27796464; PMID:40065414; labels' "ergoline structure" attribution; PMID:15197703; PMID:19346455 |
| **Postmarketing MI / arrhythmia risk** 上市后 MI/心律失常风险 | ✅ no MI signal (label list is old free-text; FAERS MI ROR 0.85) | ✅ no MI signal (FAERS MI ROR 0.85) | §6.1–6.2 (this report); openFDA labels |
| **Heart rate / QTc at therapeutic dose** 治疗剂量心率/QTc | ✅ no labelled HR/QT effect | ✅ no labelled HR effect; ⚠️ *in-vitro* IKr block only at ≥1 µM (≈supra-therapeutic) | openFDA ropinirole ER label; Naunyn-Schmiedebergs Arch Pharmacol 2010, PMID:20668839 |
| **Impulse-control disorder / augmentation (RLS)** ICD / 加重（RLS） | ⚠️ counsel; strongest ICD signal in FAERS | ⚠️ counsel; augmentation 7–10%/yr | Label §5.3/§5.6/§5.9; PMID:40065414; PMID:41563785; PMID:29756335 |
| **Repositioning for MIRI (therapeutic hypothesis)** MIRI 重定位（治疗假设） | 🧪 **+1 preclinical only** — mouse MIRI infarct-size reduction via AMPK/autophagy; needs mechanistic + in-vivo confirmation | 🧪 **0/+1** — no direct MIRI data; class (rotigotine, D2 post-conditioning) supportive | PMID:27063800; PMID:30959104; PMID:26599746; PMID:23636617 |

**One-line verdict / 一句话结论:** in **stable CHD** these two non-ergot agonists are **acceptable with careful BP/orthostatic monitoring** ⚠️ (ropinirole's syncope signal and pramipexole's HF/oedema signal are the things to watch); they should be **avoided in the acute post-MI / unstable period** ⛔ and have **no human evidence** for a therapeutic MIRI effect ❓ — while the **preclinical MIRI signal for pramipexole is the genuine, and only, reason it surfaced in this screen** 🧪.
在**稳定 CHD** 中，这两个非麦角激动剂**可接受但须谨慎监测血压/体位变化** ⚠️（需关注罗匹尼罗的晕厥信号与普拉克索的心衰/水肿信号）；在**急性梗死后/不稳定期应避免** ⛔；对 MIRI 的治疗效应**无人体证据** ❓——而**普拉克索的临床前 MIRI 信号才是它在本次筛选中出现的唯一真实理由** 🧪。

---

## 8. Evidence gaps (stated honestly) 证据缺口（诚实说明）

**EN.**
1. **No RCT/cohort of pramipexole or ropinirole specifically in CHD or post-MI patients.** The pivotal pramipexole trials **excluded active cardiovascular disease** (label §5.2), so extrapolation is inference, not evidence.
2. **No human ischaemia-reperfusion (MIRI) study** of either drug; the cardioprotection evidence is **preclinical only** (PMID:27063800; PMID:30959104).
3. **No MI signal found anywhere** — but this is largely **absence of evidence**: MI was rarely a prespecified endpoint, cohorts were powered for PD, and follow-up was short.
4. **No measured 5-HT2B binding data** for pramipexole/ropinirole was retrieved; the non-ergot class distinction rests on scaffold absence + label attribution + differential epidemiology.
5. **HF evidence is contradictory** (negative RCT meta-analysis PMID:27222673 vs positive observational studies PMID:22109939/PMID:30363519 vs no-signal VigiBase PMID:27796464); the pramipexole HF signal is **not separable** from *oedema* and from *PD comorbidity*.
6. **FAERS analysis here is crude** (no adjustment, no case adjudication, notoriety-bias-prone) — a signal aid only.
7. **Ropinirole-specific mechanistic cardiac data are almost absent** — ropinirole has cardiac ion-channel data (PMID:20668839) but no MIRI model data found.

**中文。** 1）**无**针对 CHD 或梗死后患者的普拉克索/罗匹尼罗 RCT 或队列——关键试验**排除了活动性心血管疾病**（说明书 §5.2），外推属推断而非证据。2）**无**任一药物的人体缺血再灌注（MIRI）研究；心脏保护证据**仅临床前**（PMID:27063800；PMID:30959104）。3）**任何地方均未发现 MI 信号**——但这主要是**证据缺失**：MI 很少作为预设终点，队列按 PD 设计效能，随访偏短。4）**未检索到**普拉克索/罗匹尼罗的**实测 5-HT2B 结合数据**；非麦角类区分依赖骨架缺失 + 说明书归因 + 流行病学差异。5）**心衰证据自相矛盾**（RCT 荟萃阴性 PMID:27222673 vs 观察性阳性 PMID:22109939/PMID:30363519 vs VigiBase 无信号 PMID:27796464）；普拉克索心衰信号**无法与*水肿*及*PD 共病*分离**。6）本文 FAERS 分析为**粗略**（未校正、未判定病例、易受关注偏倚影响）——仅供信号参考。7）**罗匹尼罗特异的心脏机制数据几乎缺如**——仅有心脏离子通道数据（PMID:20668839），未发现 MIRI 模型数据。

---

## References 参考文献

**Labels / regulatory & pharmacovigilance databases 说明书与药物警戒数据库**
1. Pramipexole dihydrochloride tablets — US FDA label (openFDA drug-label API, `openfda.generic_name:pramipexole`; 37 labels inspected), effective 2026-09-09; older free-text postmarketing list, effective 2016-10-13. https://api.fda.gov/drug/label.json
2. Ropinirole tablets — US FDA label (openFDA drug-label API, `openfda.generic_name:ropinirole`; 45 labels inspected), effective 2026-09-04. https://api.fda.gov/drug/label.json
3. Mirapex SPL, DailyMed setid `3736d691-e8b9-4fdd-aeca-9e96cf907c51`; Requip XL SPL, DailyMed setid `2deb758c-659e-4e5a-a087-0e9c422f2fe1`. https://dailymed.nlm.nih.gov
4. **FAERS** disproportionality analysis computed for this report — openFDA `/drug/event.json`, data cut `last_updated 2026-10-13`; total DB 20,692,305 reports. https://api.fda.gov/drug/event.json
5. **ChEMBL** ligand binding data: pramipexole CHEMBL301265; ropinirole CHEMBL589. https://www.ebi.ac.uk/chembl/

**Receptor pharmacology & ergot/5-HT2B mechanism 受体药理与麦角/5-HT2B 机制**
6. The Serotonin 2B (5-HT2B) Receptor: Narrative Review of Preclinical and Clinical Evidence. Clin Pharmacol Ther 2026. PMID:42209433.
7. Severe multivalvular heart disease: a new complication of the ergot derivative dopamine agonists. Mov Disord 2004. PMID:15197703.
8. Cyproheptadine prevents pergolide-induced valvulopathy in rats. Am J Physiol Heart Circ Physiol 2009. PMID:19346455.
9. Lisuride, a dopamine receptor agonist with 5-HT2B receptor antagonist properties: absence of cardiac valvulopathy. Clin Neuropharmacol 2006. PMID:16614540.
10. The bulky N6 substituent of cabergoline is responsible for 5-HT2A/5-HT2B agonism. J Pharmacol Exp Ther 2011. PMID:21518772.
11. Assessment of valvulopathy in PD patients on pergolide and/or cabergoline. Clin Neurol Neurosurg 2007. PMID:17307289. (see also Neurology 2006, PMID:17030757; Can J Neurol Sci 2008, PMID:18574930; Can J Neurol Sci 2006, PMID:16583718)
11b. Risks of cardiac valve regurgitation and heart failure associated with ergot- and non-ergot-derived dopamine agonist use in PD: systematic review of observational studies. CNS Drugs 2015. PMID:26585874.
12. Ergot and non-ergot dopamine agonists and heart failure in PD (VigiBase). Eur J Clin Pharmacol 2017. PMID:27796464.
13. Comparative safety signals of dopamine agonists (FAERS). BMC Pharmacol Toxicol 2025. PMID:40065414.
14. Dopamine and vascular dynamics control. Curr Neurovasc Res 2011. PMID:21722093.
15. D1 receptor-mediated coronary vasodilation in the anaesthetized dog. J Pharmacol Exp Ther 1989. PMID:2563287.
16. Vasodilator action of dopamine and dopamine agonists in renal and coronary beds of the dog. Br J Pharmacol 1982. PMID:6982091.
17. Effect of dopamine and blockade of dopaminergic and adrenergic receptors on coronary blood flow in ischaemic rabbit myocardium. J Cardiovasc Pharmacol 1985. PMID:2410711.
18. Systemic, pulmonary and coronary haemodynamic actions of a novel dopamine receptor agonist in awake pigs. Br J Pharmacol 1997. PMID:9134223.
19. Pharmacokinetic evaluation of pramipexole. Expert Opin Drug Metab Toxicol 2011. PMID:21892895.

**Heart failure / CV safety epidemiology 心衰与心血管安全性流行病学**
20. Dopamine agonist use and the risk of heart failure (UK GPRD). Pharmacoepidemiol Drug Saf 2012. PMID:22109939.
21. Risk of heart failure associated with dopamine agonists: nested case-control (Taiwan NHIRD). Drugs Aging 2013. PMID:23881697.
22. Associations between cardiovascular events and nonergot dopamine agonists in PD (Cerner). Mov Disord Clin Pract 2016. PMID:30363519.
23. Non-ergot dopamine agonists do not increase the risk of heart failure in PD: meta-analysis of RCTs. J Clin Med Res 2016. PMID:27222673.
24. Non-ergot dopamine agonists and the risk of heart failure and other adverse cardiovascular reactions in PD. Brain Sci 2024. PMID:39199470.
25. Risk of heart failure following treatment with dopamine agonists in PD. Expert Opin Drug Saf 2014. PMID:24547918.
26. The link between non-ergot-derived dopamine agonists and heart failure: how strong is it? Consult Pharm 2015. PMID:25760663.
27. Pramipexole: heart failure. Prescrire Int 2013. PMID:24171217.
28. Evaluation of the effect of non-ergot dopamine agonists on left ventricular systolic function (speckle tracking). Anatol J Cardiol 2018. PMID:30297579.
29. Cardiovascular adverse events reported in the placebo arm of RCTs in PD. J Parkinsons Dis 2020. PMID:32116264.
30. Treating restless legs syndrome was associated with low risk of cardiovascular disease (MarketScan cohort). J Am Heart Assoc 2021. PMID:33550813.
31. RLS and risk of incident CVD in women and men. BMJ Open 2012. PMID:22447047.
32. Risk of CVD associated with an RLS diagnosis (Kaiser Permanente). Sleep 2015. PMID:26083613.
33. Incident cardiovascular events and death in individuals with RLS or PLMS: systematic review. Sleep 2017. PMID:28329227.
34. Associations of incident cardiovascular events with RLS and PLMS in older men. Sleep 2017. PMID:28199705.
35. Restless legs syndrome: a review. JAMA 2026. PMID:41563785.
36. Treatment of restless legs syndrome: evidence-based review (revised 2017). Mov Disord 2018. PMID:29756335.
37. Efficacy and safety of pramipexole in Chinese patients with RLS. Sleep Med 2012. PMID:22137119. (see also Sleep Med 2015, PMID:25533540)
38. Pramipexole versus ropinirole: polysomnographic acute effects in RLS. Mov Disord 2011. PMID:21370262.

**Hypotension / autonomic / cardiac electrophysiology 低血压、自主神经与心脏电生理**
39. Chronic dopaminergic treatment in RLS: does it affect the autonomic nervous system? Sleep Med 2015. PMID:26298781.
40. Ropinirole-associated orthostatic hypotension as cause of a prescribing cascade. Cureus 2021. PMID:34268037.
41. Hypotension and bradycardia, a serious adverse effect of piribedil. BMC Neurol 2018. PMID:30591018.
42. Factors to consider in the selection of dopamine agonists for older persons with PD. Drugs Aging 2019. PMID:30623310.
43. Differential effects of treatment combinations on cardiovascular dysfunction in PD. Acta Neurol Scand 2004. PMID:14653849.
44. Effects of ropinirole on action potential characteristics and underlying ion currents in canine ventricular myocytes. Naunyn Schmiedebergs Arch Pharmacol 2010. PMID:20668839.
45. Falls in Parkinson's disease subtypes. Int J Environ Res Public Health 2019. PMID:31234571.
46. Pramipexole in patients with early PD (PROUD): randomised delayed-start trial. Lancet Neurol 2013. PMID:23726851.
46b. Clinical evaluation of pramipexole in advanced Parkinson's disease: double-blind, placebo-controlled, parallel-group study. Neurology 1997. PMID:9222185.

**MIRI / IRI / mechanism 心肌缺血再灌注与机制**
47. Pramipexole pretreatment attenuates myocardial ischaemia/reperfusion injury through upregulation of autophagy. Biochem Biophys Res Commun 2016. PMID:27063800.
48. Dopamine D2 receptors contribute to cardioprotection of ischaemic post-conditioning via autophagy in isolated rat hearts. Int J Cardiol 2016. PMID:26599746.
49. Dopamine D2 receptor activation in ischaemic post-conditioning cardioprotection via PKC-ε. Mol Cell Biochem 2013. PMID:23636617.
50. Mediation of dopamine D2 receptor activation in post-conditioning-attenuated cardiomyocyte apoptosis. Exp Cell Res 2014. PMID:24412422.
51. Preconditioning cardioprotection in hyperlipidaemic rat hearts: role of dopamine D2 signalling. BMC Cardiovasc Disord 2015. PMID:26216101.
52. Effects of rotigotine and rotigotine extended-release microsphere therapy on myocardial ischaemic injury in mice. Eur J Pharm Sci 2019. PMID:30959104.
53. DRD3 agonist (pramipexole) reduces morphine-induced cardiac fibrosis. Biochem Biophys Res Commun 2020. PMID:32819568.
54. Loss of function in dopamine D3 receptor attenuates LV cardiac fibroblast migration and proliferation in vitro. Front Cardiovasc Med 2021. PMID:34708087.
55. Emodin reduces MIRI by modulating the RUNX1/miR-142-3p/DRD2 pathway. Exp Ther Med 2022. PMID:36561980.
56. D2 dopamine receptor antagonist raclopride induces non-canonical autophagy in cardiac myocytes. J Cell Biochem 2013. PMID:22886761.
57. DRD3 (but not DRD2) activates autophagy through MTORC1 inhibition. Autophagy 2020. PMID:31538542.
58. DRD3 induces transient mTORC1-dependent autophagy… Cells 2025. PMID:40358175.
59. Pramipexole inhibits astrocytic NLRP3 inflammasome activation via Drd3-dependent autophagy. Acta Pharmacol Sin 2023. PMID:35896696.

**Antioxidant / mitochondrial (MPTP model) / neuroprotection 抗氧化、线粒体与神经保护**
60. Pramipexole inhibits lipid peroxidation and reduces MPP⁺-induced substantia nigra injury. Neurosci Lett 2000. PMID:10704769.
61. Neuroprotection by pramipexole against dopamine- and levodopa-induced cytotoxicity. Life Sci 1999. PMID:10227583.
62. Pramipexole protects against MPTP toxicity in non-human primates. J Neurochem 2006. PMID:16464239.
63. Neuroprotective effects of the D2/D3 agonist pramipexole against post-ischaemic nigrostriatal degeneration. Brain Res 1996. PMID:9117424.
64. Pramipexole prevents ischaemic cell death via mitochondrial pathways in ischaemic stroke. Dis Model Mech 2019. PMID:31235613.
65. Pramipexole in global cerebral I/R: inhibiting ferroptosis. J Stroke Cerebrovasc Dis 2025. PMID:39490461.
66. Pramipexole alleviates ferroptosis via Nrf2/SLC7A11/GPX4 in OGD/R. Sci Rep 2026. PMID:41741591.
67. Effects of different doses of pramipexole on mitochondrial membrane potential in cerebral I/R. J Stroke Cerebrovasc Dis 2023. PMID:37105127.
68. Pramipexole + N-acetylcysteine in global cerebral I/R in rats. Iran J Basic Med Sci 2018. PMID:29942446.
69. D2/D3 receptor activation in ischaemic preconditioning-mediated brain protection. Brain Res 2020. PMID:32919985.
70. Protective effects of piribedil (D2/D3 agonist) in global cerebral ischaemia in rats. Neurosci Lett 2018. PMID:30107200.
71. Pramipexole mitigates rotenone-evoked hepatic injury (antioxidant/NLRP3). Int Immunopharmacol 2026. PMID:42391877.
72. Pramipexole halts leflunomide nephrotoxicity. Eur J Pharmacol 2025. PMID:41109403.
73. Pramipexole reduces inflammation in experimental animal models. Immunopharmacol Immunotoxicol 2017. PMID:28162057.
74. Safety and efficacy of dexpramipexole in eosinophilic asthma (EXHALE RCT). J Allergy Clin Immunol 2023. PMID:37277072.
75. Pramipexole and dexpramipexole reverse chronic allodynia. Neurosci Lett 2023. PMID:37558176.
76. Dopamine, immunity, and disease. Pharmacol Rev 2023. PMID:36757901.

**Negative translational precedent for dopaminergic organ protection 多巴胺能器官保护的人体阴性先例**
77. Does fenoldopam protect kidney in cardiac surgery? Meta-analysis with trial sequential analysis. Shock 2019. PMID:30601331.
78. Fenoldopam for preventing and treating acute kidney injury. Cochrane Database Syst Rev 2024. PMID:39607014.
79. Fenoldopam for renal protection in cardiac surgery: pharmacology and evolving perspectives. J Clin Med 2024. PMID:39407923.
80. Low-dose (renal-dose) dopamine in the critically ill patient. Interact Cardiovasc Thorac Surg 2004. PMID:17670192. (see also Curr Surg 2006, PMID:16757377)

**Case reports / oedema 病例报告与水肿**
81. Pramipexole-induced reversible heart failure. Mov Disord Clin Pract 2014. PMID:30363854.
82. A case of severe leg oedema in a patient with PD treated with pramipexole. Postgrad Med J 2016. PMID:26941267.
83. Severe leg oedema associated with dopaminergic drugs in PD. Rev Med Chil 2017. PMID:29652961.
84. Pramipexole-induced oedema: the importance of renal dopamine. Intern Med J 2021. PMID:33890374.
85. Risk factors for the development of pedal oedema in patients using pramipexole. Arch Neurol 2007. PMID:17420306.
86. Ropinirole involved in a fatal case: blood and urinary concentrations. Forensic Toxicol 2022. PMID:36454487.

**Other / context 其他**
87. Risk of cardiovascular disease in males with prolactinoma (ergot-DA context). Clin Endocrinol (Oxf) 2018. PMID:29044586.
88. Cardiovascular risk with prokinetics (D2 antagonists) — self-controlled case series. Cardiovasc Toxicol 2025. PMID:40684411.
89. Drug safety analysis in a real-life PD cohort with polypharmacy. CNS Drugs 2017. PMID:29139041.
90. Impulse control and related behavioural disorders with different dopamine agonists. Clin Park Relat Disord 2024. PMID:38292814.
