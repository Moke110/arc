# MIRI Drug Candidate Literature Screening Report

**Generated:** 2026-06-12
**Source:** `fda_miri_scores/scores/fda_filtered_scores.csv` (top 50 drugs)
**Methods:** Web-based literature search via parallel subagents, 20 agents total (10 Phase 1 + 10 Phase 2)

---

## Executive Summary

- **50 drugs screened**, 25 passed Phase 1 risk screening (50% pass rate)
- **25 drugs excluded**: QT prolongation/TdP (8), CV class risks (7), market withdrawal (4), mortality signal (2), VTE risk (2), AF risk (2)
- **Phase 2 evidence**: 6 drugs with ≥3★ anti-MIRI evidence, 11 drugs with ≥3★ anti-IRI evidence
- **Standout candidates**: Doxycycline (5★ anti-MIRI, 4★ anti-IRI), Ramatroban (4★ anti-MIRI, 4★ anti-IRI), Terazosin (4★ anti-MIRI, 4★ anti-IRI)
- **Clinical-stage evidence**: Doxycycline has the TIPTOP trial (Phase II, n=110), Tirofiban has multiple RCTs in STEMI

---

## Methods

### Research Methodology

All research was conducted via web search (WebSearch + WebFetch tools), targeting PubMed, FDA labels, clinical trial registries, and published literature. Each subagent performed independent multi-query searches per drug.

### Grading Criteria

| Grade | Stars | Phase 2 Criteria |
|-------|-------|------------------|
| 临床对照/大规模队列 | ★★★★★ | Clinical controlled trials or large-scale cohort studies with MIRI/IRI benefit |
| 小规模临床/动物实验 | ★★★★ | Small clinical non-controlled studies, or animal IRI model studies |
| 病例报告/细胞实验 | ★★★ | Case reports, or in vitro cell experiments (hypoxia/reoxygenation) |
| 生信数据分析线索 | ★★ | Bioinformatics/data analysis hints only, no direct experiments |
| 无报道 | ★ | No published evidence found |

### Exclusion Criteria (Phase 1)

A drug is EXCLUDED from Phase 2 if it has any of:
- Confirmed cardiotoxicity (QT prolongation, TdP, arrhythmogenic, direct myocardial damage, heart failure risk)
- Market withdrawal in US, EU, or China for safety reasons
- FDA black-box warning specifically for cardiovascular mortality
- Practical unavailability in both US and China markets (for discontinued/obsolete drugs)

---

## Phase 1: Risk Screening Summary

### Complete Results Table (50 Drugs)

| # | Drug Name | Chinese Name | CN Status | FDA Status | Risk Flags | Verdict |
|---|-----------|-------------|-----------|------------|------------|---------|
| 1 | PEMETREXED DISODIUM HEMIPENTAHYDRATE | 培美曲塞二钠 | Approved (多generic) | Approved (Alimta brand discontinued, commercial) | Arrhythmia <1%, HF 1-10%, no QT | ✅ PASS |
| 2 | ERYTHROMYCIN GLUCEPTATE | 红霉素葡庚糖酸盐 | Not marketed (obsolete salt) | Discontinued (Ilotycin) | QT/TdP (CredibleMeds highest tier) | ❌ EXCLUDED |
| 3 | PENICILLIN G BENZATHINE | 苄星青霉素 | Approved (多generic) | Approved (shortage, not discontinued) | BB for IV route only (mechanical, not pharmacologic) | ✅ PASS |
| 4 | APREPITANT | 阿瑞匹坦 | Approved (多generic) | Approved (brand oral discontinued, commercial) | Clean cardiac, negative QTc study | ✅ PASS |
| 5 | DOXYCYCLINE | 多西环素/强力霉素 | Approved (多generic) | Approved (Vibramycin discontinued, commercial) | Clean cardiac profile, no QT | ✅ PASS |
| 6 | RIBOFLAVIN 5'-PHOSPHATE | 核黄素磷酸钠 | Approved | Approved (ophthalmic + IV multivitamin) | FMN form safe; FAD form has transient AV block | ✅ PASS |
| 7 | ZOLPIDEM TARTRATE | 酒石酸唑吡坦 | Approved (Schedule II controlled) | Approved (Schedule C-IV) | QTc prolongation, TdP case reports, hERG inhibition | ❌ EXCLUDED |
| 8 | ATENOLOL | 阿替洛尔 | Approved (多generic) | Approved | Negative inotrope, bradycardia, AV block, BB for abrupt cessation MI | ❌ EXCLUDED |
| 9 | IBREXAFUNGERP CITRATE | 枸橼酸艾瑞芬净 | Approved 2025-04 | Approved 2021 (supply disruption) | Clean cardiac, negative QTc. BB for embryo-fetal (non-cardiac) | ✅ PASS |
| 10 | CEFIDEROCOL SULFATE TOSYLATE | 头孢德罗/头孢地尔 | Approved 2026-01 | Approved 2019 | AF ≥4%, all-cause mortality signal (+13.3% at 49d) | ❌ EXCLUDED |
| 11 | DIRITHROMYCIN | 地红霉素 | Approved (多generic) | Discontinued (~2002, commercial) | Lowest macrolide QT risk; minimal CYP3A4 inhibition | ✅ PASS |
| 12 | LORACARBEF | 氯碳头孢 | Not approved | Discontinued (~2006, commercial) | Beta-lactam; clean cardiac profile | ✅ PASS |
| 13 | OMEGA-3-ACID ETHYL ESTERS | ω-3脂肪酸乙酯 | Approved 2022 | Approved 2004 | AF risk (dose-dependent, EMA/MHRA/FDA recognized) | ❌ EXCLUDED |
| 14 | RELEBACTAM | 瑞来巴坦 (亚胺西瑞) | Approved 2024-12 | Approved 2019 | Negative QTc study at 4.6x dose | ✅ PASS |
| 15 | ESTETROL | 雌四醇 | Not approved | Approved 2021 | BB: CV events (class CHC); VTE/arterial thromboembolism | ❌ EXCLUDED |
| 16 | DOCETAXEL | 多西他赛 | Approved (多generic) | Approved | VT (fatal with combo chemo), MI, cardiac tamponade | ❌ EXCLUDED |
| 17 | CEPHALEXIN | 头孢氨苄 | Approved (多generic) | Approved | No cardiac risk; standard 1st-gen cephalosporin | ✅ PASS |
| 18 | PIBRENTASVIR | 哌仑他韦 (艾诺全) | Approved 2019 | Approved 2017 | Negative QTc study (TQT verified). BB for HBV (non-cardiac) | ✅ PASS |
| 19 | DICLOXACILLIN | 双氯西林 | API only (limited use) | Approved | No cardiac risk; cholestatic hepatitis (rare, non-cardiac) | ✅ PASS |
| 20 | PIRBUTEROL ACETATE | 吡布特罗 | NOT approved | Discontinued | Beta-agonist: QT, arrhythmia, angina, MI. Not available US or CN | ❌ EXCLUDED |
| 21 | TRIMETHOPRIM SULFATE | 甲氧苄啶 | Approved (base form only) | Ophthalmic combo only; oral discontinued | QT signal only in co-trimoxazole combo (low incidence) | ✅ PASS |
| 22 | DISOPYRAMIDE | 丙吡胺 | Approved (3 manufacturers) | Approved | BB: CV mortality (CAST trial); QT/TdP; negative inotropy; CHF | ❌ EXCLUDED |
| 23 | PIPAMAZINE | 哌帕马嗪 | Never approved | Withdrawn 1969 | Hepatotoxicity withdrawal; profound hypotension; never in CN | ❌ EXCLUDED |
| 24 | PEFICITINIB HYDROBROMIDE | 氢溴酸吡西替尼 | Approved 2024-07 (import) | Not FDA approved | JAK inhibitor class MACE/cancer/VTE risk (PMDA Important Potential Risk) | ❌ EXCLUDED |
| 25 | TERAZOSIN | 特拉唑嗪 | Approved (多generic) | Approved (generics) | Excellent cardiac profile; first-dose hypotension only | ✅ PASS |
| 26 | THIETHYLPERAZINE MALEATE | 硫乙拉嗪 | Not approved | Discontinued | Phenothiazine: QT/TdP, sudden death. Discontinued US. Not in CN | ❌ EXCLUDED |
| 27 | HYDROXYSTILBAMIDINE ISETHIONATE | 羟芪巴脒 | Not approved | Discontinued | Obsolete antiprotozoal; hepatotoxicity/nephrotoxicity | ❌ EXCLUDED |
| 28 | LURBINECTEDIN | 芦比替定/赞必佳 | Approved 2024 | Approved 2020 | Clean cardiac; negative QTc study at therapeutic dose | ✅ PASS |
| 29 | ISOETHARINE MESYLATE | 异他林 | Withdrawn | Discontinued | Old beta-agonist: QT prolongation, discontinued US, withdrawn CN | ❌ EXCLUDED |
| 30 | RITODRINE | 利托君/安宝 | Approved (still marketed CN) | Withdrawn US | Pulmonary edema (~4%), MI, maternal deaths, unfavorable risk-benefit | ❌ EXCLUDED |
| 31 | TIROFIBAN HYDROCHLORIDE | 盐酸替罗非班 | Approved (18+ generics) | Approved 1998 | No QT/TdP. Bleeding risk only (on-target antiplatelet effect) | ✅ PASS |
| 32 | HYDROXYZINE PAMOATE | 盐酸羟嗪 | Approved (limited generics) | Approved | hERG blocker, QT/TdP (59 case reports, 8 fatal). MHRA/Health Canada advisories | ❌ EXCLUDED |
| 33 | CEFTRIAXONE SODIUM | 头孢曲松钠 | Approved (多generic) | Approved (brand Rocephin discontinued, business) | No standalone QT risk. DD-I with lansoprazole only (avoidable) | ✅ PASS |
| 34 | LEVONORGESTREL | 左炔诺孕酮 | Approved (OTC, 10+ generics) | Approved (OTC) | Excellent cardiac safety. CV risk in COCs driven by estrogen component | ✅ PASS |
| 35 | ESCITALOPRAM | 艾司西酞普兰 | Approved (20+ generics) | Approved 2002 | Dose-dependent QTc (+4.5ms at 10mg; +10.7ms at 30mg). TdP postmarketing | ❌ EXCLUDED |
| 36 | IRINOTECAN HYDROCHLORIDE | 盐酸伊立替康 | Approved (多generic) | Approved 1996 | MI ~2% (postmarketing, often confounded by 5-FU). No QT warning | ✅ PASS |
| 37 | HEXAMINOLEVULINATE HYDROCHLORIDE | 盐酸氨酮戊酸己酯 | Approved 2024/2026 | Approved 2010 | Minimal systemic absorption (~7%). No cardiac signal | ✅ PASS |
| 38 | LEUCOVORIN CALCIUM | 亚叶酸钙 | Approved (多generic) | Approved | No cardiac signal; decades of widespread use | ✅ PASS |
| 39 | ESTRAMUSTINE PHOSPHATE SODIUM | 磷酸雌莫司汀 | Approved (import only) | DISCONTINUED 2024-09 | MI 3.2%, TE ~7%, BB for fatal/nonfatal MI. US discontinued | ❌ EXCLUDED |
| 40 | THIETHYLPERAZINE MALATE | 硫乙拉嗪/马来酸硫乙拉嗪 | Not approved | Discontinued | Phenothiazine: QT/TdP. Discontinued US. Not approved CN | ❌ EXCLUDED |
| 41 | RAMATROBAN | 雷马曲班 | CTA approved, not marketed | Not FDA approved (Japan only) | No cardiotoxicity; cardioprotective evidence. Bleeding risk (antiplatelet) | ✅ PASS |
| 42 | LAZERTINIB | 兰泽替尼/拉泽替尼 | Approved 2025 | Approved 2024 | VTE 36% (fatal 0.5%), mandatory anticoagulation | ❌ EXCLUDED |
| 43 | NALIDIXIC ACID | 萘啶酸 | Approved (obsolete) | Discontinued | QTc prolongation. EU marketing authorization SUSPENDED 2018 | ❌ EXCLUDED |
| 44 | DIPIVEFRIN | 地匹福林 | Approved | Discontinued US | Epinephrine prodrug: tachycardia, hypertension, arrhythmias. CV contraindication | ❌ EXCLUDED |
| 45 | NATAMYCIN | 那他霉素 | Approved | Approved | Negligible systemic absorption. No cardiac concerns | ✅ PASS |
| 46 | ADENOSINE PHOSPHATE | 环磷腺苷 | Approved (cAMP form) | FDA withdrawal of non-cyclic forms | Parent compound: cardiac arrest/arrhythmia/MI risk; cAMP also has arrhythmia warnings | ❌ EXCLUDED |
| 47 | IBANDRONATE SODIUM | 伊班膦酸钠 | Approved (多generic) | Approved 2003 | No intrinsic proarrhythmic risk. Hypocalcemia-mediated QT (class effect, manageable) | ✅ PASS |
| 48 | DEXTROAMPHETAMINE SACCHARATE | 右旋安非他命 | Controlled substance (Category I) | Approved (Schedule II) | BB: CV events (Class II amphetamine). CN Category I (effectively banned) | ❌ EXCLUDED |
| 49 | MIGLUSTAT | 麦格司他/泽维可 | Approved 2016 | Approved 2003 | Minimal hERG (IC50 >1000 uM). Possible QTcF signal inconclusive | ✅ PASS |
| 50 | METOCURINE IODIDE | 甲筒箭毒 | Not approved | Discontinued | Discontinued decades ago; histamine-release → hemodynamic instability | ❌ EXCLUDED |

### Exclusion Summary

**25 of 50 drugs EXCLUDED (50%).** Breakdown by exclusion reason:

| Reason Category | Count | Drugs |
|-----------------|-------|-------|
| QT prolongation / TdP | 8 | Erythromycin, Zolpidem, Thiethylperazine Maleate, Hydroxyzine, Escitalopram, Thiethylperazine Malate, Nalidixic Acid, Isoetharine |
| CV class risk (BB/negative inotrope/VTE) | 7 | Atenolol, Disopyramide, Peficitinib, Lazertinib, Dextroamphetamine, Dipivefrin, Docetaxel |
| Market withdrawal / Discontinued | 4 | Pipamazine, Pirbuterol, Hydroxystilbamidine, Metocurine |
| Mortality signal | 2 | Cefiderocol, Ritodrine |
| AF risk | 2 | Omega-3, Estetrol |
| Other (MI/TE/arrhythmia) | 2 | Estramustine, Adenosine Phosphate |

---

## Phase 2: Evidence Summary Table (25 PASS Drugs)

| # | Drug Name | Anti-MIRI Grade | Anti-MIRI Key Finding | Anti-IRI Grade | Anti-IRI Organ(s) |
|---|-----------|-----------------|----------------------|----------------|-------------------|
| 1 | Pemetrexed | ★ | No evidence; cardiotoxic | ★ | None |
| 2 | Penicillin G Benzathine | ★ | No evidence | ★★★★ | Intestinal (H,K-ATPase) |
| 3 | Aprepitant | ★★★★ | Rat I/R: PI3K-Akt-GSK-3β-HIF-1α | ★★★★ | Hepatic, Pulmonary, Cerebral |
| 4 | Doxycycline | ★★★★★ | TIPTOP trial + animal I/R: MMP-2/9 | ★★★★ | Renal, Hepatic, Cerebral, Intestinal, Pulmonary |
| 5 | Riboflavin 5'-Phosphate | ★★★★ | Animal I/R: ATP preservation, mitochondrial | ★★★★ | Renal, Hepatic, Cerebral, Pulmonary |
| 6 | Ibrexafungerp Citrate | ★ | No evidence | ★ | None |
| 7 | Dirithromycin | ★★ | Macrolide class hints only | ★★★ | Cerebral (macrolide class) |
| 8 | Loracarbef | ★ | No evidence | ★ | None |
| 9 | Relebactam | ★ | No evidence | ★★ | Cerebral (sulbactam analogy) |
| 10 | Cephalexin | ★ | No evidence | ★★ | Weak hints |
| 11 | Pibrentasvir | ★ | No evidence | ★ | None |
| 12 | Dicloxacillin | ★ | No evidence | ★ | None |
| 13 | Trimethoprim Sulfate | ★ | No evidence | ★ | None |
| 14 | Terazosin | ★★★★ | PGK1 activation + stroke protection | ★★★★ | Cerebral, Renal |
| 15 | Lurbinectedin | ★ | No evidence | ★ | None |
| 16 | Tirofiban HCl | ★★★★ | Animal I/R + clinical RCTs in STEMI | ★★★★ | Renal, Hepatic, Cerebral |
| 17 | Ceftriaxone Sodium | ★★★ | Cardiac injury model, neuroprotective GLT-1 | ★★★★ | Cerebral |
| 18 | Levonorgestrel | ★ | No evidence | ★★ | Cerebral (epidemiological only) |
| 19 | Irinotecan HCl | ★ | Cardiotoxic, not protective | ★ | IRI inducer |
| 20 | Hexaminolevulinate HCl | ★★★ | Parent 5-ALA proven in mouse I/R (HO-1) | ★★★★ | Renal, Hepatic, Cerebral (via 5-ALA) |
| 21 | Leucovorin Calcium | ★★ | Inferential from folic acid I/R | ★★ | Hepatic (old study) |
| 22 | Ramatroban | ★★★★ | Rat MIRI model: 53% infarct reduction | ★★★★ | Splanchnic/Intestinal, Pulmonary |
| 23 | Natamycin | ★ | No evidence | ★ | None |
| 24 | Ibandronate Sodium | ★ | No MIRI evidence; MI risk signal | ★ | None |
| 25 | Miglustat | ★★★ | Cardiac fibrosis model (ISO, not I/R) | ★★ | Organ-protective in non-IRI models |

### Evidence Grade Distribution

| Stars | Anti-MIRI (Track A) | Anti-IRI (Track B) |
|-------|--------------------|--------------------|
| ★★★★★ | 1 (Doxycycline) | 0 |
| ★★★★ | 5 (Aprepitant, Riboflavin, Tirofiban, Terazosin, Ramatroban) | 9 |
| ★★★ | 3 (Ceftriaxone, Hexaminolevulinate, Miglustat) | 1 |
| ★★ | 2 (Dirithromycin, Leucovorin) | 6 |
| ★ | 14 | 9 |

---

## Detailed Evidence — Drugs with ≥3★ in Either Track

### 1. DOXYCYCLINE (多西环素) — ⭐⭐⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (5★ Clinical):**

Doxycycline is the most extensively studied drug in our list for MIRI, with evidence from preclinical to clinical trials. As the most potent MMP inhibitor in the tetracycline class, it inhibits MMP-2 and MMP-9 which degrade cardiac contractile proteins (troponin I, MLC1, MLC2, titin) during I/R injury.

- **TIPTOP trial** (Cerisano et al., *Eur Heart J*, 2014; n=110, open-label Phase II): Doxycycline 100mg BID for 7 days post-PCI in STEMI patients (LVEF <40%) significantly reduced adverse LV remodeling: LVEDVi change 0.4% vs 13.4% (P=0.012), reduced infarct size (5.5% vs 10.4%, P=0.052).
- **SALVAGE-MI trial** (2022, n=103, double-blind): Upstream (pre-PCI) doxycycline showed no benefit and a trend toward 28% larger infarct size, suggesting timing-dependence.
- **DOXY-STEMI trial** (NCT03508232): Ongoing, results expected late 2025.
- Multiple animal models (rat LAD, pig) consistently show infarct size reduction and MMP suppression.

**Anti-IRI (4★ Multi-Organ):**

Broadest multi-organ IRI protection evidence among all 25 PASS drugs:
- **Renal**: Low-dose doxycycline (3-10 mg/kg) prevents GFR reduction, normalizes BUN/creatinine, preserves (Na+K)-ATPase (Kucuk 2009, Cortes 2018)
- **Hepatic**: Suppresses IRI even given 1h after reperfusion. Unique mechanism: mitochondrial calcium uniporter (MCU) inhibition (Schwartz et al., 2013)
- **Cerebral**: Reduces infarct volume from ~45.5% to ~20.7% in stroke-prone rats via MMP-2/9 inhibition, BBB preservation
- **Intestinal**: Mesenteric IRI protection, including human Caco-2 cell hypoxia/reoxygenation
- **Pulmonary**: Attenuates lung IRI by decreasing MMP-2/9

**Verdict**: Strongest overall candidate. Clinical repurposing requires optimal timing (post-reperfusion, not pre-PCI based on SALVAGE-MI).

---

### 2. RAMATROBAN (雷马曲班) — ⭐⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (4★ Animal):**

- Squadrito et al. (1993, *Inflamm Res*): In a rat left main coronary artery occlusion (1h) + reperfusion (1h) model, ramatroban (30-60 mg/kg IV) dose-dependently improved survival from 45% (vehicle) to 80%, reduced myocardial necrosis and serum CPK by ~53%, reduced cardiac MPO (neutrophil infiltration) by >65%.
- Dual mechanism: TXA2 (TP) receptor antagonist + DP2/CRTH2 antagonist, reducing platelet aggregation, vasoconstriction, and neutrophil-mediated inflammation.
- Approved in Japan (Baynas) for allergic rhinitis since 2000 → established human safety data enables repurposing.

**Anti-IRI (4★ Multi-Organ):**

- Splanchnic/intestinal IRI protection demonstrated in multiple preclinical models
- Secondary pulmonary protection following splanchnic IRI
- Hindlimb IRI protection

**Verdict**: Strong preclinical MIRI data with established human safety from allergic rhinitis indication. Not yet approved in US or China. Key repurposing candidate — the mechanism (dual TP/DP2 antagonism) directly addresses inflammatory and thrombotic components of MIRI.

---

### 3. TERAZOSIN (特拉唑嗪) — ⭐⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (4★ Animal + Mechanistic):**

- Lee & Su (2008): Direct cardiac I/R study (Langendorff rat heart). Terazosin was neutral on infarct size (unlike prazosin/doxazosin which worsened it, and bunazosin which improved it). This neutrality may actually be advantageous compared to other alpha-blockers.
- Chen et al. (2015, *Nature Chemical Biology*): Landmark discovery — terazosin activates phosphoglycerate kinase 1 (PGK1), the first ATP-generating enzyme of glycolysis. PGK1 activation increases cellular ATP, fueling Hsp90 chaperone activity and enhancing multi-stress resistance. This effect is independent of alpha-1 receptor blockade (confirmed with terazosin analog with 1000× reduced alpha-1 affinity).
- PGK1-Hsp90 mechanism is directly relevant to cardiac ischemia where glycolysis becomes the critical energy source.

**Anti-IRI (4★ Multi-Organ):**

- **Cerebral IRI**: Two independent mechanisms demonstrated: (1) PGK1-Hsp90 → increased ATP → stress resistance (Chen et al., 2015); (2) Mitophagy enhancement → NLRP3 inflammasome inhibition → reduced neuronal pyroptosis. **Phase 3 clinical trial (NCT06863571) recruiting** for ischemic stroke.
- **Renal IRI**: Protected against contrast-induced AKI (CI-AKI) in mouse CKD model, reducing serum creatinine 11.3% (p=0.027) and urinary KIM-1 26-28% (Khan & Misra et al., 2021).

**Verdict**: Uniquely positioned — decades of clinical use for BPH/hypertension provides extensive human safety data, and the PGK1 mechanism discovered in 2015 opens a novel cardioprotective axis. Biphasic dose-response (inhibitory at high concentrations) requires careful dose optimization for MIRI.

---

### 4. APREPITANT (阿瑞匹坦) — ⭐⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (4★ Animal):**

- Qian & Liu (2022, *Acta Cir Bras*): In rat Langendorff I/R model, aprepitant (5-20 mg/kg) dose-dependently reduced LDH-1 and CK-MB, restored HIF-1α and p-Akt, and normalized p-GSK-3β/GSK-3β ratio. PI3K inhibitor LY294002 abolished protection, confirming PI3K-Akt-GSK-3β-HIF-1α pathway.
- Additional preclinical: cardioprotection in oxidized LDL-induced endothelial dysfunction, viral myocarditis, and doxorubicin-induced cardiomyopathy models.

**Anti-IRI (4★ Multi-Organ):**

- **Hepatic IRI** (2025): Two independent studies. Aprepitant dose-dependently reduced ALT/AST, hepatic MDA, TNF-α, IL-6, caspase-3; increased TAC; enhanced Sirt1/Nrf2. 100% survival at 20 mg/kg vs 62.5% control.
- **Pulmonary** (remote lung injury post hepatic IRI): Suppressed NLRP3 inflammasome/IL-1β, reduced lung MDA, increased GSH.
- **Cerebral IRI**: Preserved hippocampal CA1 neurons in rat BCCAO; improved learning/memory, PPAR-γ-dependent neuroprotection (Kaur et al., 2016).

**Verdict**: Extensive preclinical evidence across cardiac, hepatic, cerebral, and pulmonary IRI. NK1 antagonism is a well-validated target. Human safety established from antiemetic use. Strong repurposing candidate.

---

### 5. RIBOFLAVIN 5'-PHOSPHATE (核黄素磷酸钠) — ⭐⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (4★ Animal):**

- Kotegawa et al. (1994): In Langendorff perfused rat heart (30 min global ischemia + 30 min reperfusion), FMN preserved ATP and creatine phosphate, restored intracellular pH, and preserved mitochondrial respiratory function.
- Additional Chinese studies: Riboflavin preconditioning preserved coronary endothelial integrity (reduced vWF, increased NO/eNOS/SOD), reduced MDA/LDH, and suppressed platelet P-selectin during I/R.

**Anti-IRI (4★ Multi-Organ):**

- **Renal**: Reduced MDA, 8-OHdG, MPO, caspase-3; restored GSH (Ayaz Adakul et al., 2019)
- **Hepatic**: Reduced ALT/AST, oxidative/nitrosative damage; restored GSH; suppressed eNOS/iNOS (Sanches et al., 2014)
- **Cerebral**: Riboflavin kinase reduced infarction, edema, neuronal apoptosis via caspase-12/caspase-3/ER stress pathway (Zou et al., 2025)
- **Pulmonary**: Reduced vascular permeability 53%, alveolar hemorrhage 70%, lipid peroxidation 45% in hind-limb I/R lung injury (Seekamp et al., 1999)

**Verdict**: Strong preclinical data across all major organ systems. Mechanism (mitochondrial energy preservation + antioxidant) is fundamental to IRI pathophysiology. No human IRI trials exist. Vitamin safety profile is excellent.

---

### 6. TIROFIBAN HYDROCHLORIDE (盐酸替罗非班) — ⭐⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (4★ Clinical + Animal):**

- **Clinical RCTs**: Multiple meta-analyses demonstrate that early (upstream/pre-hospital) tirofiban during primary PCI for STEMI improves ST-segment resolution, TIMI myocardial perfusion grade, reduces microvascular obstruction (confirmed by cardiac MRI, 2025), lowers peak CK-MB/troponin, and reduces 30-day MACE by ~40%. A 2024 meta-analysis of 11 RCTs (1,802 patients): IC low-dose tirofiban reduces MACE (RR 0.595, p=0.001).
- **Animal I/R**: Activates RISK pathway — phosphorylates PKC-ε, PI3K/Akt, p38 MAPK, ERK1/2 — protecting cardiomyocytes from lethal reperfusion injury. Reduced myocardial necrosis, apoptosis, oxidative stress in rabbit I/R model.

**Anti-IRI (4★ Multi-Organ):**

- **Renal IRI**: Reduced serum creatinine, BUN, ROS, apoptosis, leukocyte infiltration via eNOS/iNOS restoration (Guan et al., 2015)
- **Hepatic IRI**: Preserved mitochondrial function, reduced AST/ALT/LDH, restored CAT/SOD/GSH; attenuated secondary pulmonary vascular permeability (Canedo et al., 2007)
- **Cerebral IRI**: Reduced infarct volume, neuronal apoptosis, neuroinflammation (IL-1, IL-6, TNF-α) via microglial regulation (2024)

**Verdict**: Already used in ACS/STEMI — closest to clinical translation. GP IIb/IIIa inhibition provides both antiplatelet and direct cardiomyocyte protection. Bleeding risk is manageable and already part of standard clinical practice.

---

### 7. CEFTRIAXONE SODIUM (头孢曲松钠) — ⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (3★ Indirect):**

- AlAsmari et al. (2024): Mitigated fentanyl-induced cardiac injury in mice (reduced troponin I, hypertrophic genes, caspase-3, IL-6, NF-κB) via xCT/GLT-1 upregulation and Nrf2/ARE antioxidant signaling. Not a direct I/R model.
- Clinical trial investigating ceftriaxone + amoxicillin in ACS/STEMI/NSTEMI, but rationale targets bacterial involvement in plaques, not direct cardiomyocyte protection.

**Anti-IRI (4★ Cerebral):**

- Thone-Reineke et al. (2008): Single ceftriaxone injection 90 min post-MCAO in rats reduced 24h mortality from 34.5% to 0%, decreased infarct size.
- Mechanism: GLT-1 (EAAT2) glutamate transporter upregulation → reduced excitotoxic glutamate during ischemia.
- Non-antibiotic beta-lactam analogs (MC-100093) in development to avoid chronic antibiotic use.

**Verdict**: Strong cerebral IRI data, but cardiac I/R evidence remains indirect. The GLT-1 mechanism is compelling for neuroprotection. Caution: lansoprazole interaction (QT). Non-antibiotic analogs may be preferable for chronic use.

---

### 8. HEXAMINOLEVULINATE HYDROCHLORIDE (盐酸氨酮戊酸己酯) — ⭐⭐⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (3★ via 5-ALA):**

- Nakanishi et al. (2025, *Heart and Vessels*): Parent compound 5-ALA + ferrous iron significantly reduced myocardial infarct size in mice (IA/AAR: 34.0% vs 51.7% control, p=0.001) via HO-1 induction. HO-1 inhibitor ZnPP abolished protection.
- Caveat: HAL is the hexyl ester prodrug designed for bladder tissue penetration. Whether HAL can recapitulate 5-ALA's cardioprotection (via esterase cleavage) is unstudied.

**Anti-IRI (4★ Renal/Hepatic/Cerebral via 5-ALA):**

- **Renal**: Phase II RCT (TALEN trial, Oxford, 2024): 5-ALA-SFC well-tolerated in cardiac surgery patients with significant renoprotective effects.
- **Hepatic**: 5-ALA ameliorated IRI in fatty liver mice by reducing TLR2/4-NF-κB and ROS (Li et al., 2016); promoted M2 macrophage polarization (Jin et al., 2023).
- **Cerebral**: Intrastriatal 5-ALA microinjection significantly reduced MCAO infarct volumes in rats (2020).

**Verdict**: The 5-ALA backbone is powerful (heme → HO-1 → CO/biliverdin cytoprotection), and HAL is already FDA-approved for a different indication. Key gap: direct HAL cardioprotection studies needed.

---

### 9. MIGLUSTAT (麦格司他) — ⭐⭐⭐ Anti-MIRI / ⭐⭐ Anti-IRI

**Anti-MIRI (3★ Indirect Cardiac Protection):**

- Liu et al. (2025, *Molecular Medicine*): Miglustat significantly ameliorated isoproterenol-induced cardiac fibrosis in mice — improved LVEF/LVFS, reduced collagen deposition and fibrosis markers, reversed ERK/STAT3/Akt/GSK3β activation in cardiac fibroblasts.
- UGCG (glucosylceramide synthase) knockdown mimicked protective effects, validating the target.
- Caveat: This is a β-adrenergic overactivation model, not a direct ischemia-reperfusion model. No I/R study published.

**Anti-IRI (2★ Bioinformatics only):**

- Organ-protective effects demonstrated in non-IRI models; no direct IRI studies.
- Andersson et al. (2021, *Eur Heart J*): Cardiac UGCG knockout is detrimental (heart failure, impaired β1-AR trafficking) — suggesting pharmacological inhibition differs from genetic deletion.

**Verdict**: Interesting glycosphingolipid-targeting mechanism with cardiac antifibrotic effects. Needs dedicated MIRI model studies. Approved for Gaucher/NPC diseases → human safety data available.

---

### 10. PENICILLIN G BENZATHINE (苄星青霉素) — ⭐ Anti-MIRI / ⭐⭐⭐⭐ Anti-IRI

**Anti-MIRI (1★):**

- No evidence of direct cardioprotection. Black-box warning for inadvertent IV administration (cardiorespiratory arrest). No mechanistic rationale for direct anti-ischemic effect.

**Anti-IRI (4★ Intestinal):**

- Multiple studies from Yale University School of Medicine (2019 Academic Surgical Congress): Intraluminal penicillin G significantly reduced ischemic cellular injury in rat small bowel and colon segments via H,K-ATPase activation.
- Protection independent of NO pathway; reversed by H,K-ATPase inhibitor SCH-28080.
- CNS penetration extremely poor (CSF 1-2% of plasma) → no neuroprotection.

**Verdict**: Organ-specific intestinal IRI protection is well-documented. However, MIRI relevance is nil. Interesting for intestinal transplant/ischemia applications.

---

## Overall Assessment

### Summary Statistics

| Metric | Count |
|--------|-------|
| Drugs screened (Phase 1) | 50 |
| Passed Phase 1 | 25 (50%) |
| Excluded Phase 1 | 25 (50%) |
| ≥4★ Anti-MIRI | 6 (24% of PASS) |
| ≥4★ Anti-IRI | 9 (36% of PASS) |
| ≥4★ BOTH tracks | 5: Doxycycline, Ramatroban, Terazosin, Aprepitant, Riboflavin |

### Tiered Recommendations

**Tier 1 — Strongest Repurposing Candidates (≥4★ in both tracks, clinical data):**

| Drug | MIRI | IRI | Strengths |
|------|------|-----|-----------|
| **Doxycycline** | 5★ | 4★ | Clinical trial evidence (TIPTOP); broadest organ coverage; MMP mechanism well-characterized |
| **Tirofiban** | 4★ | 4★ | Already used in ACS/STEMI; 11 RCTs; RISK pathway activation |

**Tier 2 — Strong Preclinical Candidates (≥4★ in both tracks, no clinical MIRI trial):**

| Drug | MIRI | IRI | Strengths |
|------|------|-----|-----------|
| **Ramatroban** | 4★ | 4★ | Direct rat MIRI model; dual TP/DP2 mechanism; Japan-approved |
| **Terazosin** | 4★ | 4★ | PGK1 activation mechanism (Nature 2015); Phase 3 stroke trial recruiting |
| **Aprepitant** | 4★ | 4★ | PI3K-Akt pathway; extensive preclinical; established human safety |
| **Riboflavin** | 4★ | 4★ | Multi-organ preclinical; mitochondrial mechanism; vitamin safety profile |

**Tier 3 — Niche Candidates (high evidence in one track):**

| Drug | MIRI | IRI | Strengths |
|------|------|-----|-----------|
| **Ceftriaxone** | 3★ | 4★ | Cerebral IRI: 0% mortality vs 34.5%; GLT-1 mechanism |
| **Penicillin G** | 1★ | 4★ | Intestinal IRI via H,K-ATPase (niche application) |
| **Hexaminolevulinate** | 3★ | 4★ | 5-ALA backbone with Phase II renal trial |
| **Miglustat** | 3★ | 2★ | Novel glycosphingolipid mechanism; cardiac antifibrotic |

**Tier 4 — Low Priority (≤2★ in both tracks):**

14 drugs with minimal or no evidence of anti-MIRI/anti-IRI activity. These are predominantly anti-infective agents (cephalexin, dicloxacillin, pibrentasvir, relebactam, etc.) whose mechanisms are unrelated to ischemia-reperfusion pathophysiology.

### Key Patterns

1. **Anti-infectives dominate low-scoring drugs**: Most antibiotics and antivirals in the PASS list have no published IRI-protective evidence. Their inclusion in the MIRI pipeline likely reflects XPert-predicted transcriptional effects that may not translate to in vivo protection.

2. **MMP inhibition is the strongest common mechanism**: Doxycycline's MMP-2/9 inhibition is the most validated cardioprotective mechanism in our list, with clinical trial evidence.

3. **Anti-inflammatory/immunomodulatory drugs outperform**: Aprepitant (NK1), Ramatroban (TP/DP2), and Ceftriaxone (GLT-1) all target inflammation/excitotoxicity pathways that are central to IRI.

4. **Many drugs lack any IRI research**: 14/25 PASS drugs (56%) scored ★ in both Track A and Track B, representing a significant evidence gap between computational prediction and published validation.

### Limitations

- **Web search depth**: PubMed indexing may miss non-English literature, conference abstracts, or very recent preprints.
- **No full-text analysis**: Evidence grading is based on abstracts and key excerpts; full-text confirmation is recommended for Tier 1-2 candidates.
- **Search bias**: Well-studied drugs (doxycycline, tirofiban) are more likely to return positive results simply because they have been studied more extensively.
- **XPert model caveat**: The computational model predicts 97.2% down-regulated DEGs with low inter-drug variance. Many low-scoring drugs may represent false positives from the computational screen.

### Next Steps

1. **In vitro validation** of Tier 1-2 candidates in cardiomyocyte hypoxia/reoxygenation models
2. **Full-text review** of key papers for Doxycycline (TIPTOP trial methodology), Ramatroban (Squadrito 1993 dose-response)
3. **Mechanism-network mapping** for Terazosin PGK1 activation in cardiac ischemia context
4. **Clinical feasibility assessment** for Tirofiban (already used in ACS — could be repurposed for MIRI prevention during PCI)
5. **Patent landscape** analysis for Ramatroban and Terazosin repurposing

---

*Report generated via 20 parallel AI research agents. All cited findings should be verified against original publications before clinical decision-making.*

