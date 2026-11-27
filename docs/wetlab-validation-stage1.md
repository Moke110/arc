<!-- wetlab-validation-stage1
stage: 1
candidates: TRIMETHOPRIM (MP / 甲氧苄啶), MIRTAZAPINE (MT / 米氮平)
raw_data: wet_data/ (gitignored: stage1数据.pptx, stage1动物实验计划.xlsx)
-->

# Stage 1 wet-lab validation — trimethoprim (MP) & mirtazapine (MT)

First experimental (wet-lab) check of Stage-1 drug-repositioning candidates. Two candidates
were advanced from the in-silico MIRI screens — **trimethoprim (MP, 甲氧苄啶)** and
**mirtazapine (MT, 米氮平)** — and tested for (i) cardiomyocyte safety and hypoxia/reoxygenation
efficacy in vitro (CCK8) and (ii) safety and efficacy in a mouse ischaemia–reperfusion (IR)
model (%EF / %FS / % infarct area).

> **Headline.** Both drugs were non-toxic in vitro and in vivo. **Mirtazapine (MT) significantly
> reduced infarct area** in the mouse IR model (IR+MT vs IR, p<0.01); trimethoprim (MP) did not.
> Neither drug reached significance on echocardiographic function (%EF/%FS), and neither rescued
> cardiomyocyte viability in the in-vitro hypoxia/reoxygenation CCK8 assay. This is the first
> experimental readout for the Stage-1 screen: it **partially supports MT** and **does not support
> MP**.

The raw data live in `wet_data/` (gitignored): `stage1数据.pptx` (result figures) and
`stage1动物实验计划.xlsx` (animal-experiment plan). Values below are read from the slide
figures; no raw numeric tables or statistical output were provided.

## 1. Stage 1 provenance

Both drugs are among the 1,004 L1000-contained FDA drugs screened by all four Stage-1
workflows (v1 / v1.1 / v2 / v3). Their standing by the step-11 integrated MIRI-reversal score:

| workflow | TRIMETHOPRIM (MP) | MIRTAZAPINE (MT) |
|---|---|---|
| `xpert-fda-miri-v1` | rank 692/1004 (pct 0.69) | **rank 70/1004**; finalist **#15/42** (integrated +0.7867) |
| `xpert-fda-miri-v1.1` | **rank 68/1004** (pct 0.07) | rank 91/1004 (pct 0.09) |
| `xpert-fda-miri-v2` | **rank 93/1004** (pct 0.09) | **rank 46/1004**; finalist **#18/64** (integrated +0.9608) |
| `xpert-fda-miri-v3` | rank 856/1004 (pct 0.85) | rank 884/1004 (pct 0.88) |

MT was a **filtered finalist** in v1 (rank 15) and v2 (rank 18). MP was a **high-integrated-score
candidate** in v1.1/v2 (top ~7–9% by integrated score) but was excluded from the filtered final
list because its WCS was not < 0 (the finalist gate also requires ≥6/8 critical-pathway reversal).
The two candidates therefore represent two different Stage-1 signal classes: a filtered finalist
(MT) and an integrated-score near-miss (MP).

Note (v2 evidence survey): MT's dossier notes the human post-MI evidence is neutral and the only
direct preclinical cardiac result is a single mouse MI→HF cardioprotection study
(BBA Mol Basis Dis 2025, PMID 40398828; GRB14/RAB3A). No direct MIRI/cardioprotection study
exists for MP.

## 2. Candidate drugs — identity & cardiac-relevant biology

| | **Trimethoprim (MP)** | **Mirtazapine (MT)** |
|---|---|---|
| Chinese | 甲氧苄啶 | 米氮平 |
| Class | diaminopyrimidine antibacterial; DHFR inhibitor | tetracyclic antidepressant; NaSSA |
| ATC | J01EA01 | N06AX11 |
| Target / MoA | bacterial DHFR (folate → THF block); CYP2C8 off-target | α₂-adrenergic antagonist (↑NA/5-HT); 5-HT₂/5-HT₃ antagonist; H₁ antagonist |
| FDA | approved 1973 (now mainly co-trimoxazole TMP/SMX); UTI / bacterial infection | approved 1996; major depressive disorder |
| Human dose | TMP 100 mg q12h / 200 mg qd; co-trimoxazole 160/800 q12h | 15 mg qhs, 15–45 mg/day |
| MIRI/AMI evidence | **no direct study found**; only an in-silico drug–disease link (PMID 35716513) | mouse MI→HF cardioprotection (PMID 40398828); rat renal IRI (PMID 22032749) |
| CV safety | ENaC-mediated **hyperkalaemia**; ↑sudden death with RAAS inhibitors (PMID 8328738, 10738846, 25359996) | small dose-dependent QTc ↑ (~+2.4 ms at 45 mg, PMID 30052325); "possible" TdP class; possible MACE signal vs SSRI (PMID 30900189) |

Mechanistic notes for the MIRI hypothesis: MT has a plausible **anti-oxidant/anti-inflammatory**
profile (↑GSH/SOD, ↓MDA/MPO; ↓NF-κB/TLR4, ↑Nrf2 across tissues) and a concrete cardiac axis
(GRB14/RAB3A). MP's MIRI rationale is weak — DHFR selectivity is bacterial, its main CV action
(hyperkalaemia) is a liability, and any benefit would be indirect (anti-inflammatory/lipid effects
of co-trimoxazole, PMID 41029321).

## 3. Experimental design

### 3.1 In vitro — cardiomyocytes, CCK8 viability
- **Safety:** vehicle vs MP / MT at 10 nM, 100 nM, 1 µM, 10 µM; viability normalised to vehicle; n=3.
- **Efficacy:** vehicle vs MP / MT at 5, 10, 20 µM under **normoxia (Nor)** vs **hypoxia/reoxygenation
  (“H2R2”, inferred as 2 h hypoxia + 2 h reoxygenation)**; viability normalised to the normoxia
  vehicle; n=3.

### 3.2 In vivo — mouse model (from `stage1动物实验计划.xlsx`)

| group | dosing | n | freq | days | vol | conc | BW | per-dose (derived) |
|---|---|---:|---|---:|---:|---:|---:|---:|
| Sham | vehicle | 3 | 1×/d | 3 | 10 mL/kg | — | 0.02 kg | — |
| MP | MP | 3 | 2×/d | 3 | 5 mL/kg | 5 mg/mL | 0.02 kg | ≈25 mg/kg/dose (≈50 mg/kg/day) |
| MT | MT | 3 | 1×/d | 3 | 10 mL/kg | 0.75 mg/mL | 0.02 kg | ≈7.5 mg/kg/dose (≈7.5 mg/kg/day) |
| IR | vehicle | 4 | 1×/d | 3 | 10 mL/kg | — | 0.02 kg | — |
| IR+MP | MP | 4 | 2×/d | 3 | 5 mL/kg | 5 mg/mL | 0.02 kg | ≈25 mg/kg/dose |
| IR+MT | MT | 4 | 1×/d | 3 | 10 mL/kg | 0.75 mg/mL | 0.02 kg | ≈7.5 mg/kg/dose |

- **Safety arm:** Sham / MP / MT in non-IR animals → cardiac function (%EF, %FS).
- **Efficacy arm:** Sham / IR / IR+MP / IR+MT → %EF, %FS, % infarct area.
- “Derived” dose = concentration × volume-per-kg; the plan sheet does not state mg/kg directly.

## 4. Results

### 4.1 In vitro
| assay | readout | observation |
|---|---|---|
| safety | relative viability (vehicle = 1.0) | MP and MT, 10 nM–10 µM: no reduction; viability flat-to-higher (MP 100 nM–1 µM ≈1.6–2.1; MT 100 nM–1 µM ≈1.6–1.9). **No cytotoxicity.** |
| efficacy | normalised viability | H/R lowered viability to ≈0.5 vs ≈1.0 normoxia (p<0.001). MP/MT at 5/10/20 µM trended up (≈0.65–0.83) but were **not significant vs H/R vehicle (ns)**. **No significant in-vitro rescue.** |

Relaxed (non-normalised) CCK8 plots (slide 1) also show the MP/MT dose series rising above the
vehicle column, consistent with a growth/viability-promoting rather than toxic effect at these
concentrations.

### 4.2 In vivo
Approximate values read from the figures (mean; n as above):

| readout | Sham | IR | IR+MP | IR+MT | significance |
|---|---:|---:|---:|---:|---|
| %EF | ≈78 | ≈32 | ≈50 | ≈54 | Sham vs IR **; IR vs MP ns; IR vs MT ns |
| %FS | ≈43 | ≈15 | ≈25 | ≈27 | Sham vs IR **; IR vs MP ns; IR vs MT ns |
| % infarct area | ≈1 | ≈18 | ≈18 | ≈12.5 | Sham vs IR ****; **IR vs MP ns**; **IR vs MT \*\*** |
| %EF (safety arm) | Con ≈76 | — | MP ≈76 | MT ≈68 | no drop vs control |
| %FS (safety arm) | Con ≈43 | — | MP ≈44 | MT ≈37 | no drop vs control |

**Verdict.** MT produced a significant reduction in infarct area (**, p<0.01) and consistent
(non-significant) upward trends in %EF/%FS; it was non-toxic in non-IR animals. MP produced no
significant change on any in-vivo readout and was non-toxic.

## 5. Interpretation

- **MT (mirtazapine):** the animal efficacy endpoint supports the Stage-1 prediction — significant
  infarct-size reduction with a safety signal. Consistent with the single prior mouse MI→HF study
  (PMID 40398828). This is the first wet-lab readout that **partially supports** a Stage-1 candidate.
- **MP (trimethoprim):** the animal efficacy endpoint does **not** support the Stage-1 prediction —
  no significant effect on infarct or function, despite being safe. The integrated-score signal was
  not borne out at the infarct endpoint.
- **In-vitro vs in-vivo discordance:** the H/R CCK8 viability assay detected no rescue by either
  drug, while MT reduced infarct area in vivo. CCK8 viability is a coarse, single-endpoint readout;
  it may be insensitive to the protective mechanism, the H/R window may be too short/severe, and/or
  the in-vitro doses (5–20 µM) may not match the effective in-vivo exposure. The two arms are not
  contradictory, but the in-vitro assay currently does not track the in-vivo benefit.

### 中文摘要
首次对 Stage-1 候选药 **甲氧苄啶 (MP)** 与 **米氮平 (MT)** 做湿实验验证。细胞 CCK8：两药
10 nM–10 µM 均无毒性；缺氧/复氧 (H2R2) 模型中 5–20 µM 有上升趋势但无统计学差异 (ns)。
小鼠 IR 模型：两药在非 IR 动物中不影响心功能（安全）；**MT 显著减小梗死面积**（IR+MT vs IR,
p<0.01），MP 无显著作用；两者对 %EF/%FS 的改善均未达显著。结论：**MT 得到部分支持，MP
未获支持**，证明回归模型能筛出有体内获益信号的候选药，但样本量小、单剂量、且体外 CCK8
未能复现 MT 的体内获益。

## 6. Limitations

- Small animal n (3–4/group); echo endpoints underpowered → ns despite directional trends.
- **Figure-read values only** — no raw numeric tables, no full statistical output were supplied.
- `H2R2` decoded as hypoxia/reoxygenation by inference (paired with `Nor` = normoxia); not confirmed
  against a protocol sheet.
- Single in-vivo dose per drug; no dose–response, no mechanism, no biomarker, no blinding/randomisation
  documented.
- MP was not a Stage-1 filtered finalist; the Stage-1 provenance differs between the two candidates
  (see §1).
- Drug biology/citations in §2 are from a literature/database briefing, not from this experiment.

## 7. Data & files

- Raw results figures: `wet_data/stage1数据.pptx` (gitignored).
- Animal-experiment plan: `wet_data/stage1动物实验计划.xlsx` (gitignored).
- This record: `docs/wetlab-validation-stage1.md`.

## References

Drug biology / safety (literature & database briefing; see `reference/` for the background library):
ChEMBL CHEMBL22 / CHEMBL654; KEGG D00145 / D00563; PMIDs 8328738, 10738846, 25359996, 41029321
(trimethoprim); 40398828, 22032749, 30052325, 30900189 (mirtazapine). In-silico MIRI link for
trimethoprim: PMID 35716513. Full mirtazapine evidence dossier:
`stage1/12_drug_survey/drugs/MIRTAZAPINE.md`.
