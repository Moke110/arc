# Step-12 per-drug evidence-survey brief (12_drug_survey)

This template is the instruction given to one biomni agent per drug from the
top-20 candidates of `11_integrate/output/ranked_drugs.csv` (step-11 ranking,
integrated score descending). The agent researches the drug with the biomni
tool suite (PubMed/literature, ClinicalTrials.gov, openFDA/DailyMed, ChEMBL,
KEGG drug, OpenTargets, TxGNN repurposing, web search), proposes two evidence
scores, and returns a structured bilingual dossier. The step owner then
sanity-checks the proposed scores against the cited evidence, possibly adjusts
them, and writes the dossier to `drugs/{DRUG}.md`. Scoring is deliberately
evidence-based — report honestly what is and is not found; never fabricate a
report or a citation.

---

## Drug under investigation

- Drug name (author L1000 / INN): **{DRUG}**
- Canonical SMILES: `{SMILES}`
- Step-11 rank: **{RANK}** of 20 (integrated score {INTEGRATED})

## Evidence-score scale (apply to BOTH report sets below)

Two independent scores per drug, each drawn from this rubric:

| Score | Meaning |
|---:|---|
| +2 | 临床有效 — human clinical evidence shows efficacy/benefit (trials, guideline-endorsed post-MI/IRI use, outcome improvement). |
| +1 | 临床无报道且临床前有效 — no human clinical report, but positive preclinical evidence (animal models / mechanism / cardioprotection). |
| 0 | 完全无报道，或临床前报道之间存在矛盾 — nothing found at all, or preclinical reports conflict / no clear benefit. |
| −1 | 临床无报道且临床前无效 — no human clinical report, and the preclinical evidence is negative (no benefit / worsening in models). |
| −2 | 临床无效/有害 — clinical evidence shows ineffectiveness or harm (contraindicated, no outcome benefit, or increased mortality/events in humans). |

## Report sets to research

### A. 心梗 / 心肌缺血再灌注 (MI / MIRI)
- **临床报道** (clinical reports): trials, post-MI / peri-PCI use, epidemiology, guideline status, contraindication warnings.
- **临床前报道** (preclinical reports): animal MI / coronary occlusion-reperfusion / Langendorff models, mechanism of cardioprotection.
- Then: **proposed `score_miri`** with a short justification referencing the actual findings.

### B. 其他缺血再灌注类疾病 (other ischemia-reperfusion diseases)
Cerebral / renal / hepatic / intestinal / skeletal-muscle / pulmonary / spinal / retinal IRI or transplantation-reperfusion, stroke-reperfusion, etc.
- **临床报道** and **临床前报道** as above (explicitly name the disease/organ model).
- Then: **proposed `score_iri`** with justification.

Every individual report you cite must carry its source **inline at the point it
is mentioned** — journal name + PMID (or DOI/URL/ClinicalTrials.gov ID). If a
search genuinely returns nothing for a section, write "no evidence found"
rather than guessing.

## Required fields in the returned dossier

1. **中文名** — standard Chinese name (药典/transliteration if any).
2. **NMPA 审批** — China NMPA approval status (approved / not found / withdrawn; indication if known).
3. **FDA 审批** — US FDA approval status (approved / not found / withdrawn; year + indication if known).
4. **药理学分类** — pharmacological class, mechanism, ATC if known.
5. **A. MI/MIRI 临床 + 临床前** (each with inline sources), takeaway line.
6. **B. 其他缺血再灌注 临床 + 临床前** (each with inline sources), takeaway line.
7. **proposed `score_miri`** (integer −2…+2) + rationale.
8. **proposed `score_iri`** (integer −2…+2) + rationale.

Return the dossier as bilingual markdown (English + 中文), with the inline
citations preserved so the step owner can verify each claim.
