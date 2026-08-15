# Per-drug fact-sheet agent query (12_report)

This template is loaded by `run_drug_agents.py`, which substitutes `{DRUG}` and
`{SMILES}` and sends the result to one biomni A1 agent per drug. The agent must
output a `<solution>` block with the exact field-prefixed lines below so the
assemblers can parse it.

---

Investigate the drug {DRUG} (canonical SMILES: {SMILES}) for a drug-repositioning
fact sheet on MIRI (myocardial ischemia-reperfusion injury) / AMI (acute
myocardial infarction). Use the biomni tools (KEGG drug, ChEMBL, openFDA,
OpenTargets, PubMed/literature, TxGNN repurposing) to gather evidence for each
field. If a field cannot be determined, write 'not found' rather than guessing.

Determine exactly these fields:
- ENGLISH_NAME: English name (with common brand names)
- CHINESE_NAME: Chinese name (中文名, standard transliteration if any)
- US_FDA_APPROVAL: US FDA approval status (approved / withdrawn / not found; year and indication if known)
- CHINA_NMPA_APPROVAL: China NMPA approval status (approved / not found)
- CLASSIFICATION: Pharmacological classification (drug class, mechanism, ATC if known)
- CLINICAL_SCOPE: Clinical application scope (approved indications and main clinical uses)
- MIRI_AMI_CLINICAL: Connection to MIRI/AMI: clinical evidence (trials, post-MI use, epidemiology)
- MIRI_AMI_PRECLINICAL: Connection to MIRI/AMI: preclinical evidence (animal models, mechanism, cardioprotection)

For the two MIRI/AMI fields you MUST make a genuine effort: run at least one
PubMed/literature search and one TxGNN repurposing prediction for
MIRI/AMI/ischemia-reperfusion, and report what you find (a negative result is
fine — report it as the evidence you actually found), before writing 'not found'.

After gathering evidence, output a <solution> block that is a plain fact sheet
with exactly these lines, one field per line, prefixing each value with the
field key and ': ' (keep values on a single line; do not use markdown lists
inside a field):
ENGLISH_NAME:
CHINESE_NAME:
US_FDA_APPROVAL:
CHINA_NMPA_APPROVAL:
CLASSIFICATION:
CLINICAL_SCOPE:
MIRI_AMI_CLINICAL:
MIRI_AMI_PRECLINICAL:
