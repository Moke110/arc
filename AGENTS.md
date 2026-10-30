# AGENTS.md

Routing index for the ARC repo root: AI-driven MIRI drug-repositioning pipeline (SMILES → HG/UniMol embeddings → XPert DEG prediction → efficacy ranking). Research code, not a software project.

- `pert_data/` — CM drug-perturbation data lake and XPert pretrain h5ads (gitignored; own AGENTS.md).
- `drug_emb_model/` — HG HeteroGNN drug-embedding model library (pure code).
- `pert_model/` — XPert drug-perturbation prediction library; vendored candidates: deployed `prnet`, `cpa`, `chemcpa`, `biolord`, `perturbnet`, `scvidr`, `squidiff`, `cmonge`, `decant`; reference-only `wave`, `expo`, `state`, `TranSiGen`, `xpert`.
- `stage1/` — Stage 1 MIRI workflows: frozen v0, finalized v1/v1.1/v2/v3 screens.
- `stage2/` — Stage 2 model setup and selection: `01_setup/` run cards for the nine deployed models, `02_cm_data/` (frozen CM corpus + `datasets-1/`, ADR 0026), `03_models_cmbench/` (executed 60-unit cold-drug benchmark as two substeps: `03.1_cold_drug_18k/` primary, `03.2_hvg_sensitivity_5k/` sensitivity fork), `_envs/` (shared per-model bench environments), `audit/`, and `_archive/cm_bench_v1/` (old benchmark).
- `disease_data/` — shared per-GSE MIRI disease inputs (raw data + DEG lists).
- `wet_data/` — Stage-1 wet-lab validation raw data (result figures pptx + animal-experiment plan; gitignored); record in `docs/wetlab-validation-stage1.md`.
- `reference/` — background literature library (PDFs, markdowns, index; gitignored).
- `presentation/` — thesis proposal deck (`开题.pptx`, 20 slides + notes) and its background literature review (`my_review.md`).
- `docs/` — ADRs, the AICM roadmap (`roadmap.md`), the wet-lab record (`wetlab-validation-stage1.md`), agent docs, and teaching material.
- `archived/` — gitignored provenance archive of orphaned HG FDA-extension data/weights.
- `logs/` — runtime logs (e.g. unimol_tools).
- `GLOSSARY.md` — domain glossary: shared vocabulary and category taxonomy.
- `environment.yml` — shared conda environment `arc`.
- `skills-lock.json` — pinned agent-skills lockfile.
- `AGENTS.md` — this routing index; details live in subdirectory AGENTS.md files.
- `.agents/`, `.claude/` — agent skills and local config.
- `.scratch/`, `.uploads/`, `.vscode/`, `.serena/`, `.env`, `.gitignore` — workspace scratch, tooling, and config.

## Agent skills

### Issue tracker

Issues and specs live as local markdown files under `.scratch/<feature-slug>/`. See `docs/agents/issue-tracker.md`.

### Triage labels

Default canonical labels (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: root `GLOSSARY.md` + `docs/adr/`. See `docs/agents/domain.md`.
