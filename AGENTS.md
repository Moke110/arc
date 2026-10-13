# AGENTS.md

Routing index for the ARC repo root: AI-driven MIRI drug-repositioning pipeline (SMILES → HG/UniMol embeddings → XPert DEG prediction → efficacy ranking). Research code, not a software project.

- `pert_data/` — CM drug-perturbation data lake and XPert pretrain h5ads (gitignored; own AGENTS.md).
- `drug_emb_model/` — HG HeteroGNN drug-embedding model library (pure code).
- `pert_model/` — XPert drug-perturbation prediction library; plus non-implemented new-model design docs.
- `workflow/` — MIRI workflows: frozen v0, finalized v1/v1.1/v2/v3 screens, plus bench HTML; and the Stage 2 whole-transcriptome benchmark `workflow/cm_bench_v1/`.
- `disease_data/` — shared per-GSE MIRI disease inputs (raw data + DEG lists).
- `wet_data/` — Stage-1 wet-lab validation raw data (result figures pptx + animal-experiment plan; gitignored); record in `docs/wetlab-validation-stage1.md`.
- `reference/` — background literature library (PDFs, markdowns, index; gitignored).
- `presentation/` — thesis proposal deck (`开题.pptx`, 20 slides + notes) and its background literature review (`my_review.md`).
- `docs/` — ADRs, the AICM roadmap (`roadmap.md`), the wet-lab record (`wetlab-validation-stage1.md`), agent docs, and teaching material.
- `archived/` — gitignored provenance archive of orphaned HG FDA-extension data/weights.
- `logs/` — runtime logs (e.g. unimol_tools).
- `CONTEXT.md` — workspace overview and category taxonomy.
- `environment.yml` — shared conda environment `arc`.
- `skills-lock.json` — pinned agent-skills lockfile.
- `AGENTS.md` — this routing index; details live in subdirectory AGENTS.md files.
- `.agents/`, `.claude/` — agent skills and local config.
- `.scratch/`, `.uploads/`, `.vscode/`, `.serena/`, `.env`, `.gitignore` — workspace scratch, tooling, and config.
