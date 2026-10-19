# AGENTS.md

Routing index for `stage1/xpert-fda-miri-v1/`, the canonical 11-step MIRI FDA screening workflow (author-aligned, L1000-contained screen).

- `01_drug_embed/` — merged L1000 drug list + author-reference HG embeddings.
- `02_cm_baseline/` — consensus cardiomyocyte baseline with authoritative gene mapping.
- `03_cm_finetune_data/` — workflow-owned frozen DToxS finetune h5ad.
- `04_xpert_pretrain/` — workflow-owned XPert L1000 pretrain weight + config.
- `05_xpert_finetune_bench/` — DTOXS fold generation + multi-fold finetune benchmark.
- `06_xpert_finetune/` — all-train DTOXS XPert finetune + config.
- `07_xpert_predict/` — perturbed-transcriptome prediction and z-scoring.
- `08_miri_phenotype/` — MIRI reference gene sets + rank reference.
- `09_miri_efficacy_scores/` — three MIRI efficacy scorers, z-components, plots.
- `10_filter/` — WCS and critical-pathway candidate filters.
- `11_integrate/` — integrated ranking over all drugs + deterministic report.
- `report/` — archived agentic bilingual final drug report stage.
- `input/` — frozen drug arrays, FDA list, MIRI pathways, author HG reference.
- `_lib/` — shared paths, scoring, gene mapping, gene gates, notebook builder.
- `screen.ipynb` — assembled screening notebook.
- `AGENTS.md` — this routing index.
