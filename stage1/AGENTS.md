# AGENTS.md

Routing index for the retained Stage 1 union-CM MIRI screen (ADR 0105). See README.md for reproduction.

- `screen.ipynb` — end-to-end notebook walkthrough of the 11-step screen.
- `01_drug_embed/` — merged L1000 drug list plus author-reference HG embeddings.
- `02_cm_baseline/` — control consensus baseline for the screen.
- `03_cm_finetune_data/` — builds the union finetune h5ad, gene maps, manifest, reports, gates.
- `04_xpert_pretrain/` — workflow-owned L1000 pretrain weight and config.
- `05_xpert_finetune_bench/` — three finetune benchmarks: `full/` (official folds), `data-ablation/`, and `warmstart_pr/` (direction-aware PR curves + AUPRC on the official warm-start folds).
- `06_xpert_finetune/` — all-train union finetune producing the final prediction model.
- `07_xpert_predict/` — perturbed-transcriptome prediction and z-scoring.
- `08_miri_phenotype/` — integrated MIRI reference gene sets and rank reference.
- `09_miri_efficacy_scores/` — three efficacy scorers with z-components and plots.
- `10_filter/` — critical-pathway and WCS filters over drug candidates.
- `11_integrate/` — integrated score, ranked candidates, deterministic report.
- `12_drug_survey/` — agentic bilingual evidence survey of the top-20 candidates.
- `13_limits/` — read-only limitations analysis: landmark coverage, critical-pathway overlap, screened-drug embedding PCA (+ `report.md`).
- `_other/` — ad-hoc analyses outside the numbered screen (GO study, structure gallery).
- `_lib/` — shared utils, gene mapping, drug resolution, scoring, verification helpers.
- `input/` — frozen screen assets: HG reference, UniMol array, SMILES, drug list, pathways, and predecessor schema/fold H5ADs in references/.
- `AGENTS.md` — this routing index.
