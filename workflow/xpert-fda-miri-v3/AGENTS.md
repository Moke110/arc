# AGENTS.md

Routing index for `workflow/xpert-fda-miri-v3/` — all-CM finetune sibling of v1.

- `screen.ipynb` — top-level screening notebook.
- `01_drug_embed/` — merged L1000-contained drug list + author-reference HG embeddings.
- `02_cm_baseline/` — DToxS-primary control consensus baseline.
- `03_cm_finetune_data/` — builds per-dataset and union finetune h5ads, reports, schema gate.
- `04_xpert_pretrain/` — frozen L1000 pretrain weight + config.
- `05_xpert_finetune_bench/` — benchmarks; `dtoxs-primary/` official folds, `data-ablation/` scale study.
- `06_xpert_finetune/` — all-train union finetune to the final model.
- `07_xpert_predict/` — perturbed-transcriptome prediction + z-scores.
- `08_miri_phenotype/` — MIRI reference gene sets and rank reference.
- `09_miri_efficacy_scores/` — efficacy scorers, z-components, plots.
- `10_filter/` — critical-pathway filters.
- `11_integrate/` — integrated scores, ranked candidates, report.
- `input/` — frozen drug/gene assets (HG reference, UniMol, SMILES, pathways).
- `_lib/` — shared path utils, gene mapping, drug resolution, gene gates.
- `config/` — empty placeholder directory.
