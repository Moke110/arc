# AGENTS.md

Terse index of `stage1/xpert-fda-miri-v1.1/` (DToxS-LINCS testbed, an 11-step mirror of v1).

- `screen.ipynb` — notebook view of the screen; regenerate via `_lib/build_screen_notebook.py`.
- `01_drug_embed/` — drug list + author-reference HG embeddings.
- `02_cm_baseline/` — DToxS-primary consensus baseline.
- `03_cm_finetune_data/` — rebuilds the DToxS-LINCS finetune h5ad, parity delta, composition report.
- `04_xpert_pretrain/` — carried-over L1000 pretrain weight + config.
- `05_xpert_finetune_bench/` — unified warm-start/cold-drug folds + finetune benchmark.
- `06_xpert_finetune/` — all-train DToxS-LINCS finetune producing the prediction model.
- `07_xpert_predict/` — perturbed-transcriptome prediction + z-scoring.
- `08_miri_phenotype/` — MIRI reference gene sets + rank reference.
- `09_miri_efficacy_scores/` — efficacy scorers, z-components, plots.
- `10_filter/` — WCS and critical-pathway filters.
- `11_integrate/` — integrated score, ranked candidates, deterministic report.
- `input/` — frozen assets: HG reference, UniMol array, idx→SMILES, FDA list, MIRI pathways.
- `_lib/` — shared helpers: paths, gene mapping, drug resolution, scoring, gene gates, notebook builder.
- `AGENTS.md` — this index.
