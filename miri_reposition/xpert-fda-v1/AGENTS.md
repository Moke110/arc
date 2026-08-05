# AGENTS.md

Active MIRI drug repositioning workflow. Ranks 3,127 FDA-approved small-molecule drugs by predicted reversal of the myocardial ischemia-reperfusion injury (MIRI) transcriptional signature. Refactored successor of the frozen `miri_reposition/xpert-fda-v1_history/` (gitignored; parity reference). Cardiotoxicity scoring deferred from v1.

## Structure

Step-wise: each numbered step dir has its script, `output/` artifacts, and plots, reading previous steps' `output/` plus static files in `input/`.

```
01_drug_list/          # 3,127-drug FDA list -> output/drug_list.csv
02_drug_embeddings/    # verifies frozen inputs; HG inference -> output/all_drugs_idx2HG.npy
03_baseline_transcriptome/  # DTOXS consensus baseline -> output/cm_ctl_landmark.csv
04_finetune_xpert/     # full-DToxS finetune (final model) -> output/dtoxs_dtoxs/
05_benchmark_xpert/    # split finetunes + eval -> output/benchmark.csv + benchmark.png
06_perturbed_transcriptome/  # XPert predict -> output/xpert_fda-cm_deg.csv
07_zscore/             # per-gene z-score -> output/xpert_fda-cm_deg-zs.csv
08_miri_ref_path/      # integrated MIRI reference -> output/MIRI_{up-down,rank}.json
09_efficacy_scores/    # 3 scorers -> output/scores.csv + plots
10_filter/             # WCS<0 + >=6/8 critical pathways -> output/ + pyramid png
11_integrate/          # integrated score, final rank -> output/report.md + heatmaps
```

`_lib/` holds shared helpers: `utils.py` (path resolution + DEG/z-score), `scoring.py` (WCS+RGES, Zhang/ssCMap, pathway RGES, integration), `verify_parity.py`, `build_screen_notebook.py`. `screen.ipynb` orchestrates all steps; `input/` holds frozen external inputs (dtoxs.h5ad, UniMol/HG embeddings, MIRI_pathways.json).

## Steps

- 01–03: drug list, embedding verification + HG inference, baseline transcriptome — all **idempotent** (`--force` to regenerate).
- 04 finetune XPert on full DToxS (final model). Idempotent (`--force` to retrain). Pins history hyperparameters (batch 32, init_epoch 0, finetune_lr 0.002). XPert stages into `<workflow>/finetune/output/`; step moves results to its own `output/`. For bit-exact weights set `XPERT_DETERMINISTIC=1` + `CUBLAS_WORKSPACE_CONFIG=:4096:8`.
- 05 benchmark on splits (warm-start, cold-drug, cold-cate). Idempotent per split (`--force` retrain+reeval, `--force-eval` re-evaluate only). Downstream uses the step-04 model, not benchmark weights.
- 06 predict perturbed DEG for all FDA drugs. Consumes step-04 model; indexes the precomputed UniMol array by drug index (feature consistency). Idempotent (`--force`). Current table is the archived Jun-10 prediction (byte-identical to history).
- 07 z-score the predicted DEG (per-gene, across-drug) before scoring — always use it.
- 08 MIRI reference from `disease_data/miri/{GSE}/deg/`; byte-identical to the frozen reference.
- 09 three scorers → z-components for all drugs (does NOT compute integrated score).
- 10 filter WCS<0 + ≥6/8 critical pathways.
- 11 integrated score over ALL drugs (z-components are population stats), rank filtered candidates.

## Key files

- `_lib/scoring.py` — `merge_components` (step 09) + `add_integrated` (step 11); byte-identical to frozen history scripts on identical inputs.
- `input/all_drugs_unimol_arr.npy` — frozen UniMol features; finetune and predict consume it by drug index.

## Dependencies

- Model libraries: `pert_model/xpert/` (XPert), `drug_emb_model/hg/` (HG inference). Disease inputs: `disease_data/miri/{GSE}/deg/`. Env `arc`.

## Entry points

```bash
cd miri_reposition/xpert-fda-v1/ && conda activate arc
python 01_drug_list/01_drug_list.py
# ... each step in order ...
python 11_integrate/11_integrate.py
# or run all in order:
jupyter notebook screen.ipynb
# parity vs the frozen pipeline:
python _lib/verify_parity.py
```

## Notes

- Heavy steps (02, 04, 05, 06) require a GPU.
- The frozen `xpert-fda-v1_history/` snapshot is the parity reference for this refactor.
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- The 5.2 GiB float64 UniMol array is memory-mapped by the XPert loaders (see `pert_model/xpert/AGENTS.md`), so concurrent finetunes no longer OOM WSL (per-finetune peak RSS ~2.1 GB).
