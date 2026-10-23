# CPA fidelity ledger (ticket 13)

> Derived from `stage2/03_bench_models/temp_env_all_18k/models/cpa/FIDELITY.md`; the gene-space rows are re-declared for the 5,000-HVG sensitivity axis (ADR 0028).

| component | author setting | run setting | class | reason | expected effect | ticket |
|---|---|---|---|---|---|---|
| gene space | 5,000 HVG (Kang/Norman) / 978 LINCS | 5,000 rank-sum HVG (corpus subset) | recipe_substitution | width matches the authors; the benchmark fixes one cross-view axis (ADR 0028) | small: same width, different HVG statistic and gene identity | 09 |
| corpus scale | 1,000s–100,000s of cells | 3,600 wells/fold, 276 train drugs | recipe_substitution | curated CM view is the benchmark | negative: fewer cells and contexts | 09 |
| initialisation | released Kang checkpoint for the setup smoke; from scratch for new data | random init, no released weights | recipe_omission | ticket 09 forbids public-weight fine-tuning | negative at fixed epochs | 09 |
| HVG selection | `highly_variable_genes(5000)` on the whole object | cross-view rank-sum HVG (5,000) | recipe_substitution | ADR 0028 replaces the per-object Seurat HVG with a fixed multi-source rank sum | small: both keep 5,000 genes; only the selection rule differs | 09 |
| reconstruction | Gaussian NLL on author-normalised data (paper) or NB on counts (tutorial) | Gaussian NLL on `log1p_cp10k` | recipe_substitution | benchmark Δ space is cp10k log1p | small; same loss family as the paper | 09 |
| cold-drug encoder | `use_rdkit_embeddings` frozen 2,048-bit Morgan + shared Linear | identical, fingerprints from parent-form SMILES | adapter_extension | CPA's OOD path, enabled at the adapter seam | enables cold drugs; no score penalty expected | 09 |
| cold-drug dose scaler | per-drug `beta`/`bias` learned | unseen drug's row at `β=1`, `bias=0` | adapter_extension | no treated rows exist for test drugs | negative: dose response for cold drugs is not calibrated | 09 |
| setup data | AnnData of the split under study | full fold AnnData (train+valid+test metadata + controls) | adapter_extension | gives every test drug an embedding row | none: test rows are never optimised or validated | 13 |
| controls | ordinary training rows, never held out | 272 DMSO cells in the training split | adapter_extension | controls are shared across folds by contract | none expected | 04 |
| epochs | author heuristic 400 for small data; tutorial 2,000 | explicit 2,000 | config_only | ticket 09 fixes the author budget | none if early stopping fires first | 09 |
| adversarial plan | 3 Adam, lr 5e-4/3e-4/3e-4, wd 1e-6/4e-7/4e-7, StepLR(45, 0.9) | identical | config_only | author defaults reused | none | 09 |
| architecture | 512x3 encoder/decoder, 128 latent, dropout 0.1, logsigm doser | identical | config_only | author tutorial settings | none | 09 |
| prediction | cell-level `CPA.predict`, mean over samples | basal = matched plate DMSO mean, one deterministic sample | config_only | one-shot test scoring | slightly lower variance; no expected bias | 09 |
