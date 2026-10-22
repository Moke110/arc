# CPA fidelity ledger (ticket 13)

| component | author setting | run setting | class | reason | expected effect | ticket |
|---|---|---|---|---|---|---|
| gene space | 5,000 HVG (Kang/Norman) / 978 LINCS | 18,164 corpus genes | recipe_substitution | fixed benchmark gene universe; no HVG selection | large negative: 3.6x–18x wider decoder | 09 |
| corpus scale | 1,000s–100,000s of cells | 3,600 wells/fold, 276 train drugs | recipe_substitution | curated CM view is the benchmark | negative: fewer cells and contexts | 09 |
| initialisation | released Kang checkpoint for the setup smoke; from scratch for new data | random init, no released weights | recipe_omission | ticket 09 forbids public-weight fine-tuning | negative at fixed epochs | 09 |
| HVG selection | `highly_variable_genes(5000)` on the whole object | none | recipe_omission | ticket 09 fixes all 18,164 genes | negative: noise genes in the decoder | 09 |
| reconstruction | Gaussian NLL on author-normalised data (paper) or NB on counts (tutorial) | Gaussian NLL on `log1p_cp10k` | recipe_substitution | benchmark Δ space is cp10k log1p | small; same loss family as the paper | 09 |
| cold-drug encoder | `use_rdkit_embeddings` frozen 2,048-bit Morgan + shared Linear | identical, fingerprints from parent-form SMILES | adapter_extension | CPA's OOD path, enabled at the adapter seam | enables cold drugs; no score penalty expected | 09 |
| cold-drug dose scaler | per-drug `beta`/`bias` learned | unseen drug's row at `β=1`, `bias=0` | adapter_extension | no treated rows exist for test drugs | negative: dose response for cold drugs is not calibrated | 09 |
| setup data | AnnData of the split under study | full fold AnnData (train+valid+test metadata + controls) | adapter_extension | gives every test drug an embedding row | none: test rows are never optimised or validated | 13 |
| controls | ordinary training rows, never held out | 272 DMSO cells in the training split | adapter_extension | controls are shared across folds by contract | none expected | 04 |
| epochs | author heuristic 400 for small data; tutorial 2,000 | explicit 2,000 | config_only | ticket 09 fixes the author budget | none if early stopping fires first | 09 |
| adversarial plan | 3 Adam, lr 5e-4/3e-4/3e-4, wd 1e-6/4e-7/4e-7, StepLR(45, 0.9) | identical | config_only | author defaults reused | none | 09 |
| architecture | 512x3 encoder/decoder, 128 latent, dropout 0.1, logsigm doser | identical | config_only | author tutorial settings | none | 09 |
| prediction | cell-level `CPA.predict`, mean over samples | basal = matched plate DMSO mean, one deterministic sample | config_only | one-shot test scoring | slightly lower variance; no expected bias | 09 |
