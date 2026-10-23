# Squidiff fidelity ledger (ticket 15)

> Derived from `stage2/03_bench_models/temp_env_all_18k/models/squidiff/FIDELITY.md`; the gene-space rows are re-declared for the 5,000-HVG sensitivity axis (ADR 0028).

Classes follow ticket 09: `config_only`, `adapter_extension`,
`recipe_omission`, `recipe_substitution`. Every row states an expected effect;
the target is zero *undeclared* deviations.

Author settings are cited to the vendored source in
`.scratch/stage2-03-bench-models/research/03-squidiff-cmonge-decant.md` §1
(`pert_model/squidiff` @ `abdfc27d`).

| component | author setting | run setting | class | reason | expected effect | ticket |
|---|---|---|---|---|---|---|
| gene space | 200 HVGs (authors' SciPlex run; `prep_sciplex_data.ipynb` cell 2) | 5,000 rank-sum HVG (`--gene_size 5000 --output_dim 5000`) | recipe_substitution | fixed benchmark gene universe (ADR 0028); the code has no cap (`MLPModel.py:134-135`) | large negative; 25x gene-space extrapolation (5,000 / 200), never validated by the authors above 200 genes | 09 |
| corpus scale | full SciPlex atlas (single-cell, all compounds/doses) | 2,523 train wells = 276 drugs x 3 doses per fold, 846 cold-valid wells | recipe_substitution | the curated CM view is the benchmark; fold isolation | large negative: far fewer perturbations and cell contexts, 3 wells per (drug, dose) group | 09 |
| initialisation | FigShare checkpoints exist but are a manual download; authors train from scratch | random init, no released weights loaded | recipe_omission | ticket 09 forbids starting from public weights | large negative at fixed steps; the smoke is 20 steps from random init | 09 |
| paired control | `paired_control_index` pairs every treated cell with a chosen DMSO cell from the same plate | one matched DMSO cell per well via `Controls.cell_for(control_group, replicate)`; no `paired_control_index` exists in the corpus | adapter_extension | the corpus has no paired-control column; the replicate index is the only well-level plate link | none expected: same basal state, one deterministic DMSO cell per treated well | 15 |
| preprocessing | `normalize_total` (median library) + `log1p` + top-200 HVG, HVG fitted on the whole dataset | frozen `log1p_cp10k` on the fixed 5,000-gene axis; no per-run HVG step | recipe_substitution | the benchmark data contract is cp10k log1p; test wells expose no library size so the authors' transform cannot be inverted exactly | small monotone per-cell rescaling; no effect on direction ranking; the cross-view HVG axis is declared in ADR 0028 | 04, 09 |
| training data | treated cells + all paired controls from the full atlas | fold's training-role treated wells only; controls enter as the paired encoder input and the Δ baseline | adapter_extension | fold isolation; valid/test expression never enters training | none expected | 04, 09 |
| step budget / early stopping | fixed 100,000 optimizer steps, no epochs, no early stopping (`train_squidiff.py:115`, `train_util.py:188-211`) | fixed 100,000 steps, no early stopping (20 in smoke) | config_only | author recipe | none | 09, 10 |
| checkpoint selection | last saved checkpoint / EMA weights (`save_interval=1e4`, `model_0.9999.pt`) | author epsilon-MSE on the cold valid drugs at the fixed 10% grid (steps 10k..100k), raw and EMA both evaluated, best kept | config_only | ticket 09 fixes a checkpoint grid for fixed-step models; no test signal is used | small positive or neutral: picks the best generalisation point instead of an arbitrary one | 09, 10 |
| EMA sampling | EMA 0.9999 weights are saved and can be loaded by `sample_squidiff` | EMA 0.9999 kept as a checkpoint candidate at every grid point; the smoke selected the raw state (EMA valid 1.1196 vs raw 1.0039 at step 20) | config_only | ticket 09 keeps the EMA endpoint in the candidate set | neutral at 10^5 steps; the EMA lags badly in a 20-step smoke | 09 |
| drug condition | `rFCFP4(1024)` (Morgan, radius 2, `useFeatures=True`) x `log10(dose+1)`, `use_drug_structure=True` | identical; the shared `common.morgan_fingerprints` (PRnet's helper) reproduces `Drug_dose_encoder` bit-for-bit | config_only | authors' cold-drug path; the CLI default is `use_drug_structure=False`, which feeds the clean target to the encoder and is not a cold-drug model | none: the only author-faithful cold-drug setting | 09 |
| dose transform | continuous dose folded into the fingerprint as `log10(dose_uM + 1)` | identical (doses 0.2 / 1 / 10 uM) | config_only | author function reused verbatim | none; drug identity and dose are collinear by design | 09 |
| prediction space | absolute denoised expression; authors compare group means | absolute `log1p_cp10k` sample minus the matched plate DMSO mean = benchmark Δ | adapter_extension | ticket 04 defines the shared Δ space | none expected; the counterfactual baseline is the plate DMSO mean | 04 |
| sampling | `z_sem` precomputed once, DDIM from pure noise, 1000 steps (`sample_squidiff.py:148-158`) | identical (10 DDIM steps in smoke), one draw per well instead of N draws per group | config_only | smoke budget; the benchmark scores per-well rows and averages 3 wells per group | more sampling variance than the authors' group-level average; no expected bias | 15 |
| architecture / optimizer | MLP 2048x3, time embed 2048, latent 60; AdamW 1e-4, wd 0, linear anneal, batch 64, no grad clip, no dropout | identical via `create_model` / `create_gaussian_diffusion` / `AdamW` | config_only | authors' code reused; only dimensions and steps are configured | none | 09 |
| time | no time input (single 24 h timepoint) | no time input; `time_h` kept in metadata only | config_only | author design | none | 09 |
