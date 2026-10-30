# PerturbNet adapter — run card (ticket 14)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.1861** (95% CI 0.1699–0.2031), rank 3/12, `ideal_fit: partial`. This is the best entrant on the primary metric, but it remains below the `delta_mean` baseline; it is `partial` because it takes no control transcriptome. The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/PerturbNet/source/` (vendored, commit `e1375b2c`,
  GPL-3.0). No file under `pert_model/` is modified or generated
  (`PYTHONDONTWRITEBYTECODE=1` + `sys.dont_write_bytecode`).
- **Environment:** conda `arc-stage2-bench-perturbnet`, the authors' faithful
  pin: **Python 3.7.16**, `torch 1.13.1+cu117`, `scvi-tools 0.7.1`,
  `numpy 1.19.5`, `anndata 0.7.8`, `pandas 1.3.5`, `scipy 1.7.3`,
  `rdkit-pypi 2022.9.5`, `scanpy 1.8.2`.
  - **TensorFlow 1.15 caveat:** the pin is only needed by the
    normalised-expression cell VAE (`data_vae/vae.py`), which this corpus does
    not use — counts take the scVI/PyTorch branch (`cinn/flow.py:3289-3356`).
    TF 1.15 is **not installed**: its CUDA 10 runtime segfaults beside
    torch 1.13/cu117 on this host. `perturbnet.util` imports that module
    transitively, so the adapter installs a two-line `tensorflow` stub before
    importing the vendored package. No code under `pert_model/` is changed.
  - The ticket-10 py3.10 fallback was **not** used; no
    `reproducibility_break` is recorded.
  - `pip freeze` lock: [`../../envs/perturbnet.lock`](../../envs/perturbnet.lock)
    (sha256 `b559d36a251dc9059012c4e5971c1c682245d64e73269dc9274fd157b0df0920`).
- **Config:** [`config.json`](config.json) records every author default for the
  three stages and the 18,164-gene output setting.
- **Code:** `stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/perturbnet/run.py`; shared
  plumbing in `adapters/common.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md) — 20 rows, including the 18,164
  extrapolation, the ZINC → train-drug ChemicalVAE substitution, the omitted
  ZINC corpus, random init, controls missing from the training rows, and
  `ideal_fit: partial`.

## Command

```bash
# one-fold smoke (one epoch per stage). Holds the shared GPU lock.
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock "${TMPDIR:-/tmp}/stage2_bench_gpu.lock" \
  bash stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/run.sh --model perturbnet --fold 0 --timeout 1800

# full recipe (ticket 16 runs folds 0..4)
bash stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/run.sh --model perturbnet --fold 0
```

Artifacts land in `stage2/03_models_cmbench/03.1_cold_drug_18k/output/units/perturbnet/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Recipe as run

| stage | data | settings |
|---|---|---|
| ChemicalVAE | 276 unique training-drug representative SMILES (sanitized; author 80/20 split → 220 trained) | 120x35 one-hot, 196-d latent, 525 epochs, batch 128, Adam 1e-4; mu/std recomputed on the 276 structures |
| scVI cell VAE | training-role treated wells, raw counts (2,523 x 18,164) | ZINB, 700 epochs, `n_latent=10`, `n_hidden=128`, `n_layers=2`, `n_batch=1`, frequency 20 |
| cINN | every training well, condition = 196-d chem latent + 3-bin dose one-hot | `ConditionalFlatCouplingFlow` (20 flows, 10-d latent, hidden 1024 x 2, BN conditioner, `activation=none`) wrapped by `Net2NetFlow_scVIChemStdStatesFlow`; 100 epochs, batch 128, Adam lr 4.5e-6, betas (0.5, 0.9), amsgrad; best epoch by author NLL on the cold valid drugs |

No validation/test expression enters any stage; the valid bundle is used only
for the cINN NLL checkpoint selection (scVI encodes it, never trains on it).

## Prediction

`SCVIZ_CheckNet2Net` semantics: the cINN samples a 10-d latent for each
(drug, dose) condition; the scVI decoder maps it to the ZINB rate
(`px_rate`), which is the per-draw predictive mean. `flow_draws = 8` draws are
averaged per well, counts are expressed as `log1p_cp10k` at the median training
library size, and the matched plate DMSO mean (`--controls`) is subtracted.
Dose is the 3-bin one-hot frozen by ticket 09. Drugs whose SMILES cannot be
encoded by the 35-symbol / 120-char ChemicalVAE (2 of 92 test drugs in fold 0)
get the mean (zero) standardized chemical latent.

## Smoke result (fold 0, one epoch per stage)

| quantity | value |
|---|---|
| fit wall clock | 20.7 s |
| predict wall clock (valid / test) | 16.0 s / 9.7 s |
| unit wall clock | 63.0 s |
| peak RSS | 6.3 GB |
| flow draws | 2 (smoke) |
| test AUPRC k=200 (headline) | 0.0968 |
| test AUPRC k=100 / k=500 | 0.0692 / 0.1453 |
| test PCC / Spearman | 0.1570 / 0.2134 |
| test R² / RMSE | −6.30 / 0.3085 |
| valid AUPRC k=200 | 0.0925 |
| checkpoint | 742.8 MB (`cinn.pt` 679 MB, 169.8 M-param flow; `chemicalvae.pt` 16 MB; `scvi_vae.pt` 9 MB) |

Two smokes from the same code and seed (the second forced) produced
bit-identical `pred_valid.npy` / `pred_test.npy` (max diff 0.0), so the seed
and the deterministic flags hold.

All 837 test wells produced finite 18,164-vectors; the runner aggregated the
276 test groups (92 drugs × 3 doses) and wrote `pred_test.npy` `(837, 18164)`.
A one-epoch model is far from the author schedule (525/700/100 epochs); the
low AUPRC is expected and is **not** evidence about PerturbNet.

### Full-recipe scratch check (not a benchmark unit)

Fold 0 run once outside the runner with the real budgets: ChemicalVAE 525
epochs in 32 s, scVI 700 epochs in 184 s, cINN 100 epochs in 295 s (best epoch
54 by cold-valid NLL), 8 flow draws. Final checkpoint 742.8 MB; test AUPRC
k=200 **0.1136** (k=100 0.0955 / k=500 0.1424), PCC 0.089, Spearman 0.084,
R² −1.50, RMSE 0.195. That is ~9 min of fit per fold, so the five-fold budget
is ~1 h — well inside the 16 h/model cap. The number is for wiring/behaviour
only and is not a leaderboard entry.

## Deviations (see `FIDELITY.md` for all 20)

- 18,164 genes vs the authors' 5,087/978 (`recipe_substitution`).
- ChemicalVAE trained from scratch on the fold's 276 training structures
  instead of ZINC, so the ZINC corpus and its 50k-structure z-score are omitted
  (`recipe_substitution` + `recipe_omission`).
- No released checkpoint anywhere (`recipe_omission`, random init).
- 3-bin one-hot dose (`conditioning_dim=199`) instead of raw continuous dose
  (`config_only`, ticket 09).
- Deterministic prediction: mean of the ZINB rate over 8 flow draws instead of
  sampled counts (`config_only`).
- No control transcriptome enters the model; `ideal_fit: partial`.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 543 | 506.0 | 12.8 | 9.8 | 8955 | 100 | 54 |
| 1 | completed | 515 | 479.1 | 12.2 | 10.3 | 9191 | 100 | 64 |
| 2 | completed | 556 | 518.9 | 13.1 | 10.0 | 9457 | 100 | 43 |
| 3 | completed | 546 | 508.8 | 13.1 | 10.4 | 9133 | 100 | 65 |
| 4 | completed | 530 | 494.8 | 12.1 | 10.3 | 9149 | 100 | 55 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
