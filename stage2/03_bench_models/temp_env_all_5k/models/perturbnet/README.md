# PerturbNet adapter — 5k sensitivity run card

> Derived from `stage2/03_bench_models/temp_env_all_18k/models/perturbnet/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.1929** (95% CI 0.1790–0.2083), rank 4/12, `ideal_fit: partial`. Paired gap vs `delta_mean`: **-0.0563** (paired CI -0.0688..-0.0434); was the best entrant at 18k (-0.0128); the HVG axis moves it behind the baseline. This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/perturbnet/` (vendored, commit `e1375b2c`,
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
  three stages and the 5,000-gene output setting.
- **Code:** `stage2/03_bench_models/temp_env_all_5k/scripts/adapters/perturbnet/run.py`; shared
  plumbing in `adapters/common.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md) — 20 rows, including the 5,000-gene
  output, the ZINC → train-drug ChemicalVAE substitution, the omitted
  ZINC corpus, random init, controls missing from the training rows, and
  `ideal_fit: partial`.

## Command

```bash
# one-fold smoke (one epoch per stage). Holds the shared GPU lock.
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock /tmp/opencode/stage2_bench_gpu.lock \
  bash stage2/03_bench_models/temp_env_all_5k/scripts/run.sh --model perturbnet --fold 0 --timeout 1800

# full recipe (ticket 16 runs folds 0..4)
bash stage2/03_bench_models/temp_env_all_5k/scripts/run.sh --model perturbnet --fold 0
```

Artifacts land in `stage2/03_bench_models/temp_env_all_5k/output/perturbnet/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Recipe as run

| stage | data | settings |
|---|---|---|
| ChemicalVAE | 276 unique training-drug representative SMILES (sanitized; author 80/20 split → 220 trained) | 120x35 one-hot, 196-d latent, 525 epochs, batch 128, Adam 1e-4; mu/std recomputed on the 276 structures |
| scVI cell VAE | training-role treated wells, raw counts (2,523 x 5,000) | ZINB, 700 epochs, `n_latent=10`, `n_hidden=128`, `n_layers=2`, `n_batch=1`, frequency 20 |
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

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 16.5 s |
| total wall | 39 s |
| peak RSS | 5.6 GB |
| test AUPRC k=200 | 0.1548 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 467 | 442 | 10.9 | 9.3 | 6.8 | 100 / 56 |
| 1 | completed | 457 | 432 | 11.3 | 9.5 | 6.9 | 100 / 39 |
| 2 | completed | 457 | 432 | 11.4 | 9.4 | 6.9 | 100 / 37 |
| 3 | completed | 465 | 439 | 11.2 | 9.3 | 7.6 | 100 / 58 |
| 4 | completed | 459 | 434 | 11.4 | 9.5 | 7.6 | 100 / 51 |

Machine-readable copy: `results/run_cards.csv`.
