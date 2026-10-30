# Squidiff adapter — run card (ticket 15)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.0212** (95% CI 0.0207–0.0218), rank 11/12, `ideal_fit: full`. **Convergence failure at 18,164 genes:** cold-valid ε-MSE stayed at 0.976 (random = 1.0), so the samples are out of distribution. The row is retained as a documented failure, not a meaningful ranking. The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/Squidiff/source/` (vendored, commit `abdfc27d`, BSD-3). No
  file under `pert_model/` is modified or generated (`__pycache__` is
  disabled via `PYTHONDONTWRITEBYTECODE=1` and `sys.dont_write_bytecode`).
- **Environment:** conda `arc-stage2-bench-squidiff` (Python 3.13.16, torch
  `2.14.1+cu130`, numpy 2.5.4, RDKit 2026.09.1), cloned from
  `arc-stage2-bringup-squidiff` and given the CUDA wheel (ticket 10). Lock:
  [`stage2/_envs/squidiff.lock`](../../../../../_envs/squidiff.lock).
- **Config:** [`config.json`](config.json) records the author architecture,
  optimizer, schedule, step budget, and the 18,164-gene setting.
- **Code:** `stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/squidiff/run.py`; shared plumbing
  in `adapters/common.py`; model code is imported through `PYTHONPATH`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (20 optimization steps, 10 DDIM sampling steps, fold 0)
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock "${TMPDIR:-/tmp}/stage2_bench_gpu.lock" \
  bash stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/run.sh --model squidiff --fold 0 --timeout 1800

# full recipe (ticket 16: folds 0..4, fixed 100,000 steps, ~57 min/fold)
bash stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/run.sh --model squidiff --fold 0
```

## Full-run projection

A 50-step timing probe of the same loop at batch 64 / 18,164 genes measured
**34.3 ms/step** (peak 3.8 GB GPU), i.e. ≈ **57 min per fold** and ≈ **4.8 h
for the five folds** — inside the ticket-10 budget, a little above its
43 min/fold estimate because each step copies the batch from host memory.

Artifacts land in `stage2/03_models_cmbench/03.1_cold_drug_18k/output/units/squidiff/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Smoke result (fold 0, 20 steps / 10 sampling steps)

| quantity | value |
|---|---|
| fit wall clock | 15.5 s (9.4 s training loop + 10 grid evaluations + 1 checkpoint write) |
| predict wall clock (valid / test) | 4.3 s / 4.2 s |
| peak GPU / RSS | 4.6 GB / 2.1 GB |
| parameters | 164,609,840 at 18,164 genes (author `num_layers=3`, width 2048) |
| checkpoint | `checkpoint/model.pt`, 658,490,081 B (658 MB, raw weights @ step 20) |
| valid epsilon-MSE (author loss) | raw 1.00386, EMA 1.11960 |
| test AUPRC k=200 (headline) | 0.0224 |
| test AUPRC k=100 / k=500 | 0.0152 / 0.0408 |
| test PCC / Spearman | −0.00002 / 0.0033 |
| test R² / RMSE | −77,215 / 31.68 |
| groups | 837 test wells → 276 group means (92 drugs × 3 doses) |

All 837 test wells produced finite 18,164-vectors and were averaged to the 276
test groups. The valid loss is the author epsilon-MSE and still decreasing
(1.068 at step 2 → 1.004 at step 20); with 20 steps the denoiser predicts ≈ 0,
so the DDIM chain from pure noise does **not** reach expression scale and the
smoke predictions are noise-scale (mean 24, max 347) — hence the very negative
R². This is the expected behaviour of an essentially untrained diffusion model,
not a scoring bug; it is **not** evidence about Squidiff.

## How the adapter works

- Each training well is one row. The encoder input is
  `[matched DMSO cell ; rFCFP4(1024) × log10(dose_uM + 1)]` — the authors'
  `EncoderMLPModel` replaces `x_start` with exactly that concatenation
  (`MLPModel.py:200-201`); the drug vector is built with the shared PRnet
  fingerprint helper and is byte-for-byte the authors' `Drug_dose_encoder`.
- **Control pairing rule:** every well takes the DMSO cell on its own plate
  selected by the well's 0-based `replicate` index (`Controls.cell_for`), i.e.
  one matched control cell per treated well. The corpus has no
  `paired_control_index`, so this is the reproducible equivalent; controls are
  the 16 global DMSO groups and are never held out.
- Training uses the authors' `GaussianDiffusion.training_losses` (epsilon MSE,
  uniform timesteps), `AdamW(lr 1e-4, wd 0)`, the linear anneal to 0, batch 64,
  EMA 0.9999 and no gradient clipping. The step budget is fixed (100,000; 20 in
  smoke) with **no early stopping**.
- At each 10 % grid step the author loss is evaluated on the fold's cold valid
  drugs for both the raw and the EMA parameters; the best state is kept
  (`checkpoint/model.pt`). Test expression never enters training or selection.
- Prediction precomputes the semantic latent once per batch (as
  `sample_squidiff.py` does), runs DDIM from seeded pure noise at the
  checkpoint's sampling budget, and returns `generated − matched plate DMSO
  mean` in the shared `log1p_cp10k` Δ space.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 3394 | 3294.9 | 44.5 | 42.4 | 7627 | 100000 | 80000 |
| 1 | completed | 3390 | 3291.9 | 44.1 | 41.9 | 7627 | 100000 | 80000 |
| 2 | completed | 3393 | 3294.6 | 42.5 | 43.6 | 7614 | 100000 | 80000 |
| 3 | completed | 3398 | 3300.9 | 42.5 | 44.0 | 7613 | 100000 | 100000 |
| 4 | completed | 3390 | 3292.2 | 42.5 | 44.0 | 7368 | 100000 | 70000 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
