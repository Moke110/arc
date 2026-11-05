# CMonge adapter — 5k sensitivity run card

> Derived from `stage2/03_models_cmbench/03.1_cold_drug_18k/models/cmonge/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.1153** (95% CI 0.1108–0.1200), rank 10/12, `ideal_fit: full`. Paired gap vs `delta_mean`: **-0.1339** (paired CI -0.1538..-0.1150); small worsening of -0.0082 vs 18k. This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/CMonge/source/` (vendored, commit `195f5cd1`, MIT). No file
  under `pert_model/` is modified or generated (`PYTHONDONTWRITEBYTECODE=1` and
  `sys.dont_write_bytecode = True`).
- **Environment:** `stage2/_envs/cmonge` — a `cp -a` copy of the setup
  bringup venv (the model-local `.venv`, removed in ticket 04). Python 3.11.15, jax/jaxlib
  `0.4.38`, flax 0.10.4, optax 0.2.5, ott-jax 0.5.2, rdkit 2023.09.6,
  scanpy 1.11.5). JAX backend: **CPU** (`jax.default_backend() == "cpu"`); the
  CUDA 12 jaxlib install did not complete inside the 1 h environment budget, so
  the row is recorded as a `reproducibility_break` in
  [`FIDELITY.md`](FIDELITY.md). The 18k probe measured 32 min for a 50-step AE + 10,000-iteration fit on
  fold 0 (CPU, shared machine); the 5k axis is re-timed by the fold-0 smoke. Lock:
  `envs/cmonge.lock`.
- **Config:** [`config.json`](config.json) — AE `data_dim=5000`, latent 50,
  2x512; conditional Monge map 4x64, 50-d OT latent, **10,000 iterations** (the
  `conditional-monge-sciplex` setting, not the 500,000 all-drug setting);
  Sinkhorn `epsilon_fitting=1` + Monge gap `epsilon_regularizer=1e-2`; AdamW
  1e-4 / wd 1e-5; batch 512; 200 RDKit props + log-dose.
- **Code:** `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/adapters/cmonge/run.py`; shared plumbing in
  `adapters/common.py`; model code is imported from `pert_model/CMonge/source`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (fold 0, shared GPU/queue lock)
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock "${TMPDIR:-/tmp}/stage2_bench_gpu.lock" \
  bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model cmonge --fold 0 --timeout 1800

# full recipe (ticket 16 runs folds 0..4)
bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model cmonge --fold 0 --timeout 21600
```

Artifacts land in `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/units/cmonge/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 71.2 s |
| total wall | 89 s |
| peak RSS | 1.4 GB |
| test AUPRC k=200 | 0.1222 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 2709 | 2693 | 5.7 | 5.7 | 1.2 | 10000 / 3000 |
| 1 | completed | 5740 | 5724 | 5.7 | 6.3 | 1.5 | 10000 / 1000 |
| 2 | completed | 6773 | 6757 | 5.8 | 5.6 | 1.5 | 10000 / 1000 |
| 3 | completed | 5946 | 5931 | 5.7 | 5.6 | 1.5 | 10000 / 2000 |
| 4 | completed | 971 | 954 | 5.8 | 5.9 | 1.1 | 10000 / 1000 |

Machine-readable copy: `results/run_cards.csv`.
