# CPA adapter — 5k sensitivity run card

> Derived from `stage2/03_models_cmbench/03.1_cold_drug_18k/models/cpa/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.1559** (95% CI 0.1513–0.1605), rank 5/12, `ideal_fit: full`. Paired gap vs `delta_mean`: **-0.0932** (paired CI -0.1120..-0.0746); gap improved by +0.0147 vs 18k but remains below the baseline. This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/CPA/source/` (vendored, commit `fbd7c025`, BSD-3-Clause).
  No file under `pert_model/` is modified or generated.
- **Environment:** conda `arc-stage2-bench-cpa` (Python 3.9, torch
  `2.0.1+cu118`, scvi-tools 0.20.3, anndata 0.9.2), cloned from the setup
  bringup env and given the CUDA wheel (ticket 10).
- **Config:** [`config.json`](config.json) records the author architecture and
  the three-optimizer training plan.
- **Code:** `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/adapters/cpa/run.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, fold 0)
BENCH_SMOKE=1 bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model cpa --fold 0 --timeout 7200

# full recipe (ticket 16)
bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model cpa --fold 0
```

Artifacts land in `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/units/cpa/fold_0/`.

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 16.7 s |
| total wall | 32 s |
| peak RSS | 1.6 GB |
| test AUPRC k=200 | 0.1228 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 1920 | 1906 | 5.8 | 4.3 | 1.5 | 340 / 340 |
| 1 | completed | 2599 | 2583 | 7.2 | 4.5 | 1.6 | 460 / 400 |
| 2 | completed | 2090 | 2078 | 4.3 | 4.2 | 1.6 | 370 / 340 |
| 3 | completed | 2359 | 2343 | 4.4 | 7.0 | 1.6 | 420 / 410 |
| 4 | completed | 1577 | 1561 | 7.2 | 4.4 | 1.6 | 280 / 270 |

Machine-readable copy: `results/run_cards.csv`.
