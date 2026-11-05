# chemCPA adapter — 5k sensitivity run card

> Derived from `stage2/03_models_cmbench/03.1_cold_drug_18k/models/chemcpa/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.1165** (95% CI 0.1129–0.1203), rank 9/12, `ideal_fit: partial`. Paired gap vs `delta_mean`: **-0.1327** (paired CI -0.1498..-0.1152); small improvement of +0.0054 vs 18k. This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/chemCPA/source/` (vendored, MIT). No file under
  `pert_model/` is modified or generated.
- **Environment:** conda `arc-stage2-bench-chemcpa` (Python 3.9, torch
  `2.8.0+cu128`, Lightning 2.2.5, RDKit 2025.09.2), cloned from the setup
  bringup env (already CUDA; no rebuild needed).
- **Config:** [`config.json`](config.json) records the paper's hyperparameters,
  the three-optimizer schedule, and the 5,000-gene from-scratch setting.
- **Code:** `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/adapters/chemcpa/run.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, fold 0)
BENCH_SMOKE=1 bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model chemcpa --fold 0 --timeout 5400

# full recipe (ticket 16)
bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model chemcpa --fold 0
```

Artifacts land in `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/units/chemcpa/fold_0/`.

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 9.4 s |
| total wall | 18 s |
| peak RSS | 1.7 GB |
| test AUPRC k=200 | 0.1493 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 50 | 42 | 2.6 | 2.2 | 1.7 | 201 / 201 |
| 1 | completed | 49 | 40 | 2.3 | 2.4 | 1.8 | 201 / 201 |
| 2 | completed | 52 | 42 | 2.3 | 2.2 | 1.8 | 201 / 199 |
| 3 | completed | 47 | 39 | 2.3 | 2.2 | 1.8 | 201 / 200 |
| 4 | completed | 52 | 40 | 5.0 | 2.3 | 1.8 | 201 / 201 |

Machine-readable copy: `results/run_cards.csv`.
