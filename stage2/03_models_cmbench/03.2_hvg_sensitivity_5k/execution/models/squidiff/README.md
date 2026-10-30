# Squidiff adapter — 5k sensitivity run card

> Derived from `stage2/03_models_cmbench/03.1_cold_drug_18k/models/squidiff/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.0523** (95% CI 0.0517–0.0530), rank 11/12, `ideal_fit: full`. Paired gap vs `delta_mean`: **-0.1968** (paired CI -0.2150..-0.1790); convergence failure persists: R² -2.98e5 and RMSE 74.3. This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/Squidiff/source/` (vendored, commit `abdfc27d`, BSD-3). No
  file under `pert_model/` is modified or generated (`__pycache__` is
  disabled via `PYTHONDONTWRITEBYTECODE=1` and `sys.dont_write_bytecode`).
- **Environment:** conda `arc-stage2-bench-squidiff` (Python 3.13.16, torch
  `2.14.1+cu130`, numpy 2.5.4, RDKit 2026.09.1), cloned from
  `arc-stage2-bringup-squidiff` and given the CUDA wheel (ticket 10). Lock:
  [`stage2/_envs/squidiff.lock`](../../../../../_envs/squidiff.lock).
- **Config:** [`config.json`](config.json) records the author architecture,
  optimizer, schedule, step budget, and the 5,000-gene setting.
- **Code:** `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/adapters/squidiff/run.py`; shared plumbing
  in `adapters/common.py`; model code is imported through `PYTHONPATH`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (20 optimization steps, 10 DDIM sampling steps, fold 0)
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock "${TMPDIR:-/tmp}/stage2_bench_gpu.lock" \
  bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model squidiff --fold 0 --timeout 1800

# full recipe (ticket 16: folds 0..4, fixed 100,000 steps)
bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model squidiff --fold 0
```

Artifacts land in `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/units/squidiff/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 10.0 s |
| total wall | 21 s |
| peak RSS | 1.6 GB |
| test AUPRC k=200 | 0.0550 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 1880 | 1799 | 37.3 | 39.7 | 1.6 | 100000 / 10000 |
| 1 | completed | 1877 | 1800 | 37.0 | 36.3 | 1.6 | 100000 / 70000 |
| 2 | completed | 1887 | 1810 | 37.0 | 37.0 | 1.6 | 100000 / 70000 |
| 3 | completed | 1884 | 1806 | 36.7 | 37.2 | 1.6 | 100000 / 90000 |
| 4 | completed | 1881 | 1804 | 35.8 | 36.8 | 1.6 | 100000 / 70000 |

Machine-readable copy: `results/run_cards.csv`.
