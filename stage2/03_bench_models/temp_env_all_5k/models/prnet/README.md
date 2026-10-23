# PRnet adapter — 5k sensitivity run card

> Derived from `stage2/03_bench_models/temp_env_all_18k/models/prnet/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.2735** (95% CI 0.2590–0.2891), rank 1/12, `ideal_fit: full`. Paired gap vs `delta_mean`: **+0.0243** (paired drug CI +0.0113..+0.0367); the only entrant above `delta_mean`; it was -0.0345 at 18k. This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/prnet/` (vendored, commit `f19174bd`, Apache-2.0).
  No file under `pert_model/` is modified or generated (`__pycache__` is
  disabled).
- **Environment:** conda `arc-stage2-bench-prnet` (Python 3.9, torch
  `2.0.1+cu118`, RDKit 2025.09.2), cloned from the setup bringup env and given
  the CUDA wheel (ticket 10).
- **Config:** [`config.json`](config.json) records the author architecture,
  optimizer, schedule, epoch budget, and the 5,000-gene output setting.
- **Code:** `stage2/03_bench_models/temp_env_all_5k/scripts/adapters/prnet/run.py`; shared plumbing in
  `adapters/common.py`; model code is imported through `PYTHONPATH`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, fold 0)
BENCH_SMOKE=1 bash stage2/03_bench_models/temp_env_all_5k/scripts/run.sh --model prnet --fold 0 --timeout 3600

# full recipe (ticket 16 runs folds 0..4 under the 16 h/model budget)
bash stage2/03_bench_models/temp_env_all_5k/scripts/run.sh --model prnet --fold 0
```

Artifacts land in `stage2/03_bench_models/temp_env_all_5k/output/prnet/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 5.0 s |
| total wall | 16 s |
| peak RSS | 1.1 GB |
| test AUPRC k=200 | 0.1236 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 20 | 7 | 3.7 | 5.4 | 1.1 | 54 / 33 |
| 1 | completed | 14 | 5 | 2.8 | 2.8 | 1.1 | 56 / 35 |
| 2 | completed | 18 | 5 | 3.0 | 5.5 | 1.1 | 63 / 42 |
| 3 | completed | 15 | 5 | 2.8 | 2.8 | 1.1 | 57 / 36 |
| 4 | completed | 18 | 5 | 2.9 | 3.0 | 1.1 | 55 / 34 |

Machine-readable copy: `results/run_cards.csv`.
