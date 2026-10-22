# PRnet adapter — run card (ticket 13)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.1643** (95% CI 0.1545–0.1750), rank 4/12, `ideal_fit: full`. This is the second-best entrant on the primary metric and the best entrant on Spearman (0.384); it is `full`. The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/prnet/` (vendored, commit `f19174bd`, Apache-2.0).
  No file under `pert_model/` is modified or generated (`__pycache__` is
  disabled).
- **Environment:** conda `arc-stage2-bench-prnet` (Python 3.9, torch
  `2.0.1+cu118`, RDKit 2025.09.2), cloned from the setup bringup env and given
  the CUDA wheel (ticket 10).
- **Config:** [`config.json`](config.json) records the author architecture,
  optimizer, schedule, epoch budget, and the 18,164-gene output setting.
- **Code:** `stage2/03_bench_models/temp_env_all_18k/scripts/adapters/prnet/run.py`; shared plumbing in
  `adapters/common.py`; model code is imported through `PYTHONPATH`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, fold 0)
BENCH_SMOKE=1 bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model prnet --fold 0 --timeout 3600

# full recipe (ticket 16 runs folds 0..4 under the 16 h/model budget)
bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model prnet --fold 0
```

Artifacts land in `stage2/03_bench_models/temp_env_all_18k/output/prnet/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Smoke result (fold 0, 1 epoch)

| quantity | value |
|---|---|
| fit wall clock | 4.3 s |
| predict wall clock (valid / test) | 3.5 s / 3.4 s |
| peak RSS | 1.6 GB |
| valid Gaussian NLL | −0.0520 |
| test AUPRC k=200 (headline) | 0.0761 |
| test AUPRC k=100 / k=500 | 0.0562 / 0.1107 |
| test PCC / Spearman | 0.0206 / 0.1936 |
| test R² / RMSE | −13.64 / 0.4266 |
| checkpoint | 4.7 M params, 18,164-gene output |

All 837 test wells produced finite 18,164-vectors and were averaged to the 276
test groups. A one-epoch model is far from the author schedule (100 epochs);
the low AUPRC is expected and is **not** evidence about PRnet.

## How the adapter works

- Each training well is one row: input is the matched plate DMSO mean, the
  drug label is the author's `rFCFP4(1024) x log10(dose_uM + 1)`, and the
  target is the well's `log1p_cp10k` expression.
- The author's `PRnet.get_PGM()` is trained with Gaussian NLL, Adam
  (1e-3, wd 1e-8) and `ReduceLROnPlateau(0.5, patience 10)`; the validation
  loss on the cold valid drugs selects the checkpoint (patience 20).
- Prediction uses the Gaussian mean head with zero noise and subtracts the
  same DMSO mean, i.e. Δ in the shared benchmark space.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 32 | 10.8 | 3.5 | 4.1 | 3456 | 33 | 12 |
| 1 | completed | 35 | 13.7 | 4.2 | 3.2 | 3882 | 61 | 40 |
| 2 | completed | 35 | 14.1 | 4.3 | 3.3 | 3882 | 67 | 46 |
| 3 | completed | 32 | 11.1 | 4.2 | 3.3 | 3882 | 46 | 25 |
| 4 | completed | 35 | 14.1 | 4.3 | 3.3 | 2811 | 67 | 46 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
