# CPA adapter — run card (ticket 13)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.0910** (95% CI 0.0863–0.0960), rank 7/12, `ideal_fit: full`. This is the best entrant on R² (0.038) and RMSE (0.140), but 7th on the primary AUPRC; it is `full`. The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/CPA/source/` (vendored, commit `fbd7c025`, BSD-3-Clause).
  No file under `pert_model/` is modified or generated.
- **Environment:** conda `arc-stage2-bench-cpa` (Python 3.9, torch
  `2.0.1+cu118`, scvi-tools 0.20.3, anndata 0.9.2), cloned from the setup
  bringup env and given the CUDA wheel (ticket 10).
- **Config:** [`config.json`](config.json) records the author architecture and
  the three-optimizer training plan.
- **Code:** `stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/cpa/run.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, fold 0)
BENCH_SMOKE=1 bash stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/run.sh --model cpa --fold 0 --timeout 7200

# full recipe (ticket 16)
bash stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/run.sh --model cpa --fold 0
```

Artifacts land in `stage2/03_models_cmbench/03.1_cold_drug_18k/output/units/cpa/fold_0/`.

## Smoke result (fold 0, 1 epoch)

| quantity | value |
|---|---|
| fit wall clock | 14.0 s |
| predict wall clock (valid / test) | 5.5 s / 5.3 s |
| peak RSS | 2.1 GB |
| valid reconstruction loss (`val_recon`) | 1.85e4 |
| test AUPRC k=200 (headline) | 0.0759 |
| test AUPRC k=100 / k=500 | 0.0561 / 0.1096 |
| test PCC / Spearman | 0.0194 / 0.1488 |
| test R² / RMSE | −15.14 / 0.4497 |
| checkpoint | 30.3 M params, 18,164-gene output |

All 837 test wells produced finite vectors and were averaged to 276 groups.
The one-epoch reconstruction loss is still far above the author's converged
value; treat the score as a plumbing result only.

## Cold-drug design

`setup_anndata` runs on the full fold AnnData — train, valid, **and** the 837
test metadata rows (zero expression, `split='test'`) — so the 92 cold test
drugs get a `pert_encoder` row and a frozen Morgan-2048 fingerprint before
training. Only `split=='train'` rows (2,523 treated + 272 controls) are
optimised. The cold drug's per-drug `beta`/`bias` dose-scaler row is never
trained and stays at its init (`β=1`, `bias=0`), as declared in the ledger.

Prediction builds a metadata AnnData whose `X` is the matched plate DMSO mean
and whose condition is the test well's drug/dose; `CPA.predict` returns the
absolute expression and the adapter subtracts that same DMSO mean.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 1956 | 1932.0 | 8.1 | 4.7 | 3962 | 340 | 320 |
| 1 | completed | 3175 | 3151.2 | 5.1 | 4.9 | 3962 | 550 | 450 |
| 2 | completed | 1056 | 1032.3 | 5.1 | 6.4 | 3958 | 180 | 180 |
| 3 | completed | 2901 | 2878.0 | 5.0 | 4.8 | 3958 | 500 | 480 |
| 4 | completed | 2503 | 2478.9 | 7.7 | 4.7 | 3958 | 430 | 420 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
