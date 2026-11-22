# TranSiGen drug representation comparison

Seven representations were compared on all 352 Tahoe parents across three cold-drug
folds at seed 42. Each arm used FiLM, the same G_5000 input, and the same fold rows.
The first drug-encoder layer accepted the native feature width for each arm and
produced 400 outputs. Later layers shared architecture and initialization. Each
response model trained from scratch. Inner-validation MSE selected checkpoints.

| Representation | Parent-equal OOF Δ PCC | Paired Δ PCC vs ECFP4C | Defined parents |
|---|---:|---:|---:|
| Mean_cell | 0.46100 | -0.07616 | 352/352 |
| ECFP4C | 0.53715 | +0.00000 | 352/352 |
| ECFP6C | 0.53459 | -0.00256 | 352/352 |
| FCFP4C | 0.52729 | -0.00986 | 352/352 |
| FCFP6C | 0.52900 | -0.00815 | 352/352 |
| E3FP | 0.53453 | -0.00262 | 352/352 |
| KPGT | 0.52187 | -0.01529 | 352/352 |
| UniMol | 0.51491 | -0.02224 | 352/352 |

The primary result averages finite per-profile Δ PCC within each held-out parent
and gives each parent equal weight. The paired result uses profiles with defined
scores for both arms before averaging within parent. Fold SD describes dispersion
across the three folds for this single seed.

- [Per-fold metrics](../../figure/data/tahoe/runs.csv) and [equal-fold summary](../../figure/tahoe_metrics.png)
- [Native training loss and validation MSE](../../figure/tahoe_learning_curves.png)
- [Test metrics by measured dose](../../figure/tahoe_metrics_cold_drug_dose.png)
- [Primary pooled result](../../figure/data/tahoe/pooled.csv) and [paired effects](../../figure/data/tahoe/paired.csv)
- [Run inventory](run_inventory.json) and [acceptance](acceptance.json)

Run `python figure/script/render.py --source tahoe_bulk --prefixed --output <new-directory>` from this step to
reproduce the plots using the checked-in CSV files.
