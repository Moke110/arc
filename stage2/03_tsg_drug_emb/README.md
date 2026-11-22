# TranSiGen drug representation comparison

Seven representations were compared on all 459 tenv parents across three cold-drug
folds at seed 42. Each arm used FiLM, the same G_5000 input, and the same fold rows.
The first drug-encoder layer accepted the native feature width for each arm and
produced 400 outputs. Later layers shared architecture and initialization. Each
response model trained from scratch. Inner-validation MSE selected checkpoints.

| Representation | Parent-equal OOF Δ PCC | Paired Δ PCC vs ECFP4C | Defined parents |
|---|---:|---:|---:|
| Mean_cell | 0.23164 | -0.03746 | 459/459 |
| ECFP4C | 0.26910 | +0.00000 | 459/459 |
| ECFP6C | 0.26745 | -0.00165 | 459/459 |
| FCFP4C | 0.27892 | +0.00983 | 459/459 |
| FCFP6C | 0.26913 | +0.00003 | 459/459 |
| E3FP | 0.26869 | -0.00041 | 459/459 |
| KPGT | 0.28279 | +0.01370 | 459/459 |
| UniMol | 0.27398 | +0.00488 | 459/459 |

The primary result averages finite per-profile Δ PCC within each held-out parent
and gives each parent equal weight. The paired result uses profiles with defined
scores for both arms before averaging within parent. Fold SD describes dispersion
across the three folds for this single seed.

- [Per-fold metrics](figure/data/runs.csv) and [equal-fold summary](figure/tenv_metrics.png)
- [Native training loss and validation MSE](figure/tenv_learning_curves.png)
- [Test metrics by measured dose](figure/tenv_metrics_cold_drug_dose.png)
- [Primary pooled result](figure/data/pooled.csv) and [paired effects](figure/data/paired.csv)
- [Run inventory](log/run_inventory.json) and [acceptance](log/acceptance.json)

Run `python figure/script/render.py --output <new-directory>` from this step to
reproduce the plots using the checked-in CSV files.

The [Tahoe comparison](log/tahoe/README.md) uses the three frozen Tahoe cold-drug folds.
