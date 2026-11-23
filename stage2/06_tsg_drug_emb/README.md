# TranSiGen drug representation comparison

See the [historical evidence index](../03_view/03.1_hvg-deg_views/historical/README.md) for sealed records.

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

## Raw-first inputs, evaluation and figures

New comparisons use the [raw-first dose input pair](../07_tsg_dose/README.md#raw-first-inputs-evaluation-and-figures)
prepared from the default [HVG-DEG full view](../03_view/03.1_hvg-deg_views/README.md).
The preparation command freezes the selected source from the default full
view. Use `--inputs` to reuse a verified source freeze.
A comparison receipt pins the source's pair, expression payload, view, corpus,
ordered G_5000 axis, training seed and shared representation manifest. The
representation reader aligns unchanged vectors at native width to the verified
drug-parent identities. Each source retains its original cold-drug fold protocol.

```sh
D=stage2/06_tsg_drug_emb
python "$D/execution/raw_first.py" prepare --source tenv \
  --receipt /path/to/new/comparison.json
python "$D/execution/raw_first.py" evaluate --inputs /path/to/new/comparison.inputs/tenv \
  --receipt /path/to/new/comparison.json --arm ECFP4C --fold 0 \
  --prediction /path/to/prediction.npy --prediction-ids /path/to/profile_ids.csv \
  --output /path/to/new/evaluation
python "$D/execution/raw_first.py" render --summary /path/to/new/evaluation/summary.csv \
  --output /path/to/new/metrics.png
```

Predictions contain absolute expression values in fold order, with one
`profile_id` CSV row per prediction and columns in the published axis order.
Evaluation uses the new matched-control matrix and preserves the study's
parent-equal Δ Pearson summary. Rendering reads only that recorded CSV.
For a separate materialized checkout, pass its accepted representation
directory using `--representations` on both input commands.

The reported historical results above remain tied to their original input
and feature receipts. See the [accepted study records](log/).
