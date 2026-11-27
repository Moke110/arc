# Step 05.4 · TranSiGen Delta PCC sweep

Step 05.4 tested whether adding a bounded Pearson correlation penalty to
TranSiGen improved held-out treatment response. It used the Step 05.3
objective with unit MSE weights for control reconstruction, treated
reconstruction and Δ prediction, and a Δ L1 weight of 1.75. For each
condition, it computed Pearson correlation across genes between predicted
and measured Δ using the measured DMSO control from the same plate. The
added loss was the batch mean of `(1 - clamp(PCC, -1, 1)) / 2` multiplied
by `alpha`. The sweep tested 0, 0.1, 0.2, 0.3, 0.4 and 0.5 on Denoise-DEG
5k with frozen KPGT vectors, tenv cold-drug folds 0–2 and seed 42. Held-out
scoring used frozen HVG-DEG 5k bandwidths.
[Step 05.3 report](../05.3_mse-l1_mse-x/README.md) ·
[ADR 0092](../../../docs/adr/0092-transigen-bounded-delta-similarity-sweeps.md) ·
[Denoise-DEG 5k manifest](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) ·
[Frozen tenv folds](../../03_data_train/output/frozen-folds/tenv/manifest.json) ·
[Bandwidths](log/bandwidths.csv)

## Accepted results and follow-up

The two-epoch smoke test at weights 0 and 0.5 on fold 0 passed exact
prediction replay and byte-identical figure reproduction from the report
CSVs. The formal sweep completed 18 fresh CUDA fits with exact independent
replay, 162 held-out fold scores and 54 equal-fold aggregates. The three
fresh zero-weight fits exactly matched the Step 05.3 unit-reconstruction
fits on objective traces, selected and stopping epochs, validation and test
predictions, Δ arrays, fold scores and eligibility counts.
[Smoke receipt](log/338-smoke.json) · [Formal receipt](log/339-formal.json) ·
[Run ledger](log/runs.csv) · [Fold scores](log/fold-scores.csv)

Weight 0 had the highest mean Δ PCC (0.431920) and the lowest mean Δ RMSE
(0.249215). Weight 0.3 was the closest positive weight on both metrics,
with mean Δ PCC of 0.431519 and mean Δ RMSE of 0.249290. Step 05 retained
the zero-weight objective. Step 05.5 tested a cosine penalty with the same
base objective and its own fresh zero-weight fits.
[Equal-fold scores](figure/data/delta-metrics.csv) ·
[Step 05 selected objective](../README.md#selected-objective) ·
[Step 05.5 cosine sweep](../05.5_mse-l1_cos/README.md)

The current comparison figures include wJAC-2000. Paired differences
subtract the fresh zero-weight score for the same fold and seed. The
learning curves show training and validation through each fit's stopping
epoch.
[Metric scores](figure/metrics.png) ·
[Paired differences](figure/metrics-paired.png) ·
[Learning curves](figure/learning-curve.png) ·
[Plot CSVs](figure/data/) · [Render receipt](log/figure-render.json) ·
[Figure verification](AGENTS.md#verify-the-current-figures)

## Reproduce the smoke

The study guide covers the bounded two-fit run, scoring and CSV-only
smoke renders. [Smoke procedure](AGENTS.md#reproduce-the-smoke)

## Reproduce the formal sweep

The acceptance receipt pins the joint matrix and PCC report manifests.
[Sealed PCC replay](AGENTS.md#replay-a-sealed-fit) ·
[Joint formal procedure](../README.md#reproduce-the-original-formal-delta-similarity-sweeps) ·
[Registration verification record](../../util/training/verification.md)
