# Step 05.5 · TranSiGen Delta cosine sweep

Step 05.5 tested whether adding a bounded cosine penalty improved held-out
treatment response. It used the Step 05.3 objective with unit MSE weights
for control reconstruction, treated reconstruction and Δ prediction, and
a Δ L1 weight of 1.75. The penalty was the batch mean of
`(1 - clamp(cosine, -1, 1)) / 2` multiplied by `alpha`. Cosine was
computed across the ordered genes on uncentered predicted and measured Δ
from the same measured plate DMSO control. Normalization used epsilon
`1e-8` and gave zero similarity when either vector had zero norm.
[Step 05.3 report](../05.3_mse-l1_mse-x/README.md) ·
[Step 05.4 PCC sweep](../05.4_mse-l1_pcc/README.md) ·
[ADR 0092](../../../docs/adr/0092-transigen-bounded-delta-similarity-sweeps.md)

The sweep tested weights 0, 0.1, 0.2, 0.3, 0.4 and 0.5 with frozen KPGT
vectors on tenv Denoise-DEG 5k, cold-drug folds 0–2 and seed 42. Held-out
scoring used the frozen HVG-DEG 5k bandwidths.
[Denoise-DEG 5k manifest](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) ·
[Frozen tenv folds](../../03_data_train/output/frozen-folds/tenv/manifest.json) ·
[Bandwidths](log/bandwidths.csv)

## Accepted results and follow-up

The two-epoch smoke test at weights 0 and 0.5 on fold 0 passed exact
prediction replay and byte-identical figure reproduction from scoring
CSVs. The formal sweep completed 18 fresh CUDA fits with exact independent
replay, 162 held-out fold scores and 54 equal-fold aggregates. The three
fresh zero-weight fits exactly matched the Step 05.3 unit-reconstruction
fits on objective traces, selected and stopping epochs, validation and test
predictions, Δ arrays, fold scores and eligibility counts.
[Smoke receipt](log/338-smoke.json) · [Formal receipt](log/339-formal.json) ·
[Run ledger](log/runs.csv) · [Fold scores](log/fold-scores.csv)

Weight 0 gave the highest mean Δ PCC (0.431920) and lowest mean Δ RMSE
(0.249215). Weight 0.3 raised the mean Jaccard overlap among the 100
most downregulated genes (`JAC_neg100`) from 0.183801 to 0.185508, with
mean Δ PCC of 0.428635 and mean Δ RMSE of 0.250528. Step 05 retained the
zero-weight objective. Step 06.1 used the same objective in the architecture
comparison.
[Equal-fold scores](figure/data/delta-metrics.csv) ·
[Step 05 selected objective](../README.md#selected-objective) ·
[Step 06.1 architecture comparison](../../06_tsg_shape/06.1_round1/README.md)

The current comparison figures include wJAC-2000. Paired differences
subtract the fresh zero-weight score from the same fold and seed. The
learning curves show training and validation through each fit's stopping epoch.
[Metric scores](figure/metrics.png) ·
[Paired differences](figure/metrics-paired.png) ·
[Learning curves](figure/learning-curve.png) ·
[Plot CSVs](figure/data/) · [Render receipt](log/figure-render.json) ·
[Figure verification](AGENTS.md#verify-the-current-figures)

## Reproduce the smoke

The study guide covers the bounded two-fit run, scoring and CSV-only
smoke renders. [Smoke procedure](AGENTS.md#reproduce-the-smoke)

## Reproduce the formal sweep

The formal receipt pins the joint matrix, preflight and cosine report
manifests. [Sealed cosine replay](AGENTS.md#replay-a-sealed-fit) ·
[Joint formal procedure](../README.md#reproduce-the-original-formal-delta-similarity-sweeps) ·
[Registration verification record](../../util/training/verification.md)
