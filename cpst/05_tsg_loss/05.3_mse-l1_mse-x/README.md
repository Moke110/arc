# Step 05.3 · TranSiGen reconstruction MSE sweep

Step 05.3 tested how strongly TranSiGen should reconstruct control and treated
expression while learning the treatment response. It fixed the Delta L1
weight selected in Step 05.2 at 1.75 and shared one reconstruction MSE
weight, `alpha`, across both expression terms.
[Step 05.2 report](../05.2_mse_l1/README.md) ·
[Spec 11](https://github.com/Moke110/carps/issues/334) ·
[ADR 0091](../../../docs/adr/0091-transigen-reconstruction-mse-sweep.md)

The objective was `0.1 * K + MSE_Δ + 1.75 * L1_Δ + alpha * (MSE_X1 + MSE_X2)`,
where `K` sums the three native KL terms. The sweep used weights 0, 0.5,
0.75, 0.9, 1, 1.1, 1.25, 1.5 and 2 on Denoise-DEG 5k with frozen KPGT
vectors, tenv cold-drug folds 0–2 and seed 42. Measured same-plate DMSO
served as the model input and the reference for predicted and observed Δ.
Each fit selected the earliest minimum of its registered validation objective
with 30-epoch patience and a 1,000-epoch cap. Held-out scoring used the frozen
HVG-DEG 5k bandwidths.
[Denoise-DEG 5k manifest](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) ·
[Frozen tenv folds](../../03_data_train/output/frozen-folds/tenv/manifest.json) ·
[Bandwidths](log/bandwidths.csv)

## Accepted results and follow-up

A two-epoch comparison at weights 0 and 1 on fold 0 passed smoke verification.
The formal sweep completed 27 fresh fits with exact independent replay,
243 held-out fold scores and 81 equal-fold aggregates. All three unit-weight
fits preserved the Step 05.2 objective arithmetic and exactly reproduced its
traces, selected and stopping epochs, validation and test predictions, Delta
arrays, fold scores and eligibility counts.
[Smoke receipt](log/335-smoke.json) · [Formal receipt](log/336-formal.json) ·
[Run ledger](log/runs.csv) · [Fold scores](log/fold-scores.csv)

Weight 1 had the highest mean Δ PCC (0.431920) and the lowest mean Δ RMSE
(0.249215). Removing both reconstruction terms at weight 0 gave the lowest
mean Δ E-distance (0.020303) and mean Δ MMD (0.143878). At weight 0, mean
Δ PCC was 0.401073 and mean Δ RMSE was 0.266313. Step 05 retained unit
reconstruction weights in its selected objective. Steps 05.4 and 05.5 used
that objective as the reference for testing additional Delta similarity terms.
[Equal-fold scores](figure/data/delta-metrics.csv) ·
[Step 05 selected objective](../README.md#selected-objective) ·
[Step 05.4 PCC sweep](../05.4_mse-l1_pcc/README.md) ·
[Step 05.5 cosine sweep](../05.5_mse-l1_cos/README.md)

The current figures include wJAC-2000. Paired differences subtract the fresh
unit-weight score for the same fold and seed. The learning curves show training
and validation through each fit's stopping epoch.
[Metric scores](figure/metrics.png) ·
[Paired differences](figure/metrics-paired.png) ·
[Learning curves](figure/learning-curve.png) ·
[Per-epoch objectives](figure/data/epoch-objectives.csv) ·
[Render receipt](log/figure-render.json) ·
[CSV-only verification](AGENTS.md#verify-the-current-figures) ·
[Rebuild from accepted fits](AGENTS.md#rebuild-the-current-figures)

## Run the bounded comparison

Follow the [bounded comparison procedure](AGENTS.md#run-the-bounded-comparison)
to reproduce the two-epoch smoke result on fold 0.

## Reproduce the original formal sweep

The formal receipt pins the installed `output/formal/matrix` and
`output/formal/report` manifests and both resource captures. The coordinator
recorded 125 samples with a maximum gap of 10.002 seconds. The scoped
runner recorded 124 samples with a maximum gap of 13.232 seconds.
[Coordinator samples](log/336-coordinator-resources.jsonl) ·
[Scoped samples](log/336-resources.jsonl) ·
[TranSiGen source manifest](../../../pert_model/TranSiGen/manifest.sha256)

[Formal matrix and report procedure](AGENTS.md#reproduce-the-original-formal-sweep) ·
[Replay a sealed fit at its original revision](AGENTS.md#replay-a-sealed-fit) ·
[Registration verification](../../util/training/verification.md) ·
[Original publication commands](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/05_tsg_loss/05.3_mse-l1_mse-x/README.md#run-the-formal-sweep)
