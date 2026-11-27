# Step 05.2 · TranSiGen additive Delta L1 sweep

Step 05.2 measured the effect of adding a Delta L1 penalty to the native
TranSiGen objective. It swept nine weights from 0 to 2 in 0.25 increments,
using the three frozen tenv cold-drug folds, seed 42, frozen KPGT vectors and
the Denoise-DEG 5k training release. The measured same-plate DMSO profile
served as model input and as the common reference for predicted and observed
Delta. [Spec 10](https://github.com/Moke110/carps/issues/331) ·
[ADR 0090](../../../docs/adr/0090-transigen-mse-l1-sweep.md) ·
[Training release manifest](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json)

The objective was `0.1 * K + MSE_X1 + MSE_X2 + MSE_Δ + alpha * L1_Δ`, where
`K` sums the three native KL terms. At zero weight, the runner used native
arithmetic and still recorded the raw Delta L1 component. Each fit selected
the earliest minimum of the registered validation objective with 30-epoch
patience and a 1,000-epoch cap. Held-out scoring used the Step 05.1 metrics and
frozen HVG-DEG 5k bandwidths. [Step 05.1 report](../05.1_loss_metrics/README.md) ·
[Frozen bandwidths](log/bandwidths.csv)

## Accepted results and follow-up

The accepted sweep contains 27 independently replayed fits, 243 held-out
fold scores and 81 equal-fold aggregates. The three fresh zero-weight fits
exactly matched the corresponding Step 05.1 native MSE fits at coefficient
1 on native objective traces, selected and stopping epochs, validation and
test predictions, Delta arrays, fold scores and eligibility counts.
[Formal receipt](log/333-formal.json) · [Run ledger](log/runs.csv) ·
[Fold scores](log/fold-scores.csv)

Weight 1.75 had the highest mean Δ PCC (0.431920 versus 0.406610 at zero)
and lowest mean Δ RMSE (0.249215 versus 0.259575 at zero). Weight 1.5 had the
lowest mean Δ E-distance and Δ MMD. Step 05.3 fixed the Delta L1 weight at
1.75 and varied the reconstruction MSE weight, exactly reproducing the
three Step 05.2 fits at its unit-weight reference.
[Equal-fold scores](figure/data/delta-metrics.csv) ·
[Step 05.3 report](../05.3_mse-l1_mse-x/README.md) ·
[Step 05 selected objective](../README.md#selected-objective)

The current figures show eight metrics, including wJAC-2000. Bars show
means, dots mark individual folds and error bars show sample SD. Paired
scores subtract the zero-weight score within each fold and seed. The learning
curves cover training and validation through each fit's stopping epoch.
[Metric scores](figure/metrics.png) · [Paired differences](figure/metrics-paired.png) ·
[Learning curves](figure/learning-curve.png) · [Per-epoch objectives](figure/data/epoch-objectives.csv) ·
[Render receipt](log/figure-render.json) ·
[CSV-only verification](AGENTS.md#verify-the-current-figures)

## Run a bounded comparison

The smoke run fit weights 0 and 1 on fold 0. Both fits passed exact replay of
validation predictions, the validation objective and test predictions.
[Two-fit smoke receipt](log/332-smoke.json) ·
[Bounded comparison procedure](AGENTS.md#run-a-bounded-comparison)

## Reproduce the original formal sweep

The matrix registered all 27 coordinates before fitting and checked ordered
fold membership. The formal receipt pins the installed `output/formal/matrix`
and `output/formal/report` manifests and both resource captures. The matrix
recorded peak RSS by phase and sampled WSL RAM, GPU and disk. Its 137 scoped
samples had 20 intervals longer than 12 seconds and a maximum gap of 19.951
seconds. The independent coordinator recorded 137 samples across the same
period with a maximum gap of 10.147 seconds.
[Scoped resources](log/333-resources.jsonl) ·
[Coordinator resources](log/333-coordinator-resources.jsonl) ·
[TranSiGen source manifest](../../../pert_model/TranSiGen/manifest.sha256)

[Formal matrix and report procedure](AGENTS.md#reproduce-the-original-formal-sweep) ·
[Replay a sealed fit at its original revision](AGENTS.md#replay-a-sealed-fit) ·
[Registration verification](../../util/training/verification.md)
