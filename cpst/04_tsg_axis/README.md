# Step 04 TranSiGen prediction-axis comparison

Step 04 compared HVG-DEG and Denoise-DEG at 3k, 4k and 5k genes on the
same tenv plate–parent–dose conditions. The comparison reused the frozen
three cold-drug folds from Step 03. Each axis set the expression input and
output width. All 18 fits started from random weights at seed 42 with the
original TranSiGen internal widths, frozen 2,304-dimensional KPGT vectors
and a fold-fitted numeric dose transform on the latent input.
[Spec 6](https://github.com/Moke110/carps/issues/317) ·
[Step 03 tenv training axes](../03_data_train/README.md#tenv-training-axes)

## Matched-plate protocol

Training, validation and test each used the measured same-plate DMSO profile
from the release. Training and checkpoint selection used TranSiGen's native
MSE + KLD objective with β = 0.1. The selected checkpoint minimized the
validation objective summed over batches and divided by the number of
validation conditions. Validation ran in evaluation mode with a fixed seed.
The fit stopped after 30 epochs without improvement or at 1,000 epochs.
[ADR 0086](../../docs/adr/0086-tenv-axis-matched-plate-control.md)

The model produced absolute treated expression in `prediction.npy`.
Predicted and measured Delta each subtract the same matched plate control,
with rows and columns ordered by `conditions.csv` and `genes.csv`.
Measured test treatment entered only at scoring. The comparison estimates
performance for published axes constructed with tenv evidence, conditional
on having a measured same-plate control at prediction time.
[Fixed-prediction input interface](output/api/tsg_axis/inputs.py)

The report averages defined conditions within each fold and gives equal
weight to defined folds within each axis. Signed Jaccard ranks strictly
positive and strictly negative genes separately, with ties resolved by
frozen gene index. E-distance and RBF MMD use
Euclidean distance divided by the square root of axis width and include
the empirical diagonal terms. All six axes share one bandwidth per fold,
estimated as the median distance between distinct conditions in the full
HVG-DEG 5k training measured Delta. The report records the bandwidths and
their training membership in `bandwidths.csv` and `fold-memberships.csv`.
[Formal receipt](log/319-formal.json)

For each disease gene set, signed DEG F1 compares model-axis predictions
with nonzero measured directions within its DEG1000 list. The reference
directions come from raw pooled counts and matched plate controls on the
full common CM axis and are shared across all six axes. An opposite sign
counts as one false positive and one false negative. An off-axis reference
gene counts as a false negative. Two empty sets give an undefined score.
[Disease reference implementation](output/api/tsg_axis/disease_reference.py) ·
[Shared metric implementation](../util/eval/metrics.py)

## Formal results

The bounded HVG-DEG 5k smoke fit passed matched-control, scoring and figure
replay checks. The formal run then passed exact replay of the validation
objective, validation predictions and held-out predictions on all 18 fits.
Every fit reached the 30-epoch patience limit, with 127 to 243 epochs run.
The six axes shared the same test rosters of 462, 459 and 462 conditions
for folds 0–2. Their RBF bandwidths were 0.430414, 0.425762 and 0.439392.
All nine Delta metrics and all four disease scores were defined on every
fold. The sealed report also reproduced byte for byte.
[Smoke receipt](log/318-matched-control-smoke.json) ·
[Formal receipt](log/319-formal.json)

| Delta metric | `hvg-3k` | `hvg-4k` | `hvg-5k` | `denoise-3k` | `denoise-4k` | `denoise-5k` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| JAC_top200 | 0.1453 | 0.1428 | 0.1416 | 0.1454 | 0.1427 | 0.1409 |
| JAC_pos100 | 0.1375 | 0.1319 | 0.1258 | 0.1283 | 0.1227 | 0.1134 |
| JAC_neg100 | 0.1644 | 0.1650 | 0.1695 | 0.1683 | 0.1700 | 0.1766 |
| Δ PCC | 0.3928 | 0.3967 | 0.4043 | 0.4023 | 0.4003 | 0.4066 |
| Δ Spearman | 0.3896 | 0.3989 | 0.4045 | 0.3906 | 0.3888 | 0.3961 |
| Δ RMSE | 0.3586 | 0.3240 | 0.2957 | 0.3107 | 0.2893 | 0.2596 |
| Δ E-distance | 0.0370 | 0.0394 | 0.0350 | 0.0320 | 0.0315 | 0.0319 |
| Δ MMD | 0.2253 | 0.2208 | 0.1994 | 0.1955 | 0.1873 | 0.1793 |
| Δ R² | -0.1953 | -0.2407 | -0.1842 | -0.1993 | -0.2275 | -0.1587 |

| Signed disease DEG F1 | `hvg-3k` | `hvg-4k` | `hvg-5k` | `denoise-3k` | `denoise-4k` | `denoise-5k` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nHF | 0.1603 | 0.1639 | 0.1905 | 0.1197 | 0.1063 | 0.0969 |
| ICM | 0.1083 | 0.1258 | 0.1734 | 0.1211 | 0.1057 | 0.0961 |
| Genetic cardiomyopathy | 0.2136 | 0.2469 | 0.2193 | 0.1208 | 0.1057 | 0.0950 |
| TOX | 0.3026 | 0.2587 | 0.2162 | 0.1302 | 0.1148 | 0.1036 |

## Gene axis selection

**Denoise-DEG 5k (`denoise-5k`)** led mean held-out Δ PCC (0.4066) and
had the lowest mean Δ RMSE (0.2596) and Δ MMD (0.1793). HVG-DEG variants
led signed DEG F1 for all four disease reference sets. Step 04 prioritized
held-out Delta fidelity under matched controls and selected Denoise-DEG 5k
for the subsequent TranSiGen loss and architecture studies. Step 09 later
compared eight axes with CarPS-Transcriptome and UniMol and selected
Denoise-DEG 5k under that protocol as well.
[Step 03.7 training release](../03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Step 05 loss studies](../05_tsg_loss/README.md) ·
[Step 06 architecture studies](../06_tsg_module/README.md) ·
[Step 09 selection](../09_cpst_axis/README.md#formal-results)

## Figures and evidence

| Figure | Evidence |
| --- | --- |
| [Learning curves](figure/learning-curve.png) | Training and validation traces for all 18 fits.
[Plotted data](figure/data/learning-curve.csv) |
| [Held-out metrics](figure/metrics.png) | Fold scores, means and sample SDs.
[Plotted data](figure/data/metrics.csv) |
| [Paired differences](figure/metrics-paired.png) | Each other axis minus Denoise-DEG 5k within the same fold and seed.
[Plotted data](figure/data/metrics-paired.csv) |

The seven metric panels include zero in both figures. The current figure receipt binds all
three CSVs and three PNGs to the sealed fit and report manifests. The
historical publication retains the disease F1 figure.
[Current panel layout](../../docs/adr/0109-cpst-seven-metric-comparison-figures.md) ·
[Current rendering receipt](log/figure-render.json) ·
[Historical disease F1 figure](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/04_tsg_axis/figure/disease-f1.png)

The formal receipt records the checkpoint, registration and prediction
inventory under `output/formal-matched-control/`, the report tables and
shared disease reference under `output/report-matched-control/` and the
process evidence under `log/319-runs/`. These paths require the materialized
release. The sealed report retains its original axis IDs. The canonical
export rewrites the IDs in six tables while preserving all other fields.
[Formal receipt](log/319-formal.json) ·
[Canonical export receipt](log/322-canonical-axis-report.json)

## Historical comparisons

The earlier single-dose and axis-wide-control protocols, and their dose and
width successors, are documented at a fixed revision.
[Historical comparison protocols](https://github.com/Moke110/carps/blob/59b188ec0fe6386f2c916e062eaad310c20f00f8/stage2/04_tsg_axis/README.md#historical-single-dose-comparison)

## Prepare the capped runtime

Use the owner guide to prepare the capped runtime, run fits and replay the
accepted report and figures.
[Capped runtime](AGENTS.md#prepare-the-capped-runtime) ·
[Run the formal matrix](AGENTS.md#run-the-formal-matrix) ·
[Register and replay one fit](AGENTS.md#register-and-replay-one-fit) ·
[Verify the sealed report](AGENTS.md#verify-the-sealed-report) ·
[Rebuild and verify figures](AGENTS.md#rebuild-and-verify-the-comparison-figures)
