# Step 04 TranSiGen prediction-axis comparison

Step 04 compared HVG-DEG and Denoise-DEG at 3k, 4k and 5k genes on the
same tenv plate–parent–dose conditions. The comparison uses the shared five
parent-disjoint cold-drug folds and five training seeds 17, 42, 43, 2026 and
2027. Each axis sets the expression input and output width. All 150 fits started from random weights and used
the original TranSiGen internal widths, frozen 2,304-dimensional KPGT vectors
and a fold-fitted numeric dose transform on the latent input.
[Shared-protocol rerun](https://github.com/Moke110/carps/issues/399) ·
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

The report averages defined conditions within each fold and seed, then gives
equal weight to defined fold and seed coordinates within each axis. Signed
Jaccard ranks strictly positive and strictly negative genes separately, with
ties resolved by
frozen gene index. E-distance and RBF MMD use
Euclidean distance divided by the square root of axis width and include
the empirical diagonal terms. All six axes share one bandwidth per fold,
estimated as the median distance between distinct conditions in the full
HVG-DEG 5k training measured Delta. The report records the bandwidths and
their training membership in `bandwidths.csv` and `fold-memberships.csv`.
[Formal receipt](log/399-tenv-v1-formal.json)

For each disease gene set, signed DEG F1 compares model-axis predictions
with nonzero measured directions within its DEG1000 list. The reference
directions come from raw pooled counts and matched plate controls on the
full common CM axis and are shared across all six axes. An opposite sign
counts as one false positive and one false negative. An off-axis reference
gene counts as a false negative. Two empty sets give an undefined score.
[Disease reference implementation](output/api/tsg_axis/disease_reference.py) ·
[Shared metric implementation](../util/eval/metrics.py)

## Formal results

All 150 fits passed exact validation loss, validation prediction and held-out
prediction replay. The independent report run reproduced the scientific
tables byte for byte. Each axis contributed 25 coordinates to the absolute
metrics and the within-fold and seed pairing against Denoise-DEG 5k.
[Formal receipt](log/399-tenv-v1-formal.json)

Folds 0–4 test 276, 276, 279, 276, 276 conditions, respectively.
The full-axis RBF bandwidths are 0.433788, 0.435847, 0.431563, 0.427257, 0.430072.
Fits ran for 118–227 epochs.

| Delta metric | `hvg-3k` | `hvg-4k` | `hvg-5k` | `denoise-3k` | `denoise-4k` | `denoise-5k` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Δ Spearman | 0.3906 | 0.3968 | 0.3968 | 0.3887 | 0.3908 | 0.3919 |
| JAC_top200 | 0.1468 | 0.1424 | 0.1388 | 0.1480 | 0.1429 | 0.1415 |
| Δ PCC | 0.3897 | 0.3942 | 0.3921 | 0.4005 | 0.4023 | 0.4016 |
| Δ RMSE | 0.3584 | 0.3193 | 0.2954 | 0.3088 | 0.2846 | 0.2584 |
| Δ R² | -0.2301 | -0.2020 | -0.2017 | -0.1701 | -0.1646 | -0.1501 |
| Δ E-distance | 0.0359 | 0.0328 | 0.0323 | 0.0301 | 0.0297 | 0.0279 |
| Δ MMD | 0.2184 | 0.1983 | 0.1889 | 0.1866 | 0.1786 | 0.1636 |
| JAC_pos100 | 0.1412 | 0.1336 | 0.1255 | 0.1309 | 0.1225 | 0.1188 |
| JAC_neg100 | 0.1634 | 0.1625 | 0.1630 | 0.1709 | 0.1701 | 0.1716 |

| Signed disease DEG F1 | `hvg-3k` | `hvg-4k` | `hvg-5k` | `denoise-3k` | `denoise-4k` | `denoise-5k` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nHF | 0.1606 | 0.1637 | 0.1899 | 0.1205 | 0.1071 | 0.0972 |
| ICM | 0.1084 | 0.1259 | 0.1729 | 0.1217 | 0.1064 | 0.0964 |
| genetic_cardiomyopathy | 0.2136 | 0.2465 | 0.2185 | 0.1214 | 0.1065 | 0.0952 |
| TOX | 0.3027 | 0.2587 | 0.2152 | 0.1306 | 0.1154 | 0.1039 |

HVG-DEG 5k had the highest mean Δ Spearman (0.396780), close to HVG-DEG 4k
(0.396765). Denoise-DEG 3k led JAC-top200 (0.147964), and Denoise-DEG 4k led
Δ PCC (0.402325). Denoise-DEG 5k had the lowest Δ RMSE (0.258378). The
rerun uses Denoise-DEG 5k as the established reference for successor steps.

## Gene axis selection

In the historical three-fold comparison, **Denoise-DEG 5k (`denoise-5k`)** led mean held-out Δ PCC (0.4066) and
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
| [Learning curves](figure/learning-curve.png) | Training and validation Loss traces for all 150 fits.
[Plotted data](figure/data/learning-curve.csv) |
| [Held-out metrics](figure/metrics.png) | Fold and seed scores, means and sample SDs.
[Plotted data](figure/data/metrics.csv) |
| [Paired differences](figure/metrics-paired.png) | Each other axis minus Denoise-DEG 5k within the same fold and seed.
[Plotted data](figure/data/metrics-paired.csv) |

The seven metric panels include zero in both figures. The current figure receipt binds all
three CSVs and three PNGs to the sealed fit and report manifests. The
historical publication retains the disease F1 figure.
[Current panel layout](../../docs/adr/0109-cpst-seven-metric-comparison-figures.md) ·
[Current rendering receipt](log/figure-render.json) ·
[Historical disease F1 figure](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/04_tsg_axis/figure/disease-f1.png)

The successor is installed under `output/formal-tenv-v1/` and `output/report-tenv-v1/`.
The historical receipt records the checkpoint, registration and prediction
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
