# Stage 2 perturbation prediction

<a name="selected-gene-axis"></a>
<a name="transigen-delta-loss-comparison"></a>

## Selected TranSiGen configuration

Use Denoise-DEG 5k (`denoise-5k`) and measured same-plate DMSO as model
input `X1` during training, validation and test. Compute both predicted and
measured Δ by subtracting that same control from the corresponding treated
profile. The reported results assume that a measured matched control is
available when making a prediction.

`L = 0.1 * K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 * L1_Δ`

`K` sums the three native KL terms. Both expression reconstruction MSEs and
the Δ MSE retain their unit weights. This configuration led the published
Step 05 fits on mean Δ PCC and Δ RMSE across the three frozen tenv cold-drug
folds. At Delta L1 weight 1.75, mean Δ PCC reached 0.431920 and mean Δ RMSE
fell to 0.249215, compared with 0.406610 and 0.259575 at weight zero.
[Ordered genes](02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv) ·
[Tenv training release](03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Step 05 conclusion](05_tsg_loss/README.md#selected-objective)

## Active experiments and releases

| Step | Design and result |
| --- | --- |
| [01 Raw data](01_data_raw/README.md) | Publishes seven source-native count datasets with gene coverage, observation roles and matched-control provenance. |
| [02 Prediction axes](02_gene_axis/README.md) | HVG-DEG combines tenv variability with disease DEG1000. Denoise-DEG uses tenv detection and DMSO variance filters, then ranks retained genes by tenv and disease DEG presence. Both publish nested ordered axes and disease overlap evidence. |
| [03 Training data](03_data_train/README.md) | Publishes six tenv training releases with plate-pooled counts, matched controls and the same frozen folds, plus Tahoe well-level data on HVG-DEG 5k. |
| [04 Axis comparison](04_tsg_axis/README.md#formal-results) | Fits TranSiGen on six tenv axes under matched controls and the native objective. Across three cold-drug folds, Denoise-DEG 5k had the highest mean Δ PCC (0.4066) and lowest mean Δ RMSE (0.2596) and Δ MMD (0.1793). HVG-DEG variants led signed DEG F1 for the four disease references. |
| [05.1 Delta loss forms](05_tsg_loss/05.1_loss_metrics/README.md) | Compares MSE, cosine, PCC and L1 as the Δ term at three coefficients. L1 at 1.25 led this grid on Δ PCC (0.423530) and Δ RMSE (0.254063). MSE at coefficient 1 reproduced the native Step 04 result. |
| [05.2 Additive Delta L1](05_tsg_loss/05.2_mse_l1/README.md) | Adds a weighted Δ L1 term to the native objective. The selected 1.75 weight led Δ PCC and Δ RMSE. Weight 1.5 gave marginally lower Δ E-distance and Δ MMD. |
| [05.3 Reconstruction MSE](05_tsg_loss/05.3_mse-l1_mse-x/README.md) | Varies the shared reconstruction MSE weight at fixed Δ L1 weight 1.75. Unit reconstruction weights led Δ PCC and Δ RMSE. Zero reconstruction weight gave lower Δ E-distance (0.020303) and Δ MMD (0.143878). |
| [05.4 Delta PCC](05_tsg_loss/05.4_mse-l1_pcc/README.md) | Adds a bounded Δ PCC loss to the selected objective. Weight zero led Δ PCC and Δ RMSE. |
| [05.5 Delta cosine](05_tsg_loss/05.5_mse-l1_cos/README.md) | Adds a bounded Δ cosine loss to the selected objective. Weight 0.3 increased JAC_neg100 from 0.183801 to 0.185508. Weight zero led Δ PCC and Δ RMSE. |
| [06.1 Expression shape](06_tsg_shape/06.1_round1/README.md#formal-results) | Compares ten encoder, decoder and latent architectures over three frozen cold-drug folds. Six new fits extend the encoder-width series, and the accepted original fit on each fold remains the paired reference for nine held-out metrics. |
| [06.2 Prediction architecture](06_tsg_shape/06.2_round2/README.md#formal-results) | Compares eight prediction architectures on the three frozen folds. The three accepted Step 06.1 original fits supply the paired references for 21 new fits. |

## Verify the current figures

All 28 active report figures and the palette preview use the [Stage 2 palette and drawing API](util/README.md#plotting).
The verifier checks the current render receipts and reproduces each PNG
twice from its published plot data using the `arc` rendering environment.
For byte-identical replay, use the exact Python, Matplotlib, NumPy and,
where listed, pandas versions in the receipts’ `environment` fields.

```bash
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/stage2-figure-check.XXXXXX")
python stage2/util/figure/script/verify_report_figures.py --work "$WORK/replay"
```

[Step 02 receipt](02_gene_axis/log/342-figure-render.json) ·
[Step 04 receipt](04_tsg_axis/log/342-figure-render.json) ·
[Step 05 receipt](05_tsg_loss/log/343-figure-render.json) ·
[Step 06.1 receipt](06_tsg_shape/06.1_round1/log/figure-render.json) ·
[Step 06.2 receipt](06_tsg_shape/06.2_round2/log/figure-render.json) ·
[ADR 0093](../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md)

## Local history

[ADR 0087](../docs/adr/0087-stage2-active-steps-and-history.md) defines the archive
cutover. The materialized checkout keeps superseded raw and training releases,
prior comparisons, the archived Step 05–09 studies and rollback materials in
`stage2/history`. The local `archive-manifest.json` records original
paths, archived paths and replay revisions. The duplicate render bundle
for the current Step 04 report was removed from the installed archive.
The original cutover recorded the inventory at its deployment revision.
Normal readers, builds and
comparisons use the active releases. In an archive-isolated checkout with the
active inputs materialized, run `python stage2/04_tsg_axis/execution/verify_active_inputs.py`
to verify the Step 01–03 readers, projections, and fold population.

The [cutover smoke receipt](04_tsg_axis/log/323-archive-isolated-smoke.json)
records the archive-isolated acceptance.

The [pinned archive cutover instructions](https://github.com/Moke110/arc/blob/99078e6f/stage2/README.md#deploy-the-adr-0087-cutover)
record the migration procedure for its original deployment revision.
