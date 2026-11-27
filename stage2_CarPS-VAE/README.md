# Stage 2 perturbation prediction

<a name="selected-gene-axis"></a>
<a name="transigen-delta-loss-comparison"></a>

## Selected CarPS-VAE configuration

CarPS-VAE uses the Step 06.4 `XΔ-XΔ` architecture selected in
[ADR 0099](../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md).
Separate 100-dimensional Gaussian encoders map measured same-plate DMSO
(`X1`) and measured signed response (`Δ = X2 − X1`). The native Gaussian
prediction head uses the control latent and the frozen UniMol features selected
in [ADR 0100](../docs/adr/0100-carps-vae-unimol.md). Step 10 selects
`no-dose` for the model input. One signed Δ decoder reconstructs measured Δ and
decodes the predicted latent. The treated prediction is `X1 + Δ̂`.
Both measured and predicted Δ use the same matched control on the ordered
Denoise-DEG 5k (`denoise-5k`) axis selected in
[ADR 0102](../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md).
The Step 09 cold-drug axis comparison used training-fold standardized log dose
and gave mean full-axis Δ PCC of 0.4413 ± 0.0136 and Δ RMSE of 0.2386 ± 0.0037.
[Step 09 results](09_cpsv_axis/README.md#formal-results) ·
[Original Step 09 publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/09_cpsv_axis/README.md#formal-results)

Step 10 evaluated joint generalization to held-out drug parents and measured
doses. The selected `no-dose` arm had mean held-out full-axis Δ PCC of
0.3913 ± 0.0235 and Δ RMSE of 0.2472 ± 0.0046 across nine fold and seed
coordinates. [Step 10 results](10_cpsv_dose/README.md#formal-results) ·
[ADR 0104](../docs/adr/0104-carps-vae-no-dose.md)

`L = MSE_X1 + MSE_Δrec + MSE_Δpred + 1.75 × L1_Δpred + 0.1 × (KL(q1 ‖ N) + KL(qΔ ‖ N) + KL(qpred ‖ qΔ))`

`N` is the standard normal prior. The two Δ prediction terms compare
the decoded prediction with measured Δ. Across the nine matched fold
and seed coordinates, the selected arm gave the lowest mean Δ RMSE
(0.240752) and mean Δ PCC of 0.428830. The independent measured Δ decoder
had mean Δ PCC of 0.429835. The standard deviation of the paired differences
exceeded the absolute mean difference on both metrics.

[Ordered genes](02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv) ·
[Tenv training release](03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Architecture evidence](06_tsg_shape/06.4_round4/README.md#formal-results)

## Active experiments and releases

| Step | Design and result |
| --- | --- |
| [01 Raw data](01_data_raw/README.md) | Publishes seven source-native count datasets with gene coverage, observation roles and matched-control provenance. |
| [02.1 HVG-DEG axes](02_gene_axis/02.1_hvg-deg/README.md) | Combines tenv variability with disease DEG1000 and publishes four nested axes and disease overlap evidence. |
| [02.2 Denoise-DEG axes](02_gene_axis/02.2_denoise-deg/README.md) | Filters tenv genes on detection and DMSO variance. Ranks retained genes by tenv and disease DEG presence and publishes five nested axes. |
| [03 Training data](03_data_train/README.md) | Publishes eight tenv plate-pooled releases with measured matched controls and the same frozen folds, plus the Tahoe well-level release. The tenv axes are [03.1 HVG-DEG 5k](03_data_train/03.1_tenv_hvg-5k/README.md), [03.3 HVG-DEG 3k](03_data_train/03.3_tenv_hvg-3k/README.md), [03.4 HVG-DEG 4k](03_data_train/03.4_tenv_hvg-4k/README.md), [03.5 Denoise-DEG 3k](03_data_train/03.5_tenv_denoise-3k/README.md), [03.6 Denoise-DEG 4k](03_data_train/03.6_tenv_denoise-4k/README.md) [03.7 Denoise-DEG 5k](03_data_train/03.7_tenv_denoise-5k/README.md), [03.8 HVG-DEG 2k](03_data_train/03.8_tenv_hvg-2k/README.md) and [03.9 Denoise-DEG 2k](03_data_train/03.9_tenv_denoise-2k/README.md). [03.2 Tahoe](03_data_train/03.2_tahoe_hvg-5k/README.md) uses HVG-DEG 5k. |
| [04 Axis comparison](04_tsg_axis/README.md#formal-results) | Fits TranSiGen on six tenv axes under matched controls and the native objective. Across three cold-drug folds, Denoise-DEG 5k had the highest mean Δ PCC (0.4066) and lowest mean Δ RMSE (0.2596) and Δ MMD (0.1793). HVG-DEG variants led signed DEG F1 for the four disease references. |
| [05.1 Delta loss forms](05_tsg_loss/05.1_loss_metrics/README.md) | Compares MSE, cosine, PCC and L1 as the Δ term at three coefficients. L1 at 1.25 led this grid on Δ PCC (0.423530) and Δ RMSE (0.254063). MSE at coefficient 1 reproduced the native Step 04 result. |
| [05.2 Additive Delta L1](05_tsg_loss/05.2_mse_l1/README.md) | Adds a weighted Δ L1 term to the native objective. The selected 1.75 weight led Δ PCC and Δ RMSE. Weight 1.5 gave marginally lower Δ E-distance and Δ MMD. |
| [05.3 Reconstruction MSE](05_tsg_loss/05.3_mse-l1_mse-x/README.md) | Varies the shared reconstruction MSE weight at fixed Δ L1 weight 1.75. Unit reconstruction weights led Δ PCC and Δ RMSE. Zero reconstruction weight gave lower Δ E-distance (0.020303) and Δ MMD (0.143878). |
| [05.4 Delta PCC](05_tsg_loss/05.4_mse-l1_pcc/README.md) | Adds a bounded Δ PCC loss to the selected objective. Weight zero led Δ PCC and Δ RMSE. |
| [05.5 Delta cosine](05_tsg_loss/05.5_mse-l1_cos/README.md) | Adds a bounded Δ cosine loss to the selected objective. Weight 0.3 increased JAC_neg100 from 0.183801 to 0.185508. Weight zero led Δ PCC and Δ RMSE. |
| [06.1 Expression shape](06_tsg_shape/06.1_round1/README.md#formal-results) | Compares ten encoder, decoder and latent shapes at seed 42. A 4,000-unit encoder gave the highest mean Δ PCC (0.433229) and lowest mean Δ RMSE (0.248576). Deepening both sides gave the lowest Δ E-distance and Δ MMD. The original widths remained the reference for the later prediction-head study. |
| [06.2 Prediction architecture](06_tsg_shape/06.2_round2/README.md#formal-results) | Compares eight prediction architectures on three folds at seed 42. FiLM gave the highest mean Δ PCC (0.439568). An independent signed Δ readout gave the lowest mean Δ RMSE (0.241952). The two modifications became candidates for the combination study. |
| [06.3 Architecture combinations](06_tsg_shape/06.3_round3/README.md#formal-results) | Compares signed Δ readouts and normalized FiLM across three folds and four seeds. The independent signed Δ readout improved both primary means relative to the original. Adding normalized FiLM gave a lower mean Δ PCC and higher mean Δ RMSE than that readout alone. Step 06.4 retained the original Gaussian prediction head. |
| [06.4 Measured encoder targets](06_tsg_shape/06.4_round4/README.md#formal-results) | Compares measured treated-expression and signed Δ encoders and signed Δ decoder sharing across three folds and three seeds. The shared signed Δ decoder had the lowest Δ RMSE, and the independent decoder had the highest Δ PCC. CarPS-VAE selects the shared decoder as recorded in ADR 0099. |
| [07 CarPS-VAE alpha](07_cpsv_alpha/README.md) | Freezes the shared signed Δ model and KPGT default. All nine selected checkpoints exactly replayed the validation objective and predictions and held-out predictions. |
| [08 Chemical representation comparison](08_cpsv_chem-emb/README.md) | Completed 63 matched fold and seed coordinates across seven representations. UniMol led both Δ PCC and Δ RMSE and is selected for subsequent CarPS-VAE studies. |
| [09 CarPS-VAE prediction axes](09_cpsv_axis/README.md#formal-results) | Completed 72 UniMol coordinates across eight axes. Denoise-DEG 5k led the full-axis mean Δ PCC (0.4413) and Δ RMSE (0.2386) and was selected. |
| [10 CarPS-VAE dose inputs](10_cpsv_dose/README.md) | Completed 90 cold-dose-drug coordinates across nine CarPS-VAE configurations and a separate TranSiGen reference. Selected `no-dose` for its highest mean Δ PCC and lowest mean Δ RMSE. |

## Verify the current figures

Each of the 13 comparisons publishes `learning-curve.png`, `metrics.png`
and `metrics-paired.png`. The two metrics figures use eight panels in
the same order. They show every defined fold–seed replicate, the mean,
the sample SD for groups with at least two scores, and zero. Paired
figures show the arm score minus its reference within each fold and seed.
[Comparison figure standard](../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md)

The verifier checks the current receipts and independently redraws each
of the 45 PNGs twice from its published plot data using the `arc` environment.
This covers the 39 comparison PNGs, five Step 02 figures and the palette preview.
For byte-identical replay, use the exact Python, Matplotlib, NumPy and,
where listed, pandas versions in the receipts’ `environment` fields.

```bash
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/stage2-figure-check.XXXXXX")
python stage2_CarPS-VAE/util/figure/script/verify_report_figures.py --work "$WORK/replay"
```

[Step 02 receipt](02_gene_axis/log/342-figure-render.json) ·
[Step 04 receipt](04_tsg_axis/log/figure-render.json) ·
[Step 05.1 receipt](05_tsg_loss/05.1_loss_metrics/log/figure-render.json) ·
[Step 05.2 receipt](05_tsg_loss/05.2_mse_l1/log/figure-render.json) ·
[Step 05.3 receipt](05_tsg_loss/05.3_mse-l1_mse-x/log/figure-render.json) ·
[Step 05.4 receipt](05_tsg_loss/05.4_mse-l1_pcc/log/figure-render.json) ·
[Step 05.5 receipt](05_tsg_loss/05.5_mse-l1_cos/log/figure-render.json) ·
[Step 06.1 receipt](06_tsg_shape/06.1_round1/log/figure-render.json) ·
[Step 06.2 receipt](06_tsg_shape/06.2_round2/log/figure-render.json) ·
[Step 06.3 receipt](06_tsg_shape/06.3_round3/log/figure-render.json) ·
[Step 06.4 receipt](06_tsg_shape/06.4_round4/log/figure-render.json) ·
[Step 08 receipt](08_cpsv_chem-emb/log/figure-render.json) ·
[Step 09 receipt](09_cpsv_axis/log/figure-render.json) ·
[Step 10 receipt](10_cpsv_dose/log/figure-render.json) ·
[ADR 0093](../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md)

The comparison tests inspect the rendered panels against their CSV
values, including panel order, replicate points and sample SD endpoints.

```bash
TMPDIR="$WORK" MPLCONFIGDIR="$WORK/replay/matplotlib" \
PYTHONPATH=stage2_CarPS-VAE/util/output/api \
  python -m unittest discover -s stage2_CarPS-VAE/util/execution/tests -p test_comparison.py
```

Step 07 links the Step 06.4 figures for its fixed alpha. Step 09 verifies
its installed full-axis report with the sealed fits and full-axis inputs.
[Step 07 evidence](07_cpsv_alpha/README.md) ·
[Step 09 reproduction](09_cpsv_axis/README.md#reproduce-and-verify) ·
[Original comparison publication](https://github.com/Moke110/carps/tree/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE)

## Historical evidence and active-input checks

[ADR 0105](../docs/adr/0105-stage-cleanup.md) retires the local Stage 2 archive.
Use the original Git revisions for tracked historical code, figures and receipts.
The [pinned ADR 0087 cutover](https://github.com/Moke110/carps/blob/99078e6f/stage2/README.md#deploy-the-adr-0087-cutover)
and [original CarPS path migration](https://github.com/Moke110/carps/blob/ac62a177a54aa99e6ac626728af12394f3e0d940/stage2_CarPS-VAE/README.md#deploy-the-carps-stage-2-path)
document their original deployments.

With Steps 01–03 materialized, run
`python stage2_CarPS-VAE/04_tsg_axis/execution/verify_active_inputs.py`
from a checkout without the local archive. It checks active corpus windows,
ordered projections, training folds and both tenv and Tahoe training readers.
Use the [pinned Step 06.4 replay](06_tsg_shape/06.4_round4/README.md#replay-a-sealed-step-064-fit)
for registrations whose source keys start with `stage2/`.
