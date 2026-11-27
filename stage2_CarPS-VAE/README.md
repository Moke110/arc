# Stage 2 perturbation prediction

<a name="selected-gene-axis"></a>
<a name="transigen-delta-loss-comparison"></a>

## Selected CarPS-VAE configuration

CarPS-VAE uses the Step 06.4 `XΔ-XΔ` architecture selected in
[ADR 0099](../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md).
Separate 100-dimensional Gaussian encoders map measured same-plate DMSO
(`X1`) and measured signed response (`Δ = X2 − X1`). The native Gaussian
prediction head uses the control latent, the frozen UniMol features selected
in [ADR 0100](../docs/adr/0100-carps-vae-unimol.md), and dose standardized
on the training fold. One signed Δ decoder reconstructs measured Δ and
decodes the predicted latent. The treated prediction is `X1 + Δ̂`.
Both measured and predicted Δ use the same matched control on the ordered
Denoise-DEG 5k (`denoise-5k`) axis selected in
[ADR 0102](../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md).
The Step 09 full-axis comparison gave mean Δ PCC of 0.4413 ± 0.0136
and mean Δ RMSE of 0.2386 ± 0.0037 for this axis.
[Step 09 results](09_cpsv_axis/README.md#formal-results) ·
[Original Step 09 publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/09_cpsv_axis/README.md#formal-results)

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

## Verify the current figures

All 45 active report figures and the palette preview use the [Stage 2 palette and drawing API](util/README.md#plotting).
The verifier checks the current render receipts and reproduces each PNG
twice from its published plot data using the `arc` rendering environment.
For byte-identical replay, use the exact Python, Matplotlib, NumPy and,
where listed, pandas versions in the receipts’ `environment` fields.

```bash
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/stage2-figure-check.XXXXXX")
python stage2_CarPS-VAE/util/figure/script/verify_report_figures.py --work "$WORK/replay"
```

[Step 02 receipt](02_gene_axis/log/342-figure-render.json) ·
[Step 04 receipt](04_tsg_axis/log/342-figure-render.json) ·
[Step 05 receipt](05_tsg_loss/log/343-figure-render.json) ·
[Step 06.1 receipt](06_tsg_shape/06.1_round1/log/figure-render.json) ·
[Step 06.2 receipt](06_tsg_shape/06.2_round2/log/figure-render.json) ·
[Step 06.3 receipt](06_tsg_shape/06.3_round3/log/figure-render.json) ·
[Step 06.4 receipt](06_tsg_shape/06.4_round4/log/figure-render.json) ·
[Step 08 receipt](08_cpsv_chem-emb/log/figure-render.json) ·
[Step 09 receipt](09_cpsv_axis/log/figure-render.json) ·
[ADR 0093](../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md)

## Local history

[ADR 0087](../docs/adr/0087-stage2-active-steps-and-history.md) defines the archive
cutover. The materialized checkout keeps superseded raw and training releases,
prior comparisons, the archived Step 05–09 studies and rollback materials in
`stage2_CarPS-VAE/history`. The local `archive-manifest.json` records original
paths, archived paths and replay revisions. The duplicate render bundle
for the current Step 04 report was removed from the installed archive.
The manifest retains the original `stage2/` archive-relative paths. Normal readers, builds and comparisons use the active releases. In an archive-isolated checkout with the
active inputs materialized, run `python stage2_CarPS-VAE/04_tsg_axis/execution/verify_active_inputs.py`
to verify the Step 01–03 readers, projections, and fold population.

The [cutover smoke receipt](04_tsg_axis/log/323-archive-isolated-smoke.json)
records the archive-isolated acceptance.

The [pinned archive cutover instructions](https://github.com/Moke110/carps/blob/99078e6f/stage2/README.md#deploy-the-adr-0087-cutover)
record the migration procedure for its original deployment revision.

## Deploy the CarPS Stage 2 path

[ADR 0099](../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md)
moves the materialized `stage2/` tree to `stage2_CarPS-VAE/`. Set the Git
remote to the canonical repository with
`git remote set-url origin https://github.com/Moke110/carps`.

For the first fast-forward across the rename, execute the migration script
from an integration worktree that contains the successor code. Run from the
materialized `main` checkout after the integrated deployment passes
acceptance. Set `SUCCESSOR_WT` to the integration worktree path and
`TARGET` to the accepted `dev` tip.

```bash
ROOT=$(git rev-parse --show-toplevel)
python "$SUCCESSOR_WT/script/migrate_carps_stage2.py" --root "$ROOT" \
  --work "$ROOT/tmp/carps-stage2-cutover" deploy --target "$TARGET"
```

The script parks the old directory on the checkout filesystem while Git
fast-forwards, renames it to the canonical path and installs changed tracked
files from the target tree. It records the pre/post inventory under the
specified work directory. It verifies the identities of preserved files and
directories, hashes small files and checks each external symlink target.
Use the same arguments with `resume` if interrupted after a plan has been
recorded. Run `verify` to compare the just-installed checkout against that
plan. Use `resume` and `verify` before installing the new Step 08
`output/formal` bundle. After installing the bundle, run the
[Step 08 scientific verifier](08_cpsv_chem-emb/README.md#formal-reproduction).
For a later checkout-root move, repair the Git metadata for each
relocated linked worktree with `git worktree repair <linked-worktree-path>`
from the moved primary checkout. Working-tree payload locations and Git
worktree registrations are separate.

Use the [pinned Step 06.4 replay](06_tsg_shape/06.4_round4/README.md#replay-a-sealed-step-064-fit)
for registrations whose recorded code keys start with `stage2/`.
