# CarPS-Transcriptome (CPST)

CarPS-Transcriptome predicts the signed gene-expression response of
cardiomyocytes to a chemical perturbation relative to a measured control.
The accepted tenv studies establish the model architecture, objective, gene
axis, chemical representation and dose-input policy. The [CPST guide](AGENTS.md)
routes data preparation, fitting, report reconstruction and verification.

<a name="selected-gene-axis"></a>
<a name="transigen-delta-loss-comparison"></a>

## Selected CarPS-Transcriptome configuration

The current selection uses the `XΔ-XΔ` architecture, frozen UniMol features,
the ordered Denoise-DEG 5k axis (`denoise-5k`) and the `no-dose` prediction
head. Separate 100-dimensional Gaussian encoders model measured same-plate
DMSO (`X1`) and measured signed response (`Δ = X2 − X1`). The prediction
head uses the control latent and a learned projection of the frozen UniMol
features. One signed Δ decoder reconstructs the measured response and decodes
the predicted response. Predicted treated expression equals `X1 + Δ̂`.
[Architecture decision](../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) ·
[Ordered genes](02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv) ·
[Tenv training release](03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Selected no-dose results](10_cpst_dose/README.md#formal-results)

`L = MSE_X1 + MSE_Δrec + MSE_Δpred + 1.75 × L1_Δpred + 0.1 × (KL(q1 ‖ N) + KL(qΔ ‖ N) + KL(qpred ‖ qΔ))`

The first two MSE terms measure control and measured Δ reconstruction. The
prediction MSE and L1 compare predicted Δ with measured Δ. The KL terms
regularize the two encoder distributions against the standard normal prior
`N` and align the predicted Gaussian with the measured Δ Gaussian.
[Objective lineage](05_tsg_loss/README.md#selected-objective) ·
[Shared CPST model](util/training/cpst_model.py) · [Step 07 loss design](07_cpst_loss/README.md)

## Evidence in experiment order

1. **Steps 01–03 established the data and prediction axes.** Step 01 released
   native integer counts, gene coverage and matched-control provenance for
   seven datasets. Step 02 published nested HVG-DEG and Denoise-DEG axes.
   Step 03 materialized eight axis-local, plate-pooled tenv views with measured
   same-plate DMSO and shared frozen folds, plus a separate Tahoe well-level
   view on HVG-DEG 5k. The formal axis comparisons used tenv.
   [Native counts](01_data_raw/README.md#accepted-corpus) ·
   [Axis selection evidence](02_gene_axis/README.md#axis-selection-evidence) ·
   [Training releases](03_data_train/README.md)
2. **Steps 04–05 established the TranSiGen reference axis and loss.** Step 04
   chose Denoise-DEG 5k from six tenv axes because it led held-out Δ PCC,
   Δ RMSE and Δ MMD under the matched-control protocol. HVG-DEG variants
   led signed DEG F1 for the four disease references. The loss studies
   selected predicted Δ MSE plus Δ L1 at weight 1.75, with unit reconstruction
   weights. Zero added weight led both primary metrics in the PCC and cosine
   sweeps.
   [Step 04 selection](04_tsg_axis/README.md#gene-axis-selection) ·
   [Step 05 results](05_tsg_loss/README.md#loss-studies)
3. **Step 06 selected the shared signed Δ architecture.** The four rounds
   tested network shapes, signed Δ readouts and chemical conditioning. Round 3
   accepted normalized FiLM after non-finite prediction KL blocked the
   unnormalized combinations. Round 4 kept the native Gaussian prediction head.
   Its `XΔ-XΔ` arm had the lowest mean Δ RMSE (0.240752). The arm with an
   independent Δ decoder had the highest mean Δ PCC (0.429835 versus 0.428830
   for the shared arm). The sample SD of the paired differences
   exceeded the absolute mean difference on both metrics. The architecture
   decision chose the shared decoder for CarPS-Transcriptome.
   [Four-round evidence](06_tsg_module/README.md) ·
   [Step 06.4 formal results](06_tsg_module/06.4_round4/README.md#formal-results) ·
   [Architecture decision](../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md)
4. **Step 07 compares losses on the shared signed Δ model.** The published
   alpha predecessor retains the
   Step 06.4 `XΔ-XΔ` backbone, Denoise-DEG 5k, frozen KPGT features and
   training-fold standardized log dose. Loading the nine sealed Step 06.4
   checkpoints through the alpha interface exactly reproduced the validation
   objective and predictions, the held-out treated-expression predictions
   and the measured and predicted held-out Δ arrays. This parity evidence
   applies to cold-drug folds 0–2, training seeds 17, 42 and 2026, and
   sampling seeds 22617 for validation and 100042 for the held-out set.
   The pinned alpha release preserves those settings for replay. The current
   Step 07 study registers additive loss terms on the shared tenv protocol.
   [Step 07 loss study](07_cpst_loss/README.md) ·
   [Alpha replay evidence](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/README.md#checkpoint-and-replay-evidence) ·
   [Release manifest](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/output/release.json) ·
   [Parity receipt](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/log/362-alpha.json)
5. **Step 08 compared seven representations, then ADR 0100 adopted UniMol for
   the successor configuration.** UniMol led Δ PCC and Δ RMSE across the
   63 fold and seed coordinates. The prespecified rule for changing the
   working default required the mean paired improvement to exceed its sample
   SD on both metrics. The UniMol Δ PCC
   gain over KPGT was 0.012432 with sample SD 0.013493, so the formal
   comparison retained KPGT. ADR 0100 separately adopted frozen UniMol from
   Step 09 onward.
   [Step 08 formal results](08_cpst_chem-emb/README.md#formal-results) ·
   [UniMol adoption](../docs/adr/0100-carps-vae-unimol.md)
6. **Step 09 selected Denoise-DEG 5k with UniMol.** The comparison covered
   72 coordinates across eight tenv axes using training-fold standardized log
   dose. Denoise-DEG 5k led the mean full-axis Δ PCC (0.4413 ± 0.0136) and
   Δ RMSE (0.2386 ± 0.0037). The selection applies to normalization on each
   released axis. The study reused the nine Step 08 UniMol Denoise-DEG 5k fits
   after parity checks.
   [Step 09 formal results](09_cpst_axis/README.md#formal-results) ·
   [Axis decision](../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md)
7. **Step 10 selected `no-dose` for the current configuration.** On the
   three frozen cold-dose-drug folds, `no-dose` led all comparison arms in
   mean full-axis Δ PCC (0.3913 ± 0.0235) and Δ RMSE (0.2472 ± 0.0046)
   across nine fold and seed coordinates. Dose remains part of measured
   condition identity and the holdout protocol. The Step 10 release
   contains the accepted formal fits and audit evidence for this selection.
   [Step 10 formal results](10_cpst_dose/README.md#formal-results) ·
   [No-dose decision](../docs/adr/0104-carps-vae-no-dose.md)

## Reporting for the tenv rerun

Read the primary held-out comparison metrics in order of priority. First
read Δ Spearman, then signed Jaccard at 200 genes (JAC-top200), then Δ PCC. Use each step's
established baseline model when running the next step. Each step report
documents the protocol used for its release.

## Prospective use in CarPS-Phenotype

The planned phenotype predictor would take ordered signed CarPS-Transcriptome
Δ and assay dose as separate inputs. The completed Morgan benchmark uses
fingerprints and assay dose. Selecting and documenting a Heliste prediction
release, checking coverage by parent and dose, and integrating it with the
frozen phenotype benchmark remain pending.
[Prospective CPSP interface](../cpsp/README.md#prospective-carps-transcriptome-interface)

## Verify the current figures

The current verifier checks the render receipts and independently redraws all
45 published PNGs twice from their plot data. The inventory includes 13
comparison sets, five Step 02 figures and the palette preview. Each comparison
publishes learning curves, held-out metrics and same-fold, same-seed
differences from its reference. The metric panels show defined replicates,
means, sample SD when at least two values are defined, and zero.
[Current figure verification guide](AGENTS.md#verify-the-current-figures) ·
[Shared figure verifier](util/figure/script/verify_report_figures.py) ·
[Figure standard](../docs/adr/0109-cpst-seven-metric-comparison-figures.md)

## Historical evidence and active-input checks

The acceptance records document how the active-input verifier checked native
corpus windows, ordered projections, frozen folds and both tenv and Tahoe
readers in a checkout without the retired archive. The verification inventories
bind current formal releases to their tracked receipts. For Step 06.4
registrations with original `stage2/` source keys, follow the pinned replay
route.
[Active-input verification guide](AGENTS.md#verify-the-active-inputs) ·
[CarPS-Transcriptome acceptance](util/training/cpst-verification.md) ·
[TranSiGen acceptance](util/training/verification.md) ·
[Pinned Step 06.4 replay](06_tsg_module/06.4_round4/README.md#replay-a-sealed-step-064-fit)

The public boundary retires the local Stage 2 archive. The original tracked
code, figures and receipts are available at their pinned Git revisions.
[Archive retirement](../docs/adr/0105-stage-cleanup.md) ·
[Public release boundary](../docs/adr/0106-cpst-and-cpsp-public-release-boundaries.md) ·
[Original ADR 0087 cutover](https://github.com/Moke110/carps/blob/99078e6f/stage2/README.md#deploy-the-adr-0087-cutover) ·
[Original CarPS path migration](https://github.com/Moke110/carps/blob/ac62a177a54aa99e6ac626728af12394f3e0d940/stage2_CarPS-VAE/README.md#deploy-the-carps-stage-2-path)
