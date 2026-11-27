# Step 09 CarPS-Transcriptome prediction axes

Step 09 compared HVG-DEG and Denoise-DEG at 2k, 3k, 4k and 5k genes to select the CarPS-Transcriptome prediction axis with frozen UniMol. The study used the eight plate-pooled tenv training views, including both new 2k views. The views shared condition keys, raw membership, exclusions and frozen cold-drug partitions. Each view normalized counts on its own axis, so a gene shared by two axes can have different normalized values. [Prediction axes](../02_gene_axis/README.md) · [Training views](../03_data_train/README.md) · [Study design](../../docs/adr/0101-carps-vae-prediction-axes.md)

Only the input projections of the two expression encoders and the output projections of the control and shared Δ decoders changed size with the axis. All arms used the same internal architecture, UniMol representation and projection, matched-plate DMSO controls, training-fold standardized log dose, objective and fit policy. The study paired same-shaped parameter initialization and training random streams by fold and seed, with the Step 08 UniMol Denoise-DEG 5k fits as the reference. [Step 08 comparison](../08_cpst_chem-emb/README.md#formal-results) · [UniMol decision](../../docs/adr/0100-carps-vae-unimol.md)

## Bounded CUDA verification

All eight axes passed the bounded CUDA gate on cold-drug fold 0 with seed 17. The receipt records independent checkpoint replays, initialization and per-epoch random-stream hashes, scores and figure replays for the original study. The preflight registered all 72 formal coordinates, including the nine Step 08 UniMol reference candidates. The gate also verified all nine Step 07 alpha replays. [Bounded receipt](log/367-bounded.json) · [Preflight](log/367-preflight.json) · [Step 07 replay](log/367-step07-replay.json)

## Formal results

**Step 09 selected Denoise-DEG 5k with UniMol.** It had the highest mean held-out full-axis Δ PCC (0.4413 ± 0.0136) and lowest mean Δ RMSE (0.2386 ± 0.0037). Denoise-DEG 4k was closest in Δ PCC. Its paired differences from Denoise-DEG 5k were −0.0017 ± 0.0083 in Δ PCC and +0.0249 ± 0.0012 in Δ RMSE. The decision applies to the released axis-local normalization policy. Step 10 held the selected ordered axis fixed for its dose comparison. [Selection in ADR 0102](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md) · [Step 10 dose comparison](../10_cpst_dose/README.md#formal-results)

The formal matrix covered eight axes, folds 0–2 and seeds 17, 42 and 2026. All 63 new fits and nine Step 08 reference fits exactly replayed their validation and held-out predictions. The full-axis audit verified 720 score rows and 33,192 held-out condition rows across all 72 coordinates. Table entries give the equal-weight mean ± sample SD across the nine fold and seed coordinates per axis. Paired differences subtract the Denoise-DEG 5k score for the same fold and seed where both scores are defined. [Formal acceptance receipt](log/371-formal.json) · [Reuse audit](log/368-reference.json) · [Fit lineage](log/runs.csv) · [Full-axis summaries](log/summary.csv)

| Axis | Full-axis Δ PCC ↑ | Full-axis Δ RMSE ↓ |
| --- | ---: | ---: |
| HVG-DEG 2k | 0.4209 ± 0.0123 | 0.4712 ± 0.0050 |
| HVG-DEG 3k | 0.4160 ± 0.0088 | 0.3371 ± 0.0047 |
| HVG-DEG 4k | 0.4243 ± 0.0110 | 0.3000 ± 0.0051 |
| HVG-DEG 5k | 0.4247 ± 0.0099 | 0.2764 ± 0.0043 |
| Denoise-DEG 2k | 0.4212 ± 0.0183 | 0.3210 ± 0.0050 |
| Denoise-DEG 3k | 0.4310 ± 0.0138 | 0.2865 ± 0.0039 |
| Denoise-DEG 4k | 0.4396 ± 0.0146 | 0.2635 ± 0.0043 |
| Denoise-DEG 5k | **0.4413 ± 0.0136** | **0.2386 ± 0.0037** |

The active report uses full-axis scoring under ADR 0102. The original disease-overlap results, figures and receipts remain available in the fixed Git revision. [Original dual-scope publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/09_cpsv_axis/README.md#formal-results)

## Figures

The scoring follows the Step 04 condition eligibility and foldwise distribution metrics. Step 09 verified the Step 04 full-axis MMD bandwidth for each fold and held it fixed across axes and seeds. The metric figures plot seven held-out metrics with replicate dots, means, sample SD and zero visible in every panel. Two independent renders reproduced the three figures from the plotted CSVs. [Step 04 comparison](../04_tsg_axis/README.md) · [MMD calibration](log/bandwidths.csv) · [Figure receipt](log/figure-render.json)

| Figure | Plot data |
| --- | --- |
| [Learning curves](figure/learning-curve.png) | [Training and validation objectives by epoch, including reused fits](figure/data/learning-curve.csv) |
| [Full-axis metrics](figure/metrics.png) | [Scores by axis, fold and seed](figure/data/metrics.csv) |
| [Paired full-axis metrics](figure/metrics-paired.png) | [Differences from Denoise-DEG 5k by axis, fold and seed](figure/data/metrics-paired.csv) |

## Reproduce and verify

The active figures and score tables retain the accepted values for the
current metrics. Use the pinned publication for full scientific replay.
[Sealed publication](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/09_cpst_axis/README.md)

The installed `output/formal` bundle holds the sealed fit matrix, reuse plan and audited report. Their manifests are bound by the [formal acceptance receipt](log/371-formal.json). The owning guide has the procedures to [run the bounded CUDA gate](AGENTS.md#run-the-bounded-cuda-gate), [reproduce and publish the formal study](AGENTS.md#reproduce-and-publish-the-formal-study), [verify the installed release](AGENTS.md#verify-the-installed-release), and [rebuild the report from sealed fits](AGENTS.md#rebuild-the-report-from-sealed-fits).
