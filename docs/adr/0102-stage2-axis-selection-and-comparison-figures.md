# Select Denoise-DEG 5k and standardize Stage 2 comparison figures

**Status** Accepted design, 2026-11-27.

Select the ordered Denoise-DEG 5k (`denoise-5k`) prediction axis for
CarPS-VAE with UniMol. In the Step 09 comparison, it had the highest
mean held-out full-axis Δ PCC (0.4413 ± 0.0136) and lowest mean Δ RMSE
(0.2386 ± 0.0037) over nine cold-drug fold and seed coordinates. The
closest axis on Δ PCC, Denoise-DEG 4k, differed by −0.0017 ± 0.0083
on Δ PCC and +0.0249 ± 0.0012 on Δ RMSE in paired comparisons. This
selects the prediction axis under the released axis-local normalization
policy. See the [original Step 09 results](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/09_cpsv_axis/README.md#formal-results).

This decision supersedes the dual-scope reporting requirement of
[ADR 0101](0101-carps-vae-prediction-axes.md). Step 09 reports only
full-axis scores and removes its active disease scoring, data, figures
and verification branches. The original disease results and their
receipts remain available in the fixed published Git revision above.
Step 02 owns the shared evidence used to construct the prediction axes.

Publish `learning-curve.png`, `metrics.png` and `metrics-paired.png`
for each comparison in Steps 04, 05.1–05.5, 06.1–06.4, 08 and 09.
Step 07 links the Step 06.4 comparison evidence for its fixed alpha.
Retire the extra Step 04 disease F1 and Step 06.3/06.4 contrast PNGs
from the active figures, with original results available at a fixed
revision. Both metrics figures use these rows in order.

| First | Second | Third |
| --- | --- | --- |
| ΔSpearman | wJAC-2000 | JAC-top200 |
| ΔPCC | ΔRMSE | ΔR² |
| ΔE-distance | ΔMMD | |

Keep zero in every metrics panel and retain the replicates, means and
sample SD convention of [ADR 0093](0093-stage2-shared-plotting-and-in-place-figures.md).
Pair within fold and seed and plot experiment minus reference in each
metric's units. Step 04 uses `denoise-5k`; existing reference arms
are used for the other comparisons.

For each held-out condition, apply the existing signed Jaccard definition
at each cutoff K, taking up to K/2 strictly positive and K/2 strictly
negative genes. Break ties within each direction by frozen gene index and
sum the directional intersections and unions before taking their ratio.
An empty combined union gives an undefined score. Compute wJAC-2000 as
`(0.125 × JAC-top200 + 0.25 × JAC-top500 + 0.5 × JAC-top1000 + JAC-top2000) / 1.875`
before aggregating along the study's existing condition, fold and seed
boundaries. Source the added scores and curves from verified fit artifacts.
