# Step 09 CarPS-VAE prediction axes

**Status** Accepted design. See [ADR 0101](../../docs/adr/0101-carps-vae-prediction-axes.md).

Compare CarPS-VAE across eight [Step 02 prediction axes](../02_gene_axis/README.md), using the [Step 08 UniMol representation](../../docs/adr/0100-carps-vae-unimol.md) and pairing each arm with `denoise-5k` by cold-drug fold and training seed. Use folds 0–2 and seeds 17, 42 and 2026, for 72 axis × fold × seed coordinates.

| Axis family | Axis IDs | Disease-union overlap at 2k / 3k / 4k / 5k |
| --- | --- | --- |
| HVG-DEG | `hvg-2k`, `hvg-3k`, `hvg-4k`, `hvg-5k` | 1,258 / 1,889 / 2,525 / 3,149 |
| Denoise-DEG | `denoise-2k`, `denoise-3k`, `denoise-4k`, `denoise-5k` | 906 / 1,225 / 1,375 / 1,513 |

[Step 03](../03_data_train/README.md) owns the two additional 2k training views and the six established views. Pool matched-plate counts and normalize on each selected axis with the same method. Verify common condition keys, raw membership, exclusions and frozen partition order. Shared genes can have different normalized values across axes.

Resize the input projections of both expression encoders and the output projections of the control and shared Δ decoders. Keep the internal widths and architecture, UniMol input and projection, matched-control protocol, fold-local dose transform, loss and fit policy fixed. Use common initialized values for same-shaped parameters and a paired training random stream. Anchor pairing to the Step 08 UniMol reference so equivalent reference fits can be reused. Verify each reference coordinate's provenance, protocol equivalence and exact checkpoint replay. With all nine reference coordinates accepted, 63 coordinates require new training. Freshly fit any rejected reference coordinate.

## Figures

Use the nine Stage 2 metrics in their [Step 04 report order](../04_tsg_axis/README.md#score-and-replay-an-axis-comparison): JAC_top200, JAC_pos100, JAC_neg100, Δ PCC, Δ Spearman, Δ RMSE, Δ E-distance, Δ MMD and Δ R².

Let `D` be the 3,481 distinct frozen gene indices in the eligible DEG1000 union for nHF, ICM, genetic cardiomyopathy and TOX, and let `A` be the prediction axis for a run. Define `Δpred_d` and `Δ_d` by slicing predicted and measured Δ to `A ∩ D` in axis order. The disease references define membership. The two response vectors provide the directions and magnitudes that are compared.

| Figure | Plotted values |
| --- | --- |
| `learning-curves` | Training and validation objectives by epoch, with one panel per axis and both traces for each fold and seed. |
| `axis-metrics` | Nine held-out metrics on each full prediction axis. |
| `axis-metrics-paired` | Each full-axis score minus the matching `denoise-5k` score. |
| `disease-metrics` | Nine held-out metrics between `Δpred_d` and `Δ_d`. |
| `disease-metrics-paired` | Each disease-restricted score minus the matching `denoise-5k` score. |
| `disease-correct-coverage` | Per held-out condition, `100 × G / 3481`, where `G` counts genes with both values positive or both values negative in the disease-restricted vectors. A zero in either vector contributes zero to `G`. |
| `disease-correct-coverage-paired` | Coverage minus the matching `denoise-5k` coverage, in percentage points. |

Use the Step 04 condition-level metrics, eligibility rules and foldwise distribution metrics for both scopes. Use the same full-axis MMD bandwidth as Step 04 for each fold. For the disease-restricted scope, use the same bandwidth estimator on the HVG-DEG 5k measured training Δ sliced to its disease-union intersection. Freeze one bandwidth per fold and scope across axes and seeds. Record the bandwidth, training rows and gene indices used.

Average defined condition-level scores and coverage within each coordinate. Compute distribution scores from all held-out rows of that coordinate. Summarize defined coordinates with equal weight and sample SD. Subtract the reference only for matching fold and seed coordinates with two defined scores. Retain eligible counts and undefined values. Display axis overlap counts with the disease results. The `denoise-5k` direction-correct coverage ceiling is 1,513 / 3,481, or 43.46%.

Publish all seven figures, their replicate values and plotted CSVs, and CSV-only renderers under the [Stage 2 figure contract](../AGENTS.md#step-contract). Use one dot per defined coordinate, means and sample SDs. Include zero on the paired plots and the Δ R² panels.
