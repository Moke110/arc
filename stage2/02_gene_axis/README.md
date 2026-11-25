# Stage 2 prediction axes

[HVG-DEG](02.1_hvg-deg/README.md) publishes four nested prediction axes on the fixed 18,164-gene common axis. Projection reads integer counts from the [Stage 2 raw data release](../01_data_raw/README.md) in prediction-axis order.

[Denoise-DEG](02.2_denoise-deg/README.md) publishes five nested axes from 1,000 to 5,000 genes using tenv denoising, tenv top-DEG presence and eligible disease DEG1000 membership.

Use the same prediction axis ID for selection, tenv training and Step 04.
HVG-DEG publishes `hvg-2k` through `hvg-5k`. Denoise-DEG publishes
`denoise-1k` through `denoise-5k`. With
`stage2/02_gene_axis/02.1_hvg-deg/output/api` on `PYTHONPATH`,
`prediction_axis.open_axis("denoise-5k")` verifies and returns the ordered genes.
