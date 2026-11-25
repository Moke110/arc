# Stage 2 perturbation prediction

## Selected gene axis

**Denoise-DEG 5k (`denoise-5k`)** is the Stage 2 prediction axis.
The [Step 04 matched-control comparison](04_tsg_axis/README.md#gene-axis-selection)
supports this choice based on held-out perturbation Delta prediction.
Across six axes and three cold-drug folds, Denoise-DEG 5k had the highest
mean Δ PCC (0.4066) and the lowest mean Δ RMSE (0.2596) and Δ MMD (0.1793).

Use the [ordered 5,000-gene list](02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv)
and its [Step 03.7 tenv training release](03_data_train/03.7_tenv_denoise-5k/README.md).
