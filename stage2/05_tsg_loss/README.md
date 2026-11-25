# TranSiGen loss studies

[Step 05.1](05.1_loss_metrics/README.md) compares MSE, cosine, Pearson correlation
and L1 Delta losses with coefficients 0.75, 1 and 1.25 on the denoise-5k axis.

[Step 05.2](05.2_mse_l1/README.md) varies the weight of an additive Delta
L1 term while retaining the native TranSiGen objective.

[Step 05.3](05.3_mse-l1_mse-x/README.md) varies the shared expression
reconstruction MSE weight at fixed Delta L1 weight 1.75.
