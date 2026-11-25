# Compare TranSiGen Delta losses and weights on tenv

**Status:** accepted (2026-11-25).

Reopen Step 05 (`05_tsg_loss`) as an active study. Preserve the archived
history. Use the Step 04 matched-plate protocol with the selected `denoise-5k`
axis, frozen cold-drug folds 0–2 and seed 42.

Hold the TranSiGen model and inputs fixed while crossing MSE, Cos, PCC and L1
predicted Delta losses with alpha values 1/4, 1/3 and 1/2. Fix the sum of the
three expression-term coefficients at 3. Assign the Delta term coefficient
`3 * alpha` and each native reconstruction MSE coefficient
`1.5 * (1 - alpha)`. Keep the summed native KL terms at coefficient 0.1.
The MSE coordinate at alpha 1/3 uses the author implementation directly.
This anchor distinguishes changing the relative weight of the native
prediction term from changing its functional form.

MSE and L1 sum across genes and conditions. Cos and PCC sum one minus the
per-condition similarity and multiply by axis width to align the aggregation
unit with the native reconstruction sums. Compute PCC as cosine after
centering each vector across genes. Form predicted and measured Delta
against the same measured plate control. Each coordinate trains and selects
its checkpoint on its own objective.

Use the same nine held-out Delta metrics as Step 04. Compare each metric
directly across the grid with fold-level evidence. Publish grouped bars and
fold-specific learning curves from CSVs.

## See

- [ADR 0062](0062-stage2-accepted-products-and-reproducible-figures.md)
- [ADR 0086](0086-tenv-axis-matched-plate-control.md)
- [ADR 0087](0087-stage2-active-steps-and-history.md)
- [TranSiGen paper](https://www.nature.com/articles/s41467-024-49620-3)
- [TranSiGen v1.0 loss](https://github.com/myzhengSIMM/TranSiGen/blob/v1.0/src/model.py#L178-L190)
