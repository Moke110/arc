# TranSiGen architecture studies

[Step 06.1](06.1_round1/README.md) compares ten expression architectures
on the fixed Stage 2 gene axis and training objective. The encoder-width
extension adds six fits and reuses the 24 accepted fits. The original
architecture supplies the paired reference on each frozen fold.

[Step 06.2](06.2_round2/README.md#formal-results) compares direct Δ readouts,
drug–dose conditioning, encoder sharing and latent residuals. Seven new
architectures share the three frozen folds and the three Step 06.1
original fits supply the paired references.
