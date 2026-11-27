# CarPS-Transcriptome loss study

Step 07 compares additive MSE, L1, cosine and Pearson terms for the
shared signed Δ model. The study uses matched same-plate DMSO, Denoise-DEG
5k, frozen KPGT, training-fold standardized log dose and the shared five-fold
cold-drug protocol. The grid compares each term at four weights with a
single zero-added-loss reference. [Study design](../../docs/adr/0108-cpst-loss-and-shared-evaluation.md)

The historical alpha release selected the Step 06.4 `XΔ-XΔ` model with
predicted Δ L1 at weight 1.75. Use the original checkout and materialized
inputs to reproduce its nine accepted checkpoint replays.
[Historical replay guide](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/AGENTS.md)
