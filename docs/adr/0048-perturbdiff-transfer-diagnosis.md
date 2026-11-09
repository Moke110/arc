# Diagnose the Tahoe transfer gap with paired tenv controls

**Status: accepted (2026-11-09).**

[ADR 0047](0047-perturbdiff-tahoe-transfer-to-tenv.md) found a consistent
Transfer deficit against CM on the fixed tenv splits. Step 06 isolates
positive-dose initialization and target-stage base learning rate in separate
paired interventions. Both retain the author's validation diffusion loss
for EMA checkpoint selection. Screen on the cold-drug fold 0 validation
groups. Advance an intervention at +0.010 directional mean AP@100, check
direction across the other seeds, and replicate on a predeclared cold-drug
fold with its own Tahoe exclusion. Record the remaining CM gap and the
limits of descriptive source-domain evidence. The
[completed decision map](https://github.com/Moke110/arc/issues/84) records
the designs and interpretation rules.

`stage2/06_pdiff_transfer_optimize/` owns the diagnostics and their
artifacts. For this step, the requested locations for result figures
(`output/*.png`) and all newly generated model checkpoints
(`output/checkpoints/`) take precedence over the corresponding default
placements in [ADR 0030](0030-external-data-provenance-and-step-layout.md).
Checkpoint subdirectories and filenames use readable names and numeric
ordinals. ADR 0030's step structure and provenance rules otherwise apply.
