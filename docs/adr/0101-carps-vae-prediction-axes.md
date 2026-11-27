# Compare CarPS-VAE on eight gene axes with UniMol

**Status** Accepted design, 2026-11-27.

[ADR 0102](0102-stage2-axis-selection-and-comparison-figures.md) supersedes
the reporting scope and records the selected prediction axis.

[Step 09](../../cpst/09_cpst_axis/README.md) compares the published HVG-DEG and Denoise-DEG axes at 2k, 3k, 4k and 5k. Use the frozen Step 08 UniMol representation selected in [ADR 0100](0100-carps-vae-unimol.md) and the Step 03 plate-pooled tenv views, including new 2k views for both families. Each view normalizes on its own axis. The study estimates the performance of choosing an axis under this release policy, which includes changes to the normalized values of shared genes.

Adapt only the control and measured Δ encoder input sizes and the control and shared Δ decoder output sizes to the axis width. Pair unchanged initialized parameters and training random streams by cold-drug fold and training seed. Use folds 0–2 and seeds 17, 42 and 2026, with `denoise-5k` as the reference. Accept each Step 08 UniMol reference fit only after checking input and protocol identity, paired initialization and random-stream compatibility, and exact checkpoint replay. Fit any coordinate that fails the reuse gate.

Report metrics both on each full axis and on the axis intersected with the fixed union of the four eligible [disease DEG1000](0073-full-source-tier-first-deg1000-and-four-axes.md) lists. The intersection gives the genes with a predicted value. Directional coverage instead uses the full 3,481-gene union as its denominator so that it reflects both axis coverage and direction accuracy. Treat coverage as a held-out condition metric and report paired changes in percentage points. The disease lists define membership, and the measured perturbation Δ supplies the direction.
