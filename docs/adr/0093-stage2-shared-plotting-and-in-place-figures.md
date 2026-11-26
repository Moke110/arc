# Share Stage 2 plotting and update report figures in place

**Status:** accepted (2026-11-25).

Stage 2 util owns the shared chart primitives and palette. Study renderers
bind plotted CSV fields to marks, axes and labels and compose the report
panels through that API. A palette coordinate pairs one of nine HSL hues
(`20°` through `340°` in `40°` steps) with one of four depths from dark to light. The lightness and saturation
pairs are `(0.33, 0.76)`, `(0.42, 0.72)`, `(0.52, 0.66)` and
`(0.62, 0.60)`. Both colors at a coordinate share the same RGB, with
opacity `0.90` for edges, dots, error bars and lines and `0.40` for bar
fills and range bands.

The [Stage 2 plotting contract](../../stage2_CarPS-VAE/AGENTS.md#step-contract)
defines chart defaults, replicate presentation and canonical publication.
Map the first grouping variable to hue and an optional second variable to
depth, using a consistent mapping across related panels. For learning
curves use one panel per model, with a train and validation line spanning
the observed epochs of every fold and seed. Each score group shows the
mean, the defined measurements and sample standard deviation. A singleton
has one dot and its mean bar.
Focus score y-axes near the plotted data when it helps compare groups.
Keep all finite scores, complete SD intervals and meaningful reference
values in range. Label bar figures whose y-axes omit zero.

Each publication keeps one current set of PNGs, plotting data and renderers
under `figure/`. Re-render from the complete plotted data and refresh the
rendering receipt in place. Scientific receipts continue to record the
analysis evidence from their publication.

- [ADR 0062](0062-stage2-accepted-products-and-reproducible-figures.md)
- [ADR 0087](0087-stage2-active-steps-and-history.md)
