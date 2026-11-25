# Share Stage 2 plotting and update report figures in place

**Status:** accepted (2026-11-25).

Stage 2 util owns the shared chart primitives and palette. Study renderers
bind plotted CSV fields to marks, axes and labels and compose the report
panels through that API. A palette coordinate pairs one of nine HSL hues
(`20°` through `340°` in `40°` steps) with one of five depths from dark to
light. The approved lightness and saturation pairs are
`(0.24, 0.80)`, `(0.33, 0.76)`, `(0.42, 0.72)`, `(0.52, 0.66)`,
and `(0.62, 0.60)`. Both colors at a coordinate share the same RGB,
with opacity `0.90` for edges, dots, error bars and lines and `0.40`
for bar fills and range bands.

Use consistent hue–depth assignments across panels in the same study.
The [current Stage 2 plotting contract](../../stage2/AGENTS.md#step-contract)
supersedes the original two-variable convention for categorical marks.
Each hue–depth pair identifies one level of a single grouping variable.
Single-hue numeric heatmaps use depth for ordered values. Opacity
represents mark roles.

For the presentation-only refresh of the 22 active report figures, replace
the PNGs and renderers at their existing paths. This amends the successor
directory convention for this refresh. Verify reproduction from the
unchanged plotted CSVs and publish a new rendering receipt for the same
paths. Original acceptance records keep their evidence and point to a
pinned Git revision for replay. The local archive remains the historical
boundary defined in ADR 0087.

- [ADR 0062](0062-stage2-accepted-products-and-reproducible-figures.md)
- [ADR 0087](0087-stage2-active-steps-and-history.md)
