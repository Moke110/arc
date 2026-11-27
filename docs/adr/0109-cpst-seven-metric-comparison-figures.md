# Retire wJAC-2000 and standardize CPST comparison panels

**Status:** accepted (2026-11-28).

Remove wJAC-2000 across active CPST evaluation and reporting before
freezing the new folds, moving the steps or running new fits. Remove its
calculation and API output, score tables and summaries, plots, and the
corresponding documentation and checks. Keep the historical evidence
available through the sealed releases and pinned Git revisions.

Use a two-row, four-column grid with seven panels for both `metrics.png`
and `metrics-paired.png`. Apply the grid to the shared comparison renderer
and the Step 10 renderer.

| Column 1 | Column 2 | Column 3 | Column 4 |
| --- | --- | --- | --- |
| ΔSpearman | JAC-top200 | ΔPCC | ΔRMSE |
| ΔR² | ΔE-distance | ΔMMD | |

Regenerate the active comparison score tables and metric figures from
verified existing results for this first item. Verify both layouts and
two CSV-only renders of each published metric figure. This decision amends
the implementation sequence in ADR 0108 and supersedes ADR 0102's
wJAC-2000 and panel layout requirements.

## See

- [ADR 0108](0108-cpst-loss-and-shared-evaluation.md)
- [ADR 0102](0102-stage2-axis-selection-and-comparison-figures.md)
