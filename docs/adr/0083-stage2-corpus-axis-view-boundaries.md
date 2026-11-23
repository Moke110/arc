# 0083. Separate Stage 2 corpus, axis selection and view publication

**Status:** accepted (2026-11-23), implementation pending.

## Context

Bases couples count preparation, prediction gene selection and expression views. The TSG view study owns replicate manipulations. The accepted ADR 0077 full view averages normalized vehicle profiles, while the ADR 0082 merge view aggregates raw vehicle counts. Extending the latter rule to every new view changes the control and Delta for affected profiles, so existing experiments require their original input generation.

## Decision

- `01_corpus` owns verified source counts at the finest published observation granularity, including tenv wells and existing Tahoe bulk profiles. Preserve native gene coverage and record matching status, including unmatched observations. Other cardiac perturbation count sources can be ingested with the same contract. Record exclusions for inputs without verified counts.
- `02_axis_view/02.1_hvg-deg` owns the existing HVG-DEG evidence, four ordered prediction gene axes and projection entry point. Prediction axes are subsets of the fixed 18,164-gene CM axis. Shared-axis participation requires complete source coverage of these genes.
- `03_view/03.1_hvg-deg_views` owns the published expression views and the existing view-study results. Apply the requested selection or pooling recipe to treated counts. Sum the raw counts of the matched vehicle union, counting each observation once. Normalize treated and control separately on the selected axis with CP10K and natural log1p. Record exclusions for paired views and preserve the pooled-label target semantics.
- `util` owns shared plotting and palette, gene and drug identities, drug representations, ECFP4, seeds and generic metric APIs. Each study owns its folds, metric aggregation and results.
- Renumber all downstream steps to 04–08 and migrate their active view inputs, metrics and plotting entry points. Preserve original inputs, manifests, recipes and results for historical playback. Verify new inputs, metric calls with fixed predictions and CSV-only rendering. Full downstream retraining is outside the migration.

## Consequences

Source observations can support new aggregation recipes. The shared axis has a single coverage requirement. New view releases carry new provenance, while historical findings remain tied to the input scale on which they were measured. The migration retains legacy replay support until the replacement entry points and path mapping are verified.

## See

- [ADR 0030](0030-external-data-provenance-and-step-layout.md)
- [ADR 0077](0077-five-k-only-axis-local-normalization.md)
- [ADR 0082](0082-view-variants-tenv-tahoe-replicate-structure.md)
- [Implementation spec](https://github.com/Moke110/arc/issues/286)

[Spec 2](https://github.com/Moke110/arc/issues/301) supersedes the active
variant publication in this decision. The [HVG-DEG 5k substep](../../stage2/03_view/03.1_hvg-deg-5k/README.md)
publishes full views for both sources. The previous generations and study
receipts remain under its historical boundary.
