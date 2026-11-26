# AGENTS.md — Stage 2 shared utilities

Read the [public contract](README.md), [Stage 2 routing](../AGENTS.md)
and [ADR 0083](../../docs/adr/0083-stage2-corpus-axis-view-boundaries.md).

- Import `stage2_util` from `output/api` and reuse compatible public operations.
- Keep folds, metric aggregation, and formatted result tables in each study.
- Open accepted parent and representation data through `RepresentationTable`.
- Preserve sealed legacy files when updating current defaults.
- Run the test command in README before publishing an API change.
