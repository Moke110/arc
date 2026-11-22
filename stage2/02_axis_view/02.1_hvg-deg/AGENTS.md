# HVG-DEG axis method

Read [the README](README.md), [ADR 0073](../../../docs/adr/0073-full-source-tier-first-deg1000-and-four-axes.md) and [ADR 0083](../../../docs/adr/0083-stage2-corpus-axis-view-boundaries.md).
- Use the accepted tier-first DEG1000 and tenv-HVG evidence through the unchanged selection implementation in Bases. Preserve historical inputs and receipts byte-for-byte.
- Project through `prediction_axis.open_axes` and the native corpus reader. Require complete common-gene coverage before projecting counts.
- Publish a verified candidate outside `output/`, with no expression payloads.
- Run the unit tests and the measured replay commands in README. Verify a complete Tahoe and tenv corpus for release acceptance. Large payload hashing must follow the repository resource caps.
