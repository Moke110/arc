# Keep Stage 2 history in one local archive

**Status:** accepted (2026-11-25).

Keep Steps 01–04 and Stage 2 shared utilities active. The active training
products are six tenv plate-pooled releases and the Tahoe well-level release.
The earlier paired plate-pooled release, superseded Step 01–04 material,
Step 05 and later, and retired Stage 2 rollback material go to the local,
Git-ignored `stage2/history`. It is the sole Stage 2 historical handoff.
Keep sealed bytes and record their original paths and replay provenance.

Keep live code and frozen inputs needed by the active workflows under active
owners so an unavailable archive does not block normal operation. Use
`hvg-{N}k` and `denoise-{N}k` across current interfaces. Retain the earlier
names inside sealed evidence and preserve the matched-control results and
the selection of `denoise-5k`.

## See

- [ADR 0062](0062-stage2-accepted-products-and-reproducible-figures.md) for accepted products and figure reproduction.
- [ADR 0084](0084-retire-stage2-bases-and-reserve-loss-step.md) for the former Step 05 reservation and downstream routing.
- [ADR 0085](0085-stage2-raw-and-training-data.md) for the earlier paired release and frozen folds.
- [ADR 0086](0086-tenv-axis-matched-plate-control.md) for the accepted axis-comparison protocol.
