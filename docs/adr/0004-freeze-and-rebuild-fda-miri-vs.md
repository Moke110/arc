# fda_miri_vs frozen as xpert-fda-v1_history and rebuilt

> **Status: SUPERSEDED** by [ADR 0007](0007-finalize-xpert-fda-miri.md). The
> v1_history snapshot and the v1/v2 workflow variants are removed; the
> finalized workflow is `stage1/xpert-fda-miri-v1` (formerly xpert-fda-v3)
> and the sole frozen snapshot is `stage1/xpert-fda-miri-v0`.

The old `fda_miri_vs/` is renamed wholesale to a frozen provenance snapshot (`xpert-fda-v1_history`, its `history/` subdirectory deleted) and a fresh `xpert-fda-v1/` is built against the new taxonomy, with scoring logic refactored into numbered, dependency-ordered scripts. We chose freeze-and-rebuild over in-place refactoring so the new pipeline is cleanly structured while the old one remains a runnable reference. The snapshot is gitignored because it contains large artifacts that the repo never commits.
