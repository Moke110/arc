# Stage 2 plate-matched training data

Read [README.md](README.md), ADR 0085 and the repository execution-environment
rules. This step owns the pooled training schema and bounded reader. Source
counts, the prediction axis and Step 07 study splits keep their own owners.

For [Spec 5](https://github.com/Moke110/arc/issues/310), use
[03.1 tenv](03.1_tenv_hvg-deg-5k/AGENTS.md) for the source-specific plate-pooled
product and [03.2 Tahoe](03.2_tahoe_hvg-deg-5k/AGENTS.md) for the well-level
HVG-DEG 5k product.

Build outside accepted output and verify deterministic replay before promotion.
Run real builds and replay serially under the repository RAM and thread caps.
Use `--raw`, `--axis` and `--fold-root` for linked worktrees. Keep the ADR 0085
paired release and [archived full view](03.1_hvg-deg-5k/AGENTS.md) available for
historical readers and replay.
