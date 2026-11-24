# Stage 2 plate-matched training data

Read [README.md](README.md), ADR 0085 and the repository execution-environment
rules. This step owns the pooled training schema and bounded reader. Source
counts, the prediction axis and Step 07 study splits keep their own owners.

For [Spec 5](https://github.com/Moke110/arc/issues/310), use the
[six named tenv releases](README.md#six-tenv-training-axes) for plate-pooled
data and [03.2 Tahoe](03.2_tahoe_hvg-deg-5k/AGENTS.md) for well-level HVG-DEG 5k.

Build outside accepted output and verify deterministic replay before promotion.
Run real builds and replay serially under the repository RAM and thread caps.
The tenv builder selects axes with `--family` and `--size`. Use `--raw`,
`--axis`, `--fold-root` and `--reference` as needed for linked worktrees. Keep
the ADR 0085 paired release and [archived full view](03.1_hvg-deg-5k/AGENTS.md)
available for historical readers and replay.
