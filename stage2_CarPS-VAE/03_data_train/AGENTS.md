# Stage 2 plate-matched training data

Read [README.md](README.md), ADR 0085 and the repository execution-environment
rules. This step owns the pooled training schema, bounded readers and the frozen
fold inputs under `output/frozen-folds/`. Source counts and the prediction
axis keep their own owners.

For [Spec 5](https://github.com/Moke110/carps/issues/310), use the
[six named tenv releases](README.md#six-tenv-training-axes) for plate-pooled
data and [03.2 Tahoe](03.2_tahoe_hvg-5k/AGENTS.md) for well-level HVG-DEG 5k.

Build outside accepted output and verify deterministic replay before promotion.
Run real builds and replay serially under the repository RAM and thread caps.
The tenv builder selects a prediction axis with `--axis`. Use `--raw`,
`--axis-root`, `--fold-root` and `--reference` as needed for linked worktrees.
