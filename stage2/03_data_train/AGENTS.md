# Stage 2 plate-matched training data

Read [README.md](README.md), ADR 0085 and the repository execution-environment
rules. This step owns the pooled training schema and bounded reader. Source
counts, the prediction axis and Step 07 study splits keep their own owners.

Build outside accepted output and verify deterministic replay before promotion.
Run real builds and replay serially under the repository RAM and thread caps.
Use `--raw`, `--axis` and `--fold-root` for linked worktrees.
