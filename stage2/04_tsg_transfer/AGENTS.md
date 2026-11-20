# TranSiGen Tahoe-to-tenv transfer

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md), ADR 0062,
and [ADR 0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md)
before changing or running this step.

Use the verified Step 03 axis-local freeze and immutable native TranSiGen
source. Keep the Tahoe split and checkpoint choice independent of tenv.
Use `execution/tahoe.py` to prepare folds, fit each fold, replay the best
checkpoints in fresh processes and lock the winner. Run formal fits singly
with the capped execution and 10-second resource monitoring required by
the [execution environment](../../docs/agents/execution-environment.md).

Keep candidates under `execution/checkpoints/` or an external staging
directory. Publish accepted outputs only after scientific, replay and
repository checks pass. Step 03 and Stage 3 retain their existing contracts.
