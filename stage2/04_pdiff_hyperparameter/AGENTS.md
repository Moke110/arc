# AGENTS.md — Step 04 PerturbDiff training configuration

Read [README.md](README.md) and [ADR 0046](../../docs/adr/0046-perturbdiff-standalone-training.md).
Follow the repository's [execution limits](../../docs/agents/execution-environment.md).
Only the dedicated scheduling subagent may launch or restart Step 04 heavy work.

- `execution/script/schedule.py` owns resource admission, worker monitoring and recovery.
- `execution/script/formal.py` freezes the formal plan and runs separate fit and held-out prediction jobs.
- `execution/script/formal_outcomes.py` collects terminal outcomes and publishes the complete tenv checkpoint release.
- `execution/tests/test_schedule.py` exercises the scheduler CLI with bounded fixtures.
- `execution/scheduler_record/` preserves the bounded queue rehearsal.
- `execution/script/train.py` builds and verifies the Tahoe release, prepares frozen tenv and Tahoe requests, trains, and replays EMA checkpoints.
- `execution/tests/test_train.py` checks the monitor and predictor request through the CLI.
- `execution/smoke_record/` preserves the full-axis tenv fit and reload evidence.
- `execution/script/tahoe_inputs.py` owns the Tahoe projection, controls and cold-drug folds.
- `execution/tests/test_tahoe.py` checks the Tahoe release and predictor boundaries through the CLI.
- `execution/tahoe_smoke_record/` preserves the full-axis Tahoe fit and reload evidence.
- `execution/script/pilot.py` plans paired cadence profiles and freezes per-fold formal recipes.
- `execution/pilot_record/` preserves the pilot protocol and tenv outcome.
- `execution/script/evaluate.py` scores held-out predictions and exports per-run PR curves.
- `execution/script/publish.py` renders the ADR 0046 figures and verifies the compact formal records.
- `execution/tests/test_publish.py` checks replay and incomplete publication evidence through the CLI.
- `execution/formal_record/publication.json` binds the figures and the tenv checkpoint handoff.

- Active retrieval scoring and figures use [`../evaluation/README.md`](../evaluation/README.md), with `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` on the full CM significance axis.
