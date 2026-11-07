# AGENTS.md — Step 04 PerturbDiff training configuration

Read [README.md](README.md) and [ADR 0046](../../docs/adr/0046-perturbdiff-standalone-training.md).
Follow the repository's [execution limits](../../docs/agents/execution-environment.md).
Only the dedicated scheduling subagent may launch or restart Step 04 heavy work.

- `execution/script/schedule.py` owns resource admission, worker monitoring and recovery.
- `execution/tests/test_schedule.py` exercises the scheduler CLI with bounded fixtures.
- `execution/scheduler_record/` preserves the bounded queue rehearsal.
- `execution/script/train.py` prepares frozen tenv requests, trains, and replays EMA checkpoints.
- `execution/tests/test_train.py` checks the monitor and predictor request through the CLI.
- `execution/smoke_record/` preserves the full-axis tenv fit and reload evidence.
