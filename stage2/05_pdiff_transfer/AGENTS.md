# AGENTS.md — Step 05 PerturbDiff transfer

Read [README.md](README.md), [ADR 0047](../../docs/adr/0047-perturbdiff-tahoe-transfer-to-tenv.md)
and the repository [execution limits](../../docs/agents/execution-environment.md).
Only the dedicated scheduling subagent admits or restarts heavy work.

- `execution/script/cold_tenv.py` verifies the shared frozen tenv fold, reuses
  the Step 04 CM predictions and scores both deterministic references once.
- `execution/tests/test_cold_tenv.py` checks the cold comparison CLI.
- `execution/script/cold_tahoe.py` builds and verifies the globally cold
  Tahoe fold 0 source release, then runs the bounded fit and EMA replay.
- `execution/tests/test_cold_tahoe.py` checks the source exclusion rule.
