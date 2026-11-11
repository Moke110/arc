# AGENTS.md — Step 05 PerturbDiff transfer

Read [README.md](README.md), [ADR 0047](../../docs/adr/0047-perturbdiff-tahoe-transfer-to-tenv.md)
and the repository [execution limits](../../docs/agents/execution-environment.md).
Only the dedicated scheduling subagent admits or restarts heavy work.

- `execution/script/cold_tenv.py` verifies and scores shared cold tenv inputs.
  `execution/tests/test_cold_tenv.py` checks the comparison CLI.
- `execution/script/cold_tahoe.py` builds and verifies the globally cold
  Tahoe fold 0 source release, then runs the bounded fit and EMA replay.
- `execution/tests/test_cold_tahoe.py` checks the source exclusion rule.
- `execution/script/warm_tenv.py` freezes the warm-start tenv split and Step 04
  recipes, then verifies, fits, predicts and scores through separate commands.
- `execution/tests/test_warm_tenv.py` exercises the warm-start artifact CLI.
- `execution/script/mix.py` freezes the nine directly pooled Mix slots and
  verifies, fits, predicts and scores each regime. `mix_inputs.py` dispatches
  native minibatches to each source's matched-control inputs.
- `execution/tests/test_mix.py` checks the pooled predictor and batch inputs.
- `execution/script/tahoe_zero_shot.py` verifies the historical Tahoe fold 0
  states and predicts and scores their zero-shot tenv responses.
- `execution/script/tahoe_bridge.py` seals and loads the shared initialization
  for Tahoe and Transfer. `execution/tests/test_tahoe_bridge.py` covers the
  shared loader, remapped dose rows and rebuilt recipient lookups.
- `execution/script/transfer.py` fits globally cold sources and Transfer
  recipients using the Tahoe prediction bridge, then replays both arms.
- `execution/script/formal_matrix.py` freezes and schedules the formal matrix.
  `execution/tests/test_formal_matrix.py` checks orchestration and seals.
- `execution/script/publish.py` verifies and scores the terminal inventory and
  renders the compact results. `execution/tests/test_publish.py` checks paired
  scoring coverage and offline figure replay on an incomplete fixture.

- Active retrieval scoring and figures use [`../evaluation/README.md`](../evaluation/README.md), with `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` on the full CM significance axis.
