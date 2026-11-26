# TranSiGen Delta PCC sweep

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md), and
[ADR 0092](../../../docs/adr/0092-transigen-bounded-delta-similarity-sweeps.md).

Use the shared [similarity runner](../similarity/run_similarity.py) with
`--metric pcc`. Follow the [capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
and the [execution environment](../../../docs/agents/execution-environment.md).
Keep one CUDA fit active and log WSL RAM and GPU usage every ten seconds.

Before accepting a release, verify all 18 independently replayed runs,
the three fresh zero-weight references against Step 05.3, the frozen
bandwidths and the byte-identical CSV-only figures.
