# TranSiGen Tahoe-to-tenv transfer

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md),
[ADR 0062](../../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md)
and [ADR 0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md).

Prepare new studies with `execution/raw_first.py` and the accepted full view.
Keep the Tahoe split and checkpoint choice independent of tenv. For the sealed
report, use the original Step 04 axis-local freeze and immutable native
TranSiGen source. Follow the [historical commands](README.md#historical-reproduction).
Run formal phases singly under
`execution/capped_phase.sh` with per-phase resource records and apply the
[execution environment](../../docs/agents/execution-environment.md) limits.

Work in a staging candidate. Use `execution/publish.py prepare` and a
fresh-process `verify` before installing the accepted report. Render the
figure independently from its CSVs and run the repository health gates.
Keep the formal artifacts and raw resource evidence under the ignored
`execution/checkpoints/accepted/` recovery path. Leave the sealed Step 04
and Stage 3 products intact.
