# TranSiGen width and dose encoding

Read [README.md](README.md), the [Stage 2 routing](../AGENTS.md),
[ADR 0062](../../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md),
and [ADR 0080](../../docs/adr/0080-tenv-kpgt-dose-round1.md).

Reuse the verified Step 04 tenv KPGT freeze and cold-drug-dose split.
Run fits singly with the [execution environment](../../docs/agents/execution-environment.md)
limits and 10-second RAM and GPU sampling. Write fit artifacts beneath an
ignored candidate path. Verify the complete grid and reproduce all four
figures using only the figure CSVs before publishing the report.
