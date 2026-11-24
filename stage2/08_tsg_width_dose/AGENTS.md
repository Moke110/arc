# TranSiGen width and dose encoding

Use the [historical evidence index](../03_data_train/03.1_hvg-deg-5k/historical/README.md) for sealed records.

Read [README.md](README.md), the [Stage 2 routing](../AGENTS.md),
[ADR 0062](../../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md),
and [ADR 0080](../../docs/adr/0080-tenv-kpgt-dose-round1.md).

Use `execution/raw_first.py prepare|evaluate|render` for the current
full raw-first release, verified tenv KPGT vectors and cold-drug-dose split.
Run `python -m unittest discover -s stage2/08_tsg_width_dose/execution/tests`
for fixed-prediction and CSV-only plotting verification.

Preserve the accepted figures, plotted values and study records. Run preparation
under the [execution environment limits](../../docs/agents/execution-environment.md).
