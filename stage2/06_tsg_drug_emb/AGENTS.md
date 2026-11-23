# TranSiGen drug representation comparison

Read [the report](README.md), [Stage 2 routing](../AGENTS.md),
[ADR 0079](../../docs/adr/0079-stage2-drug-representation-step.md),
and the [execution environment](../../docs/agents/execution-environment.md).
Reuse the accepted source-local cold-drug folds and dose conditioning from the
[dose comparison](../07_tsg_dose/AGENTS.md).

## Raw-first comparison

Use `execution/raw_first.py` for new comparison preparation, fixed-prediction
evaluation and CSV-only plotting. The input is the source directory emitted by
the [dose step's raw-first preparation](../07_tsg_dose/README.md#raw-first-inputs-evaluation-and-figures).
Follow the [active commands](README.md#raw-first-inputs-evaluation-and-figures).
The accepted source-specific feature tables are owned by [util](../util/README.md).
Retain the [study records](log/) and accepted figures. The report renderer
reads plotted values from `figure/data/`.
