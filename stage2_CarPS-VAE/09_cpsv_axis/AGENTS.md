# CarPS-VAE prediction axis comparison

Read the [study design](README.md), [Stage 2 routing](../AGENTS.md), [axis catalog](../02_gene_axis/README.md), [training releases](../03_data_train/README.md), [Step 08 reference](../08_cpsv_chem-emb/README.md) and [UniMol decision](../../docs/adr/0100-carps-vae-unimol.md).

Follow [ADR 0102](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md), the fit protocol in [ADR 0101](../../docs/adr/0101-carps-vae-prediction-axes.md), and the [execution limits](../../docs/agents/execution-environment.md).
Build and verify the two 2k views under Step 03 ownership. Verify the full roster with a bounded CUDA gate before fitting the formal matrix. Audit each reference reuse, pair only by fold and seed, and verify the three full-axis figures independently from their plotted CSVs.

Use the [formal reproduction commands](README.md#reproduce-and-verify) to fit in a capped scope and verify a staged comparison before installation. Keep the coordinator sampler live through the final matrix sample. The formal bundle is materialized locally at `output/formal`. The tracked acceptance receipt and figure files are published with the study code.
