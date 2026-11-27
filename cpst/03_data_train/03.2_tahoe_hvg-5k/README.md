# 03.2 Tahoe HVG-DEG 5k well-level data

This release provides matched Tahoe expression on the ordered [HVG-DEG 5k axis](../../02_gene_axis/02.1_hvg-deg/README.md). The original build projected [native Tahoe counts](../../01_data_raw/README.md) and retained one row per eligible well and cell line, including vehicle profiles. It normalized each profile and its matched vehicle pool separately with CP10K and natural `log1p` on the selected axis. The matched pool sums integer counts from distinct same-plate, same-cell-line DMSO_TF observations. The `delta` layer is `X - control`.

## Accepted release

The release contains 66,718 profiles across 5,000 genes. Its [manifest](output/release/manifest.json) pins the axis and corpus generations, the original combined-view manifest, and the hashes of three local assets. `tahoe_bulk.h5ad` contains expression and matched controls, `tahoe_bulk.members.csv` identifies each profile's expression and matched vehicle members, and `tahoe_bulk.exclusions.csv` records excluded profiles and their reasons.

The [source-specific acceptance receipt](../log/312-source-specific-releases.json) records byte-exact replay of the original combined view, digest equality for all three Tahoe assets and a bounded read of the published release. The later [archive-isolated smoke test](../../04_tsg_axis/log/323-archive-isolated-smoke.json) verified the published Tahoe release with the local archive unavailable.

## Read or verify the release

Add this substep's `output/api` and `cpst/util/output/api` to `PYTHONPATH`. Use the [bounded reader](output/api/tahoe_expression.py) with `tahoe_expression.open_tahoe()`, or pass an explicit materialized release directory. The view exposes profile and gene metadata in `obs` and `var`. Use `read(start, stop, layer)` or `batches(layer, batch_size)` to retrieve up to 256 consecutive rows of `X`, `control` or `delta` per block.

Follow the [release verification procedure](AGENTS.md#read-and-verify-the-accepted-release) from a checkout with the release installed. If you have the original combined release, see the [reconstruction procedure](AGENTS.md#reconstruct-from-a-sealed-combined-release) and the [historical evidence](../../README.md#historical-evidence-and-active-input-checks).

The native corpus publisher uses this view to cross-check the gene order and endpoint profiles before replacing the accepted corpus. The active-input smoke test checks the Tahoe reader alongside the eight [tenv training releases](../README.md#tenv-training-axes). The [Step 04 TranSiGen axis comparison](../../04_tsg_axis/README.md) used six of those tenv views, and the [Step 09 CarPS-Transcriptome axis comparison](../../09_cpst_axis/README.md) used all eight.
