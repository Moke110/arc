# 03.8 tenv training data on HVG-DEG 2k

This release pairs each tenv treatment condition with its measured same-plate
DMSO control on the ordered 2,000-gene HVG-DEG axis. The builder projects native
integer counts onto this axis before pooling wells. Each pooled profile uses
its 2k count total for CP10K followed by natural `log1p`. The released Δ
subtracts the normalized matched control from the normalized treatment.
[Step 02.1 axis](../../02_gene_axis/02.1_hvg-deg/README.md) · [Tenv file contract](../03.1_tenv_hvg-5k/README.md)

## Accepted release

The release contains 16 pooled plate controls and 1,383 pooled treatment
conditions, with 112 of 4,510 source observations excluded. The build replayed
byte for byte. Ordered condition keys, matched-control row IDs, raw membership,
exclusions, and all six frozen fold assignments and strip reasons matched the
HVG-DEG 5k reference. The independent raw-count audit verified all pooled
controls and treatments, including the ordered genes, selected-axis library
totals, normalized expression, matched Δ and dataset-wide control vector.
The manifest pins the prediction axis, native counts, frozen folds and all
eight output files.
[Release manifest](output/manifest.json) · [2k acceptance receipt](../log/366-two-tenv-axes.json)

## Downstream use and reproduction

Step 09 fitted this view across three cold-drug folds and three seeds in the
eight-axis CarPS-Transcriptome comparison with UniMol. The study selected
Denoise-DEG 5k. [Step 09 formal results](../../09_cpst_axis/README.md#formal-results)

With `cpst/03_data_train/output/api` on `PYTHONPATH`, call
`open_tenv_axis_training("hvg-2k")` from `tenv_axes`. Pass the Step 03 root
as the second argument to read a release materialized in another checkout.
Use `--axis hvg-2k` in the Step 03 publication and exact replay procedure,
with `cpst/03_data_train/03.8_tenv_hvg-2k/output` as the stable release path.
Run the independent audit after both 2k releases are materialized.
[Named reader](../output/api/tenv_axes.py) · [Build and replay procedure](../README.md#tenv-training-axes) · [2k audit](../execution/audit_tenv_2k.py)
