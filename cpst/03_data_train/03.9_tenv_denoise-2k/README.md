# 03.9 tenv training data on Denoise-DEG 2k

This release gives each pooled tenv treatment condition a measured DMSO
control from the same plate. The builder selects the ordered 2,000 Denoise-DEG
genes from the native integer counts before pooling wells. It normalizes each
pooled profile with CP10K on its selected-axis count total followed by natural
`log1p`, then subtracts the normalized control from the normalized treatment
to obtain Δ. [Step 02.2 axis](../../02_gene_axis/02.2_denoise-deg/README.md) · [Tenv file contract](../03.1_tenv_hvg-5k/README.md)

## Accepted release

The release contains 16 pooled plate controls and 1,383 pooled treatment
conditions. It excludes 112 of 4,510 source observations. Byte-exact replay
passed, and the ordered condition keys, matched control row IDs, raw
membership, exclusions, six fold assignments and strip reasons matched the
HVG-DEG 5k reference. The independent audit recomputed every pooled profile
from native counts and verified gene order, selected-axis library totals,
normalized expression, matched Δ and the dataset-wide control vector.
The manifest pins the gene axis, counts, frozen folds and all eight output
files. [Release manifest](output/manifest.json) · [2k acceptance receipt](../log/366-two-tenv-axes.json)

## Downstream use and reproduction

Step 09 evaluated this release on three cold-drug folds and three seeds in
the eight-axis CarPS-Transcriptome comparison with UniMol. The study selected
Denoise-DEG 5k. [Step 09 formal results](../../09_cpst_axis/README.md#formal-results)

With `cpst/03_data_train/output/api` on `PYTHONPATH`, call
`open_tenv_axis_training("denoise-2k")` from `tenv_axes`. Pass the Step 03 root
as the second argument to read outputs materialized in another checkout.
For publication and exact replay, use `--axis denoise-2k` and the stable release
path `cpst/03_data_train/03.9_tenv_denoise-2k/output`. Run the independent audit
after materializing both 2k releases.
[Named reader](../output/api/tenv_axes.py) · [Build and replay procedure](../README.md#tenv-training-axes) · [2k audit](../execution/audit_tenv_2k.py)
