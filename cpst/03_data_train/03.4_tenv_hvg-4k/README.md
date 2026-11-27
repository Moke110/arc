# 03.4 tenv training data on HVG-DEG 4k

This release lets axis comparisons evaluate the ordered 4,000-gene HVG-DEG
axis on the same tenv conditions and folds as the 5k reference. The builder
selects the 4,000 genes from native integer counts before pooling wells.
It applies CP10K and natural `log1p` using each pooled profile's 4k count
total. Treatment Δ subtracts the normalized DMSO control measured on the
same plate. Genes shared with the 5k axis can therefore have different
normalized expression values. [HVG-DEG axis](../../02_gene_axis/02.1_hvg-deg/README.md) · [tenv release contract](../03.1_tenv_hvg-5k/README.md)

The accepted release contains 16 plate controls and 1,383 treatment
conditions, with 112 exclusions among 4,510 raw observations. The six-axis
acceptance receipt records byte-exact replay for 4k and equality with the 5k
reference for the ordered condition keys, control matches, membership and
exclusion ledgers, and all six frozen fold assignments and strip reasons.
The manifest pins the 4k axis, native counts, fold provenance and output
digests. [Release manifest](output/manifest.json) · [Acceptance receipt](../log/313-six-tenv-axes.json) · [Frozen fold inputs](../output/frozen-folds/tenv/manifest.json)

HVG-DEG 4k was evaluated in Step 04 with TranSiGen and in Step 09 with
CarPS-Transcriptome and UniMol. Both comparisons selected Denoise-DEG 5k.
[Step 04 selection](../../04_tsg_axis/README.md#gene-axis-selection) · [Step 09 formal results](../../09_cpst_axis/README.md#formal-results)

With `cpst/03_data_train/output/api` on `PYTHONPATH`, use
`open_tenv_axis_training("hvg-4k")` from `tenv_axes` to open the verified
release. Follow the Step 03 candidate and replay procedure with
`--axis hvg-4k` and use `cpst/03_data_train/03.4_tenv_hvg-4k/output` as the
stable release path. [Reader implementation](../output/api/tenv_axes.py) · [Build and replay procedure](../README.md#tenv-training-axes) · [Builder](../execution/build_tenv_training.py)
