# 03.5 tenv training data on Denoise-DEG 3k

This release provides plate-matched tenv profiles for training on the ordered
3,000-gene Denoise-DEG axis. It selects those genes from native integer counts
before pooling DMSO wells by plate and treated wells by plate, canonical drug
parent and measured dose. Each pooled profile uses its own 3k count total for
CP10K and natural `log1p`. The released Δ subtracts the normalized same-plate
DMSO profile from the normalized treatment.
[Denoise-DEG axes](../../02_gene_axis/02.2_denoise-deg/README.md) ·
[Tenv release contract](../03.1_tenv_hvg-5k/README.md)

## Accepted release

Of 4,510 source observations, 112 were excluded. The release contains 16
plate controls and 1,383 plate–drug–dose treatment conditions. The acceptance
receipt records byte-exact replay and equality with the HVG-DEG 5k reference
for the ordered condition keys, matched controls, raw membership, exclusions,
and all six frozen fold assignments and strip reasons. The manifest pins the
gene order, raw counts, frozen folds and output file digests.
[Release manifest](output/manifest.json) ·
[Six-axis acceptance receipt](../log/313-six-tenv-axes.json)

## Downstream use and access

Step 04 evaluated Denoise-DEG 3k in the six-axis TranSiGen comparison. Step 09
evaluated it in the eight-axis CarPS-Transcriptome comparison with UniMol and
selected Denoise-DEG 5k.
[Step 04 results](../../04_tsg_axis/README.md#formal-results) ·
[Step 09 selection](../../09_cpst_axis/README.md#formal-results)

Use `tenv_axes.open_tenv_axis_training("denoise-3k")` to open the materialized
release. The shared builder uses `--axis denoise-3k` for candidate publication
and replay. Its stable release path is
`cpst/03_data_train/03.5_tenv_denoise-3k/output`.
[Reader](../output/api/tenv_axes.py) ·
[Build and replay procedure](../README.md#tenv-training-axes)
