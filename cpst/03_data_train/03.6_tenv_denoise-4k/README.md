# 03.6 tenv training data on Denoise-DEG 4k

This release provides plate-matched tenv training profiles on the ordered
4,000-gene Denoise-DEG axis. The builder selects these genes from native
integer counts, pools DMSO wells by plate and treated wells by plate, canonical
drug parent and measured dose, then applies CP10K and natural `log1p` using each
pool's 4k count total. Treatment Δ subtracts the normalized DMSO profile from
the same plate. [Denoise-DEG axes](../../02_gene_axis/02.2_denoise-deg/README.md) · [tenv release contract](../03.1_tenv_hvg-5k/README.md)

## Accepted release

The accepted release has 16 plate controls and 1,383 treatment conditions,
with 112 exclusions among 4,510 raw observations. The six-axis acceptance
receipt records byte-exact replay and equality with the HVG-DEG 5k reference
for ordered condition keys, matched control references, membership and
exclusion ledgers, and all six frozen fold assignments and strip reasons.
The release manifest records the ordered axis, native counts, frozen fold
provenance and output file digests. [Release manifest](output/manifest.json) · [Six-axis acceptance receipt](../log/313-six-tenv-axes.json) · [Frozen tenv folds](../output/frozen-folds/tenv/manifest.json)

## Downstream use

Step 04 evaluated Denoise-DEG 4k in the six-axis TranSiGen comparison and
selected Denoise-DEG 5k. Step 09 evaluated it with CarPS-Transcriptome and
UniMol. The 4k arm came closest to 5k in full-axis Δ PCC, and Step 09 retained
5k, which led in both mean full-axis Δ PCC and Δ RMSE. [Step 04 selection](../../04_tsg_axis/README.md#gene-axis-selection) · [Step 09 formal results](../../09_cpst_axis/README.md#formal-results)

With `cpst/03_data_train/output/api` on `PYTHONPATH`, call
`tenv_axes.open_tenv_axis_training("denoise-4k")` to open the materialized
release. Follow the shared build and replay procedure with `--axis denoise-4k`
and use `cpst/03_data_train/03.6_tenv_denoise-4k/output` as the stable release
path. [Reader](../output/api/tenv_axes.py) · [Build and replay procedure](../README.md#tenv-training-axes)
