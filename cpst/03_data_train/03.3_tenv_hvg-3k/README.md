# 03.3 tenv training data on HVG-DEG 3k

This release provides plate-matched tenv training profiles for comparisons
on the ordered 3,000-gene HVG-DEG axis. It pools native integer counts on
those genes, applies CP10K and natural `log1p` using each profile's
3k library total, then subtracts the normalized same-plate DMSO control
from each normalized treatment to obtain Δ. The 3k library total can change
normalized values even for genes shared with the 5k release.
[HVG-DEG axes](../../02_gene_axis/02.1_hvg-deg/README.md) · [Tenv release contract](../03.1_tenv_hvg-5k/README.md)

## Accepted release

The release contains 16 plate controls and 1,383 plate–drug–dose treatment
conditions, with 112 of 4,510 source observations excluded. It passed
byte-exact replay. The control and treatment keys, matched-control assignments,
raw membership and exclusions, and all six frozen fold assignments and strip
reasons agreed with the HVG-DEG 5k reference.
[Release manifest](output/manifest.json) · [Six-axis acceptance receipt](../log/313-six-tenv-axes.json)

## Downstream use and access

The 3k release was evaluated in the six-axis TranSiGen comparison and later
in the eight-axis CarPS-Transcriptome comparison with UniMol. The latter
study selected Denoise-DEG 5k.
[Step 04 results](../../04_tsg_axis/README.md#formal-results) · [Step 09 selection](../../09_cpst_axis/README.md#formal-results)

Use `open_tenv_axis_training("hvg-3k")` from `tenv_axes` to open a materialized
release. For candidate publication and replay, use `--axis hvg-3k` with the
shared builder. The stable destination is
`cpst/03_data_train/03.3_tenv_hvg-3k/output`.
[Reader](../output/api/tenv_axes.py) · [Build and replay procedure](../README.md#tenv-training-axes)
