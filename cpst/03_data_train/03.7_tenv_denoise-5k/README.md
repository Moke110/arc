# 03.7 tenv training data on Denoise-DEG 5k

This release supplies measured same-plate DMSO and treated expression on the
ordered 5,000-gene Denoise-DEG axis. Counts are projected onto the axis before
pooling DMSO wells by plate and treated wells by plate, canonical drug parent
and measured dose. Each profile uses its selected-axis count total for CP10K
followed by natural `log1p`. The released Δ subtracts the normalized matched
DMSO profile from the normalized treatment.
[Step 02.2 axis](../../02_gene_axis/02.2_denoise-deg/README.md) · [Tenv file contract](../03.1_tenv_hvg-5k/README.md)

## Accepted release

From 4,510 source observations, the release retains 16 pooled plate controls
and 1,383 pooled treatment conditions, with 112 source observations excluded.
The six-axis acceptance recorded byte-exact replay and equality with the
HVG-DEG 5k reference for the ordered condition keys, matched-control references,
raw membership, exclusions, and assignments and strip reasons for all three
cold-drug and three cold-dose-drug folds. The release manifest pins the axis,
raw counts, frozen folds and all eight output files.
[Release manifest](output/manifest.json) · [Acceptance receipt](../log/313-six-tenv-axes.json)

## Downstream use

Step 04 selected Denoise-DEG 5k from six axes for held-out Δ prediction with
TranSiGen. Steps 05 and 06 held the release fixed for the loss and architecture
studies. [Step 04 selection](../../04_tsg_axis/README.md#gene-axis-selection) ·
[Step 05 loss studies](../../05_tsg_loss/README.md) · [Step 06 architecture studies](../../06_tsg_shape/README.md)

Step 07 used the release for the CarPS-Transcriptome alpha, and Step 08 kept
it fixed for the chemical representation comparison. Step 09 selected
Denoise-DEG 5k again in the eight-axis comparison with UniMol. Step 10 used
this release's cold-dose-drug folds to compare dose inputs.
[Step 07 alpha](../../07_cpst_alpha/README.md) · [Step 08 representations](../../08_cpst_chem-emb/README.md#formal-results) ·
[Step 09 selection](../../09_cpst_axis/README.md#formal-results) · [Step 10 dose results](../../10_cpst_dose/README.md#formal-results)

## Read or reproduce

With `cpst/03_data_train/output/api` on `PYTHONPATH`, call
`open_tenv_axis_training("denoise-5k")` from `tenv_axes`. Pass the Step 03
root as the second argument to read a release materialized in another checkout.
For candidate publication and exact replay, use `--axis denoise-5k` and
`cpst/03_data_train/03.7_tenv_denoise-5k/output` as the stable release path.
[Named reader](../output/api/tenv_axes.py) · [Step 03 build and replay procedure](../README.md#tenv-training-axes)
