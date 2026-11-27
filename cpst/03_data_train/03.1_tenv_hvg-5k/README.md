# 03.1 tenv training data on HVG-DEG 5k

This release gives each tenv perturbation condition a measured DMSO control
from the same plate. It uses [source-native integer counts](../../01_data_raw/README.md)
and the [ordered HVG-DEG 5k axis](../../02_gene_axis/02.1_hvg-deg/README.md).
The builder selects the 5,000 genes before summing DMSO wells by plate and
treated wells by plate, canonical drug parent and measured dose. It applies
CP10K and natural `log1p` to each pooled profile using that profile's
selected-axis count total. The released Δ subtracts the normalized matched
control from the normalized treatment. The dataset-wide reference is the
equal-weight mean of the normalized plate controls.

## Accepted release

Of 4,510 source observations, 112 were excluded, including 80 incompatible
vehicle observations and 32 technical blanks. The release contains 16 pooled
plate controls and 1,383 pooled treatment conditions. Candidate publication,
independent byte-exact replay, partition checks and promotion all passed.
The candidate, replay and promoted manifest digests agree.
[Release manifest](output/manifest.json) · [Acceptance receipt](../log/312-source-specific-releases.json)

The treatment rows carry three cold-drug folds and three cold-dose-drug
folds from the [frozen tenv inputs](../output/frozen-folds/tenv/manifest.json).
Cold-dose-drug assigns their held-parent and held-dose intersection to test.
Train and validation use rows outside their union. Rows matching only one
held criterion have an empty assignment and a `held_parent` or `held_dose` value in the
corresponding `_reason` column.

| File in the materialized `output/` | Contents |
| --- | --- |
| `tenv-ref-vector.h5` | Dataset-wide control vector `C` |
| `plate-ref-mtx.h5` | Plate control matrix `M_C` and its row metadata |
| `pert-mtx.h5` | `M_T/treatment` and `M_T/delta`, condition keys, matched `control_row` and fold metadata |
| `genes.csv`, `controls.csv`, `treatments.csv` | Ordered genes and matrix row metadata |
| `members.csv`, `exclusions.csv` | Source observation membership and exclusion reasons |

## Downstream use and verification

Step 03 uses this release to check that new tenv axes preserve condition
keys, raw membership, exclusions and fold assignments. HVG-DEG 5k was
included in the Step 04 TranSiGen and Step 09 CarPS-Transcriptome axis
comparisons. Step 09 selected Denoise-DEG 5k for CarPS-Transcriptome.
[Step 03 training axes](../README.md#tenv-training-axes) · [Step 04 comparison](../../04_tsg_axis/README.md) · [Step 09 selection](../../09_cpst_axis/README.md#formal-results)

Use `open_tenv_training()` from `training_data` with
`cpst/03_data_train/output/api` on `PYTHONPATH`. Pass a materialized release
directory to override the default. `control()` returns `C`. Read up to 256
consecutive rows per call with `read_controls(start, stop)` or
`read_treatments(start, stop, layer="delta")`. The treatment reader defaults
to normalized treatment expression, and both layers align to the
`treatments` table. [Reader implementation](../output/api/training_data.py)

Follow the Step 03 candidate publication and replay procedure with
`--axis hvg-5k` and use `cpst/03_data_train/03.1_tenv_hvg-5k/output` for
the stable release. The `verify` command independently rebuilds the release
and compares the manifest and every payload digest.
[Build and replay procedure](../README.md#tenv-training-axes) · [Builder](../execution/build_tenv_training.py)
