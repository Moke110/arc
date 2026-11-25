# tenv HVG-DEG 5k training

The source-specific tenv release pools the native integer counts on the
accepted HVG-DEG 5k axis. For each eligible plate, it sums DMSO wells.
For each plate, canonical drug parent and measured dose, it sums treated
wells. After pooling, each profile uses its own selected-axis total for
CP10K and natural `log1p`. Each perturbation's Delta is treatment expression
minus the normalized control of its plate. The overall tenv reference weights
all normalized plate controls equally.

| Output | Contents |
| --- | --- |
| `tenv-ref-vector.h5` | `C`, the overall control vector |
| `plate-ref-mtx.h5` | `M_C`, its plate keys and control row metadata |
| `pert-mtx.h5` | `M_T/treatment`, `M_T/delta`, row keys, `control_row` and six fold columns |
| `genes.csv`, `controls.csv`, `treatments.csv` | Ordered genes and row metadata |
| `members.csv`, `exclusions.csv` | The native observation identities contributing to each row or excluded with a reason |
| `manifest.json` | Input and payload digests, normalization and populations |

`cold_drug_0` through `cold_drug_2` and `cold_dose_drug_0` through
`cold_dose_drug_2` carry the three frozen Step 07 folds for each family.
The test partition in cold-dose-drug holds the intersection of the held
parents and dose. Train and valid exclude their union. A row belonging to
only one of the two held sets has an empty assignment and its reason in
the corresponding `_reason` column.

Add `stage2/03_data_train/output/api` to `PYTHONPATH` and call
`open_tenv_training()` from `training_data`. Pass an output directory to
read a separately materialized release. Read `control()`, `read_controls`
or `read_treatments` with at most 256 matrix rows per call. The `treatments`
table aligns to the rows returned by `read_treatments`.

```bash
ROOT=$PWD/stage2/03_data_train
BUILDER=$ROOT/execution/build_tenv_training.py
CANDIDATE=$PWD/tmp/tenv-hvg-5k/candidate
STAGING=$PWD/tmp/tenv-hvg-5k/staging
python "$BUILDER" publish --output "$CANDIDATE" --staging "$STAGING"
python "$BUILDER" verify --output "$CANDIDATE" --staging "$STAGING"
python "$BUILDER" publish --output "$ROOT/03.1_tenv_hvg-5k/output" \
  --staging "$STAGING" --expected-sha256 "$VERIFIED_SHA256"
```

Use the manifest digest printed by the verification command for
`VERIFIED_SHA256`. Verification independently rebuilds the release and
matches its manifest and all payload digests. Use `--raw` when the accepted
native corpus is in a different checkout.
