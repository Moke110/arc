# Plate-matched training data

[03.1 tenv](03.1_tenv_hvg-5k/README.md) publishes plate-pooled tenv on the
HVG-DEG 5k axis. [03.2 Tahoe](03.2_tahoe_hvg-5k/README.md) exposes the
accepted Tahoe well-level data on the same gene axis. These source-specific
routes implement [Spec 5](https://github.com/Moke110/carps/issues/310).
The [acceptance receipt](log/312-source-specific-releases.json) records the
real tenv replay, Tahoe asset equality and sealed combined-view replay.

## Six tenv training axes

Each release follows the [03.1 contract](03.1_tenv_hvg-5k/README.md) and
uses the same native tenv counts and [frozen folds](output/frozen-folds/). Each pools integer
counts on its selected axis before normalizing on that axis. The
[six-axis acceptance receipt](log/313-six-tenv-axes.json) records exact replay
and equality of the plate–parent–dose keys, raw membership and exclusions,
and all six fold assignments and strip reasons against HVG-DEG 5k.

| tenv axis key | Release |
| --- | --- |
| `hvg-3k` | [03.3 HVG-DEG 3k](03.3_tenv_hvg-3k/README.md) |
| `hvg-4k` | [03.4 HVG-DEG 4k](03.4_tenv_hvg-4k/README.md) |
| `hvg-5k` | [03.1 HVG-DEG 5k](03.1_tenv_hvg-5k/README.md) |
| `denoise-3k` | [03.5 Denoise-DEG 3k](03.5_tenv_denoise-3k/README.md) |
| `denoise-4k` | [03.6 Denoise-DEG 4k](03.6_tenv_denoise-4k/README.md) |
| `denoise-5k` | [03.7 Denoise-DEG 5k](03.7_tenv_denoise-5k/README.md) |

With `stage2_CarPS-VAE/03_data_train/output/api` on `PYTHONPATH`, call
`open_tenv_axis_training("denoise-3k")` in `training_data`. Pass a Step 03
root as the second argument to read releases materialized in another checkout.

Build and verify a candidate using the family and width documented for its
release. This example rebuilds Denoise-DEG 3k from the original folds and raw
corpus, verifies replay, then compares the accepted stable output against the
expected digest. Builders reject population drift before installing artifacts.

```bash
D=$PWD/stage2_CarPS-VAE/03_data_train
CANDIDATE=$PWD/tmp/tenv-denoise-3k/candidate
STAGING=$PWD/tmp/tenv-denoise-3k/staging
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python "$D/execution/build_tenv_training.py" publish --axis denoise-3k \
  --output "$CANDIDATE" --staging "$STAGING"
python "$D/execution/build_tenv_training.py" verify --axis denoise-3k \
  --output "$CANDIDATE" --staging "$STAGING"
python "$D/execution/build_tenv_training.py" publish --axis denoise-3k \
  --output "$D/03.5_tenv_denoise-3k/output" --staging "$STAGING" \
  --expected-sha256 "$VERIFIED_SHA256"
```

Use the digest printed by verification for `VERIFIED_SHA256`. In a linked
worktree, `--raw` and `--reference` can point to the materialized native corpus
and HVG-DEG 5k output in another checkout.

## Frozen inputs

The [active frozen fold files](output/frozen-folds/) retain the original tenv
member assignments and Step 07 receipt. Each release binds them to the native
corpus and selected axis. The `cold_dose_drug_*` columns leave held-parent and
held-dose strips unassigned and record their separate reasons.

The earlier ADR 0085 paired release and combined well-level view are in the
[Stage 2 local archive](../README.md#local-history).
