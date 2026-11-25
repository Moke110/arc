# Plate-matched training data

[03.1 tenv](03.1_tenv_hvg-deg-5k/README.md) publishes plate-pooled tenv on the
HVG-DEG 5k axis. [03.2 Tahoe](03.2_tahoe_hvg-deg-5k/README.md) exposes the
accepted Tahoe well-level data on the same gene axis. These source-specific
routes implement [Spec 5](https://github.com/Moke110/arc/issues/310).
The [acceptance receipt](log/312-source-specific-releases.json) records the
real tenv replay, Tahoe asset equality and sealed combined-view replay.

## Six tenv training axes

Each release follows the [03.1 contract](03.1_tenv_hvg-deg-5k/README.md) and
uses the same native tenv counts and [frozen folds](output/frozen-folds/). Each pools integer
counts on its selected axis before normalizing on that axis. The
[six-axis acceptance receipt](log/313-six-tenv-axes.json) records exact replay
and equality of the plate–parent–dose keys, raw membership and exclusions,
and all six fold assignments and strip reasons against HVG-DEG 5k.

| tenv axis key | Release |
| --- | --- |
| `hvg-deg-3k` | [03.3 HVG-DEG 3k](03.3_tenv_hvg-deg-3k/README.md) |
| `hvg-deg-4k` | [03.4 HVG-DEG 4k](03.4_tenv_hvg-deg-4k/README.md) |
| `hvg-deg-5k` | [03.1 HVG-DEG 5k](03.1_tenv_hvg-deg-5k/README.md) |
| `denoise-deg-3k` | [03.5 Denoise-DEG 3k](03.5_tenv_denoise-deg-3k/README.md) |
| `denoise-deg-4k` | [03.6 Denoise-DEG 4k](03.6_tenv_denoise-deg-4k/README.md) |
| `denoise-deg-5k` | [03.7 Denoise-DEG 5k](03.7_tenv_denoise-deg-5k/README.md) |

With `stage2/03_data_train/output/api` on `PYTHONPATH`, call
`open_tenv_axis_training("denoise-deg-3k")` in `training_data`. Pass a Step 03
root as the second argument to read releases materialized in another checkout.

Build and verify a candidate using the family and width documented for its
release. This example rebuilds Denoise-DEG 3k from the original folds and raw
corpus, verifies replay, then compares the accepted stable output against the
expected digest. Builders reject population drift before installing artifacts.

```bash
D=$PWD/stage2/03_data_train
CANDIDATE=$PWD/tmp/tenv-denoise-deg-3k/candidate
STAGING=$PWD/tmp/tenv-denoise-deg-3k/staging
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python "$D/execution/build_tenv_training.py" publish --family denoise-deg --size 3000 \
  --output "$CANDIDATE" --staging "$STAGING"
python "$D/execution/build_tenv_training.py" verify --family denoise-deg --size 3000 \
  --output "$CANDIDATE" --staging "$STAGING"
python "$D/execution/build_tenv_training.py" publish --family denoise-deg --size 3000 \
  --output "$D/03.5_tenv_denoise-deg-3k/output" --staging "$STAGING" \
  --expected-sha256 "$VERIFIED_SHA256"
```

Use the digest printed by verification for `VERIFIED_SHA256`. In a linked
worktree, `--raw` and `--reference` can point to the materialized native corpus
and HVG-DEG 5k output in another checkout.

## ADR 0085 paired release

[ADR 0085](../../docs/adr/0085-stage2-raw-and-training-data.md)
defines the earlier paired tenv and Tahoe pooled release. Its axis is the ordered
5,000 genes published by [02_gene_axis](../02_gene_axis/README.md). The builder
projects source-native counts to that axis, sums compatible replicates, then
applies CP10K and natural log1p to each pooled profile.

- `C` is the equal-weight mean of all normalized plate/cell-line controls.
- `M_C` contains one normalized control per plate and cell line. It pools DMSO
  wells for tenv and DMSO_TF wells for Tahoe bulk.
- `M_T` contains normalized treatment and Delta for each plate, cell line,
  canonical drug parent and measured dose. Delta uses that row's matching
  normalized plate/cell-line control.

The release publishes `tenv/` and `tahoe-bulk/` and binds both dataset
manifests to the same raw release and axis order. Each dataset directory holds
`training.h5`, the ordered genes, control and treatment row metadata, and
membership and exclusion ledgers. The treatment metadata carries six frozen source-specific
Step 07 fold labels. Cold-drug-dose rows in either excluded strip remain
available with a separate `held_parent` or `held_dose` reason.

Add `stage2/03_data_train/output/api` to the Python import path to read the
published dataset. Matrix reads accept at most 256 rows per call.

```python
from training_data import open_training

with open_training(source="tenv") as training:
    c = training.control()
    mc = training.read_controls(0, 8)
    treatment = training.read_treatments(0, 8)
    delta = training.read_treatments(0, 8, "delta")
    labels = training.treatments.iloc[:8]
```

Use `publish --output <scratch-directory>` to build a paired candidate from
the accepted raw, axis and frozen fold inputs. Use `--source tenv` or
`--source tahoe-bulk` for an isolated candidate. The default fold input is
the [active frozen fold release](output/frozen-folds/), which preserves the
original Step 07 member ledgers and acceptance receipt. Pass `--fold-root`
and `--fold-receipt` for an explicit replay input.

```bash
D=stage2/03_data_train
python "$D/execution/build_training.py" publish \
  --output "$PWD/tmp/training-release" --staging "$PWD/tmp/training-staging"
python "$D/execution/build_training.py" verify \
  --output "$PWD/tmp/training-release" --staging "$PWD/tmp/training-staging"
```

The [paired candidate acceptance](log/308-tahoe-paired-candidate.json) records
both release bindings, the Tahoe populations and six-fold checks, and exact
replay before the stable-path cutover.

Promote a verified pair to the stable output using the manifest digest printed
by the verification run. Promotion rechecks the candidate and rebuilds the
paired release from the pinned raw, axis and original Step 07 folds under
the output lock. It installs the verified rebuild when its manifest digest
matches the candidate's and preserves the public reader.

```bash
python "$D/execution/build_training.py" promote \
  --candidate "$PWD/tmp/training-release" --expected-sha256 "$VERIFIED_SHA256" \
  --staging "$PWD/tmp/training-staging"
```

The [archived well-level full view](03.1_hvg-deg-5k/README.md) supports existing
study inputs and fixed-prediction replay. The
[cutover receipt](log/309-cutover.json) records the path mapping and verified
training and replay digests.
