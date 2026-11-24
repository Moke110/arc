# Plate-matched training data

[ADR 0085](../../docs/adr/0085-stage2-raw-and-training-data.md)
defines the pooled training release. The prediction axis remains the ordered
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
the [Step 07 frozen fold archive](../07_tsg_dose/README.md#frozen-folds-for-stage-2-training).
To use the original materialized Step 07 checkpoint, pass
`--fold-root stage2/07_tsg_dose/execution/checkpoints/round1-axis-local-inputs`.

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
