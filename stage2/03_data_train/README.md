# Plate-matched training data

[ADR 0085](../../docs/adr/0085-stage2-raw-and-training-data.md)
defines the pooled training release. The prediction axis remains the ordered
5,000 genes published by [02_gene_axis](../02_gene_axis/README.md). The builder
projects source-native counts to that axis, sums compatible replicates, then
applies CP10K and natural log1p to each pooled profile.

- `C` is the equal-weight mean of all normalized plate/cell-line controls.
- `M_C` contains one normalized DMSO control per tenv plate and cell line.
- `M_T` contains normalized treatment and Delta for each plate, cell line,
  canonical drug parent and measured dose. Delta uses that row's matching
  normalized plate/cell-line control.

The release manifest binds each dataset manifest to the same raw release and
axis order. Each dataset directory contains those matrices in `training.h5`,
the ordered genes, control and treatment row metadata, a membership ledger and
an exclusion ledger. The treatment metadata carries six frozen source-specific
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

Use the CLI with a separate destination to build and replay a candidate from
the accepted raw, axis and frozen fold inputs.

```bash
D=stage2/03_data_train
python "$D/execution/build_training.py" publish \
  --output "$PWD/tmp/training-release" --staging "$PWD/tmp/training-staging"
python "$D/execution/build_training.py" verify \
  --output "$PWD/tmp/training-release" --staging "$PWD/tmp/training-staging"
```
