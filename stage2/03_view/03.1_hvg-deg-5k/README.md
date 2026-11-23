# Raw-first HVG-DEG expression views

The versioned release uses the ordered HVG-DEG 5,000-gene prediction axis.
The full release retains eligible matched treatment profiles and individual
vehicles in native source order. For each profile, the builder sums the
matched vehicles' selected-axis raw counts, counting each vehicle once.
It applies CP10K and natural log1p separately to treatment and control using
their own selected-axis totals. Delta is treatment expression minus control.

The flat `output/release` directory contains one H5AD, membership ledger and
exclusion ledger for each of tenv and Tahoe bulk. The H5ADs contain float32
X, control and delta, with stable profile IDs, `obs.role`, and ordered
`var.gene_index`. The versioned manifest seals all six assets, records the
original H5AD schema, and binds the corpus and prediction axis.

Add `stage2/03_view/03.1_hvg-deg-5k/output/api` to the Python import path
and use the bounded public reader:

```python
from expression_view import open_view

with open_view(source="tenv") as view:
    first = view.read(0, min(64, view.shape[0]), layer="control")
```

A separate verified release path can be supplied as the first argument.
`corpus_sha256` and `axis_sha256` assert input generation identity.
`view.manifest_sha256` seals the entire release, and
`view.manifest["outputs"]` records each source's payload digests.
Batches support `X`, `control` or `delta` up to 256 rows.

```bash
STEP=stage2/03_view/03.1_hvg-deg-5k
python "$STEP/execution/build_views.py" publish --corpus /path/to/native-corpus
python "$STEP/execution/build_views.py" verify --corpus /path/to/native-corpus
python -m unittest discover -s "$STEP/execution/tests" -v
```

Publish and verify both use the same pinned corpus and axis, and verification
rebuilds separately and compares bytes. An isolated tenv rehearsal may use
`--source tenv --destination /path/to/rehearsal`.

The prior 45-fit study and its acceptance records are documented in
[the historical view-study boundary](historical/README.md).

The full-only migration carries the six accepted assets forward byte-for-byte.
The [migration receipt](log/full-only-migration.json) maps their paths and
digests. The [preceding generation](historical/raw-first-v1/manifest.json)
retains the previous release inventory.

An existing full view can be replayed from the equivalent reorganized corpus
while retaining the digest of its original derivation. For a portable checkout
where the sealed v1 archive lives elsewhere, pass
`--historical-corpus <v1-corpus-directory>` to `build_views.py verify`.
