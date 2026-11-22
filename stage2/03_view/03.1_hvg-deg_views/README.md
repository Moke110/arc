# Raw-first HVG-DEG expression views

The versioned release uses the ordered HVG-DEG 5,000-gene prediction axis.
Each treatment profile and its control pool are selected on the native corpus.
The builder sums selected raw integer counts independently for treatment and
the matched vehicle union, counting each matched vehicle once. It then applies
CP10K and natural log1p using each side's own selected-axis total. Delta is
the treatment expression minus the control expression.

| Family | Treatment and control profile |
| --- | --- |
| `full` | Each matched, eligible treated observation and each eligible individual vehicle in filtered native source order. Every row has its matched aggregate control. This is the default release. |
| `merge` | Treatment counts summed per cell line, parent and dose. The control includes the union of the vehicles matched to those treatment members. |
| `solo` | One treatment member per condition, selected with the original frozen source-aware draw. Includes that member's matched control. |
| `solo-deepest` | tenv only. The condition member with the highest selected-axis library total, tied by profile ID, and its matched control. |
| `full-pooled-label` | Full profile expression, control and row order with each treated Delta drawn from its merge condition. Individual vehicle Delta remains unchanged. |

The H5AD files for the full family contain float32 X, control and delta, with
`obs.role`, stable source profile IDs and ordered `var.gene_index`. The
pooled-label payload stores only Delta and resolves X and control from the
same release's full payload. Per-arm membership ledgers record treatment
and matched vehicles by profile identity. Per-source exclusion ledgers record
unmatched observations and source eligibility exclusions. The manifest seals
every file and binds the corpus and selected-axis manifest provenance.

Add `stage2/03_view/03.1_hvg-deg_views/output/api` to the Python import path
and use the bounded public reader:

```python
from expression_view import open_view

with open_view(source="tenv") as view:
    first = view.read(0, min(64, view.shape[0]), layer="control")
```

A separate verified release path can be supplied as the first argument.
`corpus_sha256` and `axis_sha256` assert input generation identity.
`view.manifest_sha256` seals the entire release, and
`view.manifest["outputs"]` records each arm's payload digest.
Batches support `X`, `control` or `delta` up to 256 rows.

```bash
STEP=stage2/03_view/03.1_hvg-deg_views
python "$STEP/execution/build_views.py" publish --corpus /path/to/native-corpus
python "$STEP/execution/build_views.py" verify --corpus /path/to/native-corpus
python -m unittest discover -s "$STEP/execution/tests" -v
```

Publish and verify both use the same pinned corpus and axis, and verification
rebuilds separately and compares bytes. An isolated tenv rehearsal may use
`--source tenv --destination /path/to/rehearsal`.

The prior recipes, payloads, 45-fit evidence and sealed replay remain at
[the historical view-study boundary](../../tsg_view/README.md) during migration.
