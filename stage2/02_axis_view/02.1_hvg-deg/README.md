# HVG-DEG prediction axes

The 2,000-, 3,000-, 4,000- and 5,000-gene axes retain the exact ADR 0073
ordering and symbols. Each entry is an original index on the 18,164-gene
common axis. The release contains the descriptive tenv-HVG ranking, tier-first
DEG1000 and full-union disease fusion tables, and cross-disease candidate,
composition and coverage evidence.

[The release manifest](output/release/manifest.json) binds those files to the
common gene identity and original selection method. Publication reuses the
existing method and independently compares its selected genes and symbols with
the accepted CSVs. Original Bases files and presentation artifacts remain
available for historical replay.

Add `output/api` to the Python path to use `prediction_axis.open_axes`. A
prediction axis projects a source-native corpus with complete common coverage.

```python
from native_corpus import open_corpus
from prediction_axis import open_axes

axis = open_axes().axis(5000)
with open_corpus(corpus_path) as corpus:
    projection = axis.project(corpus, "tenv")
    counts = projection.counts(0, 64)
```

`counts` preserves integer values and prediction-axis column order. Projection
also supports bounded `batches` and requires the corpus to stay open.

Publish from a checkout containing the original materialized HVG evidence.
Publication checks the canonical native corpus by default. Pass `--corpus`
when that release is materialized elsewhere.
Verification uses the bundled ranking and fusion evidence. The full verify
command also checks the native corpus and complete coverage for tenv and
Tahoe. Omit `--corpus` for an axis-only replay. Other claimed complete sources
are checked when present.

```bash
STEP=stage2/02_axis_view/02.1_hvg-deg
python "$STEP/execution/build_axes.py" publish --source-root /path/to/materialized/arc
python "$STEP/execution/build_axes.py" verify --corpus /path/to/native-corpus
python -m unittest discover -s "$STEP/execution/tests"
```
