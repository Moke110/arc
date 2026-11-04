# 01_data_corpus — experimental phenotype corpus

This step consumes frozen sources and declared curation evidence from
[`pheno_data`](../../pheno_data/README.md). Code and run records belong in
`execution/`; immutable releases belong in `output/`. The public offline
build/verification workflow is implemented under tickets #31–#32.

Each complete condition/phenotype partition keeps source, canonical drug, molar
dose and numeric endpoint in its sample identity. Readouts remain numerical
lists with aligned measurement evidence. MTT and CCK8 share `cell_viability`
while retaining assay methods. Scientifically unevaluable comparisons remain in
exclusion evidence and may yield explicit empty partitions.

The authoritative design is
[ADR 0041](../../docs/adr/0041-formal-corpus-statistics-export-and-coverage.md)
and the [formal corpus design](../../docs/research/stage3-formal-corpus-design.md).
Small files remain visible; the later model-view 500 target and training are
separate assessments.

Run with the `arc` environment (`pyarrow`, `scipy`, `xlrd`):

```bash
python stage3/01_data_corpus/execution/sirenko_tracer.py --output /tmp/sirenko-tracer.json
python stage3/01_data_corpus/execution/corpus.py build --input /tmp/sirenko-tracer.json --output /tmp/sirenko-tracer-release
python stage3/01_data_corpus/execution/corpus.py verify --release /tmp/sirenko-tracer-release
python -m unittest discover -s stage3/01_data_corpus/execution/tests -p test_workflow.py -v
```

Build and verify read frozen local assets only. Production sources must match the
active phenotype registry and include every original manifest asset; declared
local evidence is hashed separately. `--test-only` explicitly quarantines fixture
sources/releases and cannot admit them through the production workflow.

The tracer reads original F2/G2 values (26 and 27 peaks) and O2/P2 normalized
values without creating extra observations. Missing canonical identity and
matched vehicle statistics produce an explicit empty-after-QC partition and
comparison-level reasons. It demonstrates the public workflow; the full source
adapter is delivered by #33.

A release contains Parquet files plus input, schema, rule/software/code identity,
source registry, comparison/family, exclusion, coverage, manifest and checksum
sidecars. Rebuilding an identical payload is idempotent; different payloads need
a new release directory. Verification checks hashes, types, keys and evidence,
then reconstructs labels and logical outputs offline. Recorded software/code
versions are historical provenance; rule thresholds must match a supported rule
version. A numerical corpus remains separate from model-view or training status.
