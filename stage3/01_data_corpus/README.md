# 01_data_corpus — experimental phenotype corpus

This step consumes frozen sources and declared curation evidence from
[`pheno_data`](../../pheno_data/README.md). Code and run records belong in
`execution/`; immutable releases belong in `output/`. The public offline
build/verification workflow is `execution/corpus.py` (tickets #31–#32).

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
python stage3/01_data_corpus/execution/release.py build --release-id v2
python stage3/01_data_corpus/execution/release.py verify --release-id v2
python -m unittest discover -s stage3/01_data_corpus/execution/tests -v
```

For a fresh checkout, acquire the five source items with each item's
`script/download.py` first; these are the only network operations. Materialize
the frozen local additions explicitly before the offline build:

```bash
python stage3/01_data_corpus/execution/sirenko_adapter.py --freeze-curation --output /tmp/sirenko-input.json
python stage3/01_data_corpus/execution/viability_adapters.py --prepare
```

The tracked Heliste local assessment is already present in a checkout.
Releases live under `output/corpus/<release-id>/`. `release.py bundle --output
/tmp/phenotype-input.json` exports the complete frozen input without publishing;
`corpus.py build --input ... --output ...` also supports an explicit release path.

Build and verify read frozen local assets only. Production sources must match the
active phenotype registry and include every original manifest asset; declared
local evidence is hashed separately. `--test-only` explicitly quarantines fixture
sources/releases and cannot admit them through the production workflow.

The tracer reads original F2/G2 values (26 and 27 peaks) and O2/P2 normalized
values without creating extra observations. Missing canonical identity and
matched vehicle statistics produce an explicit empty-after-QC partition and
comparison-level reasons. `sirenko_tracer.py` remains a minimal example;
`sirenko_adapter.py` extracts the full frozen 24 h workbook.

A release contains Parquet files plus input, schema, rule/software/code identity,
source registry, comparison/family, exclusion, coverage, manifest and checksum
sidecars. Rebuilding an identical payload is idempotent; different payloads need
a new release directory. Verification checks hashes, types, keys and evidence,
then reconstructs labels and logical outputs offline. Recorded software/code
versions are historical provenance; rule thresholds must match a supported rule
version. A numerical corpus remains separate from model-view or training status.

Statistics run at the original comparison before sample grouping. Numerical
source p-values are preferred; matching complete author BH results are reused.
Independent experimental units support two-sided Welch, paired designs support
two-sided paired t-tests, and known independent mean/sample-SD/n supports Welch
from summaries. Technical measurements are summarized by documented unit for the
test while corpus readouts remain unchanged. Unknown or degenerate inputs retain
an unavailable-test reason.

BH includes all eligible tests in the complete source × full condition × phenotype
× assay × original batch family before magnitude filtering. Exact declared family
membership prevents selective fragments from silently becoming complete families.
FC ≥ 1.2 or ≤ 0.8 requires adjusted p < 0.05 for ±1; strict interior FC gives 0
without inventing missing p; invalid FC or missing required p produces exclusion.
The input contract is documented in
[`execution/bundle-format.md`](execution/bundle-format.md).

The first-wave result and validation evidence are recorded in
[`execution/handoff.md`](execution/handoff.md) and the step-level
[`source_coverage.csv`](source_coverage.csv). Numerical zero means the operational
joint rule failed; missing p is preserved and does not become a claim of
statistical nonsignificance. Empty partitions retain their conditions and reasons.
