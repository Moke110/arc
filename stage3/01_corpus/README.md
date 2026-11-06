# 01_corpus — Heliste source-first cell viability

This step builds separate `output/baseline/` and `output/H-R/` corpora from
unchanged Heliste 2021 Dataset S1 and frozen PubChem chemical evidence in
[`pheno_data/Heliste-2021`](../../pheno_data/README.md), under
[ADR 0044](../../docs/adr/0044-stage3-source-first-corpora.md) and spec #52.

Use the `arc` environment (Python with `openpyxl`, `pyarrow`, RDKit; frozen RDKit
version is recorded in the identity evidence):

```bash
python stage3/01_corpus/execution/corpus.py build
python stage3/01_corpus/execution/corpus.py verify
python -m unittest discover -s stage3/01_corpus/execution/tests -v
```

A fresh checkout first acquires originals with
`python pheno_data/Heliste-2021/script/download.py`. The tracked local identity
snapshot is sufficient for offline building; only an intentional identity refresh
runs `python pheno_data/Heliste-2021/script/resolve_chemicals.py` with network
access. New reference evidence may change coverage and must be reviewed.
`--output /tmp/heliste-corpus` builds an alternate output. Existing output is
verified and never silently overwritten. Remove or choose a new directory after
an intentional input/software change.

Each condition contains equivalent `corpus.parquet` and `corpus.csv`,
`source_annotation.json`, and `chemical_resolution_report.json`. The common
`output/manifest.json` hashes outputs, frozen inputs and records coverage and
software versions. Verification checks hashes and reconstructs the entire release
from the archived workbook and independently checked reference evidence offline.
It rejects changed output even if its checksum is recomputed.

One record is one supplied numerical cell. `source_label` preserves the exact
original name; `canonical_smiles` is a deterministic isomeric RDKit representation
of the reference chemical form, preserving salts and stereochemistry without
inventing an active parent. `dose_um` is micromolar, `condition` is `baseline` or
`H-R`, and `phenotype` is `cell_viability`. `source_asset` and `source_cell` identify
the workbook cell. `rep` is null; `readout` is `{"source_reported": value}` (a typed
struct in Parquet, JSON in CSV). Distinct source labels are retained even when
they resolve to the same structure. Attribute combinations are not uniqueness keys.

Dataset S1 has 689 source-label rows, four doses (0.04, 0.2, 1 and 5 µM), and
5,512 numerical cells, 2,756 per condition. CTRL denotes chemically treated
normoxia and HYP denotes hypoxia followed by reoxygenation. Every cell is admitted
with a resolved structure or accounted for in the pending-resolution report.
Chemical labels and unique structures are counted separately. No statistical or
endpoint admission filter is applied.

The assay is normalized CTG ATP-dependent luminescence in HL-1 mouse atrial
cardiomyocytes: 24 h pretreatment, 48 h hypoxia and 24 h reoxygenation, with parallel
duration-matched normoxia, approximately 96 h elapsed and 0.01% DMSO controls.
Drug persistence and the exact normalization denominator remain unresolved.
Supplied values do not establish physical replicates or group-mean semantics;
vehicle measurements, SD, n, p-values, raw RLU and endpoint labels are absent.
These limitations do not exclude resolved chemical measurements. This corpus
makes no training or model-view readiness claim.

`execution/historical/` retains unchanged endpoint-gated release evidence. Its
path strings and scientific exclusions describe the superseded workflow; the old
implementation is recoverable from Git at commit `74f7e7a`.
