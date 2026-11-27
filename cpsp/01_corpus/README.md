# 01_corpus — Heliste source-first cell viability

This step builds separate `output/baseline/` and `output/H-R/` corpora from
unchanged Heliste 2021 Dataset S1 and frozen PubChem chemical evidence in
[`pheno_data/Heliste-2021`](../../pheno_data/README.md), under
[ADR 0044](../../docs/adr/0044-stage3-source-first-corpora.md) and specs #52, #54 (desalted-parent amendment), and #56 (magnitude direction labels).

Use the `arc` environment (Python with `openpyxl`, `pyarrow`, RDKit; frozen RDKit
version is recorded in the identity evidence):

```bash
python cpsp/01_corpus/execution/corpus.py build
python cpsp/01_corpus/execution/corpus.py verify
python -m unittest discover -s cpsp/01_corpus/execution/tests -v
```

A fresh checkout first acquires originals with
`python pheno_data/Heliste-2021/script/download.py`. The tracked local identity
snapshot is sufficient for offline building; only an intentional identity refresh
runs `python pheno_data/Heliste-2021/script/desalt_chemicals.py --acquire`
with network access. The immutable v1 reference snapshot is retained; v2 freezes
the desalted parents and additional reference evidence. New reference evidence may change coverage and must be reviewed.
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
of the desalted parent: RDKit selects the largest organic fragment, removes
counterions and solvates, and uncharges removable protonation. Stereochemistry,
isotopes and permanent charges remain. Active coordination-metal drugs retain
their metal and ligand components (even when the reference SMILES is disconnected),
removing only water and alkali counterions; their charges are preserved. Bacitracin
zinc instead uses an independently confirmed, stereochemically matching Bacitracin
parent. No tautomer or stereoisomer is selected.
Original salt/form structures remain in the frozen reference evidence. `dose_um` is micromolar, `condition` is `baseline` or
`H-R`, and `phenotype` is `cell_viability`. `source_asset` and `source_cell` identify
the workbook cell. `rep` is null; `readout` is `{"source_reported": value}` (a typed
struct in Parquet, JSON in CSV). `direction` is a derived integer: -1 below 80, 0 from 80 through 120 inclusive,
and +1 above 120. This is the user-selected ±20% magnitude rule, assuming each
arm's own 0.01% DMSO vehicle reference equals 100. The assumption is not an
explicitly confirmed source formula. Class 0 means within the magnitude band,
not statistical equivalence; these labels make no significance claim. The
versioned rule and normalization assumption are recorded in each source annotation.
Schema version 2 adds this column without changing the source readout.

Distinct source labels are retained even when
they resolve to the same structure. Attribute combinations are not uniqueness keys.

Dataset S1 has 689 source-label rows, four doses (0.04, 0.2, 1 and 5 µM), and
5,512 numerical cells, 2,756 per condition. CTRL denotes chemically treated
normoxia and HYP denotes hypoxia followed by reoxygenation. Every cell is admitted
with a resolved structure or accounted for in the pending-resolution report.
Chemical labels and unique structures are counted separately. Current coverage is
686 resolved labels and 663 unique desalted parents, with 2,744 admitted values
and 12 pending values per condition. Verteporfin and Formoterol hemifumarate have
unresolved stereochemical reference ambiguity; Sevelamer is a polymer without a
supported single finite structure. Desalting does not resolve these identities. No statistical or
endpoint admission filter is applied.

The assay is normalized CTG ATP-dependent luminescence in HL-1 mouse atrial
cardiomyocytes: 24 h pretreatment, 48 h hypoxia and 24 h reoxygenation, with parallel
duration-matched normoxia, approximately 96 h elapsed and 0.01% DMSO controls.
Drug persistence and the exact normalization denominator remain unresolved.
Supplied values do not establish physical replicates or group-mean semantics;
vehicle measurements, SD, n, p-values, raw RLU and source-provided endpoint labels
are absent. Derived magnitude directions are supplied under the assumption above.
These limitations do not exclude resolved chemical measurements. The direction column supplies a magnitude target; train/test splitting and model
views remain downstream work.

`execution/historical/` retains unchanged endpoint-gated release evidence. Its
path strings and scientific exclusions describe the superseded workflow; the old
implementation is recoverable from Git at commit `74f7e7a`.
