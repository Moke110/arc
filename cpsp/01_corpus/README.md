# 01_corpus — Heliste cell viability

CarPS-Phenotype uses Heliste 2021 Dataset S1 as the source for two separate
`cell_viability` corpora. `baseline` comes from chemically treated normoxia
(`CTRL`), and hypoxia–reoxygenation (`H-R`) comes from `HYP`. The CellTiter-Glo
assay measures ATP-dependent luminescence in HL-1 mouse atrial cardiomyocytes.
The H-R schedule comprises 24 hours of pretreatment, 48 hours of hypoxia and
24 hours of reoxygenation. Parallel normoxia is duration-matched at about
96 hours. See the [source-first decision](../../docs/adr/0044-stage3-source-first-corpora.md)
and [phenotype source archive](../../pheno_data/README.md).

## Build the two arm corpora

Dataset S1 contains 689 source labels measured at 0.04, 0.2, 1 and 5 µM in
each arm. Frozen PubChem evidence supports 686 labels mapped to 663 unique
desalted parent structures. The release accounts for all 5,512 numerical
cells through the arm exports and their pending-identity reports.

| Arm | Source cells | Admitted | Pending chemical identity |
| --- | ---: | ---: | ---: |
| Baseline | 2,756 | 2,744 | 12 |
| H-R | 2,756 | 2,744 | 12 |

Each admitted row retains the original source label, dose, workbook cell and
reported value in `readout.source_reported`. Aliases remain separate rows even
when they share a desalted parent. Replicate index `rep` is null because the
workbook leaves replicate identity and individual-versus-group-mean semantics
unresolved. Admission depends on supported chemical identity. The pending
labels are Verteporfin (Visudyne) and Formoterol hemifumarate, whose evidence
leaves stereochemistry unresolved, and Sevelamer HCl, a polymer without a
supported single finite structure. Related source files and both frozen
chemical snapshots are recorded in the [Heliste provenance log](../../pheno_data/Heliste-2021/author_log.json).

## Derive magnitude direction

Schema version 2 records `direction` as −1 below 80, 0 from 80 through 120
inclusive and +1 above 120. This ±20% rule assumes each arm's own 0.01% DMSO
vehicle reference equals 100. The paper confirms DMSO normalization and leaves
the exact oxygen-specific denominator unresolved. Class 0 indicates a value
within the magnitude band. The labels describe magnitude only, and the source
annotations record the assumption and unresolved drug persistence. See the
[baseline annotation](output/baseline/source_annotation.json) and
[H-R annotation](output/H-R/source_annotation.json).

## Verify the release and use it in the Morgan baseline

Each arm provides `corpus.parquet`, `corpus.csv`, `source_annotation.json` and
`chemical_resolution_report.json`. The schema version, coverage, dependency
versions and code, input and output hashes are pinned in the
[corpus manifest](output/manifest.json). The eight arm files and source-input
hashes match the pre-migration release. CPSP acceptance verified the release
against an independent offline rebuild, and the manifest transition is pinned in
the [Morgan release record](../02_base_morgan/execution/release_record.json).
The separate [endpoint-gated release record](execution/historical/release_record.json)
describes the superseded corpus.

Use the [step instructions](AGENTS.md) to run the offline [corpus CLI](execution/corpus.py).
The [Morgan baseline](../02_base_morgan/README.md) verifies this release before
preparing both arms for cold-drug, cold-dose and warm-start evaluation. Its
direction target comes from the corpus rule above.
