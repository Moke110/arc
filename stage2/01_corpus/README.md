# Source-native perturbation counts

The source-native corpus retains tenv wells, Tahoe plate/well/cell-line bulk
profiles, four additional cardiac perturbation deposits and their observed
native gene columns. Vehicles, technical blanks and unmatched observations
remain present. No axis projection, normalization or treated–control contrast
is applied to stored counts. Each source records measured common-axis coverage.
The builder rejects an incomplete common axis for tenv and Tahoe.

| Source | Native identity | Control metadata |
|---|---|---|
| tenv | GSE262419 plate well and resolved native HGNC gene | Same-plate vehicle. The paired-study DMSO group is kept separate from Media. |
| Tahoe bulk | Plate, well, cell line and deposited Ensembl/token gene | Same plate and cell line, vehicle DMSO_TF. |
| DToxS-iPSC | Separate conventional and DGE read counts by published sample and gene | Experiment design records with cell, plate, time and control state. |
| LINCS LDS-1593 | One DToxS-primary count matrix keyed by Plate, Well, Cell, Experiment, Dish, State and Time | Each original control context is retained with its publisher table. A match requires an identified compatible solvent. |
| sctox (GSE244740) | Published ScreenSeq well barcode and Ensembl row | Optimization and Screening plate identities kept distinct with separate DMSO and untreated groups. |
| GSE247541 | Published sample replicate and versioned Ensembl row | Published DMSO replicate group. |

XPert-L1000 has no verified native raw integer count publication in the
materialized source. Annotation is reference metadata with no observations.
The release manifest inventories seven dataset directories directly under
`output/`. Each directory holds a dataset manifest and a single H5AD with
native counts, canonical observations and native gene metadata. The tenv
directory also contains the probe map and source-build manifest. The
DToxS-primary directory contains `source_observations.csv`, which traces
842 merged observations to all 956 original rows across eight tables.
Counts for shared observations must agree after alignment by native gene ID.
The [original acceptance receipt](log/acceptance.json) records the v1
populations, coverage, source replay and payload digests.

```python
from native_corpus import open_corpus

with open_corpus("stage2/01_corpus/output") as corpus:
    source = corpus.source("tenv")
    native_columns = source.native_columns([9, 20, 10])
    counts = source.counts(0, 64)[:, native_columns]
    control_rows = source.matched_control_rows(0)
```

Add `stage2/01_corpus/output/api` to the Python import path. Source keys, native
observation metadata, output digests and ordered common-axis mappings are
published in `manifest.json` and the bounded reader. `native_columns` fails
on unmeasured genes. Windows have a 256-row maximum. Read subsequent windows
through `batches()`. The reader checks a source's payload before allowing reads.

```bash
python stage2/01_corpus/execution/build_corpus.py publish
python stage2/01_corpus/execution/build_corpus.py verify
```

The publisher replays tenv and native Tahoe counts against pinned external inputs.
Tahoe aggregation checkpoints can resume after source pins are checked again.
Verification rebuilds a separate candidate, compares every sealed digest, and
retains the accepted release. A release installed with `migrate` replays its
layout from the sealed v1 generation in `historical/native-v1`. Migration
compares native matrices, observations and gene mappings for the six unchanged
datasets before installation. Existing full 5k views keep the manifest digest
of their original derivation. New view builds pin the current corpus digest. A separately staged native Tahoe
build can be read with `--tahoe-built` and goes through the same independent
source replay before it can supply the release. For linked Git worktrees,
`--materialized-root` can select the checkout holding the immutable inputs.
