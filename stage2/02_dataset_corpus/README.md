# 02_dataset_corpus — Stage 2 dataset corpus

**Status:** the five corpora (four CM + `tahoe.h5ad`), the seven crosswalk tables
and the coverage README are built (tickets #14–#16), and the `--verify` gate suite
**G1–G6 PASS** — including G4, which reproduces 5/5 frozen v1 view content hashes
(ticket #17; `execution/build_record/verify.json`). The CM-schema view is owned by
`05_prnet_pretrain/` (ticket #18), and the docs/ADR sweep is ticket #20. The
governing decision is [ADR 0033](../../docs/adr/0033-stage2-dataset-corpus.md); the
artifact record (coverage, provenance, renames, schemas, gates) is
[`output/corpus/v2/README.md`](output/corpus/v2/README.md).

## Purpose

One step owns every Stage 2 source dataset as an information-maximal h5ad
corpus: the substrate between a raw download and a model's training input. Model
adaptation (pharmacology tiers, gene-space projection, author shims) belongs to
the consumer that builds a view, not here (ADR 0033 §5). `stage2/04_tahoe_data/`
was folded into this step and removed; `05`/`06` keep their numbers, so 04 is a
gap in the sequence, not a renumber.

Dataset ids are forward-only renames (ADR 0033 §9): `sctox` (was
`screen_ctox`), `tenv` (was `temp_env`), and `dtoxs_primary`/`dtoxs_ipsc` (the
merged `dtoxs_cm` splits back into its two arms). Frozen records keep their old
strings; ticket #16's README carries the rename table.

## Layout

```
02_dataset_corpus/
├── AGENTS.md  README.md
├── execution/                        # code, intermediates, build records
│   ├── build_corpus.py               # single entry point (tickets #15/#16)
│   ├── cmdata/                       # CM library
│   │   ├── paths.py sources.py genes.py drugs.py annotation.py
│   │   ├── corpus.py qa.py primary_doses.py text.py
│   │   ├── corpus_h5ad.py            # native-axis corpus builder/writer
│   │   └── tables.py                 # seven tables + coverage README (ticket #16)
│   ├── tahoe/build_pseudobulk.py     # Tahoe corpus builder (moved from 04)
│   ├── _parts/                       # 14 per-plate CSR parts + _state (local)
│   ├── build_record/                 # tracked QA/manifests (#14/#15/#16)
│   └── tests/                        # test_cmdata.py, test_corpus_contract.py,
│                                     # test_tables_contract.py
└── output/                           # artifacts a later step reads
    └── corpus/v2/
        ├── dtoxs_primary.h5ad  dtoxs_ipsc.h5ad  sctox.h5ad  tenv.h5ad
        ├── tahoe.h5ad                # 67,018 × 62,710 CSR int32 (25.8 GB, local)
        ├── README.md                 # coverage + provenance + gates (#16)
        └── tables/                   # seven crosswalk tables (#16)
```

## Corpus contract (ADR 0033 §2–3, §10)

| corpus | rows (pert / control / blank) | chemicals | native gene axis | container |
|---|---|---|---|---|
| `dtoxs_primary` | 779 (677 / 102 / 0) | 47 | 26,328 | dense `int32` X + dense `float32` layer |
| `dtoxs_ipsc` | 312 (296 / 16 / 0) | 37 | 26,328 | same |
| `sctox` | 3,072 (2,520 / 552 / 0) | 42 | 31,880 | same |
| `tenv` | 4,510 (4,206 / 272 / 32) | 460 | 19,613 | same |
| `tahoe` | 67,018 (65,618 / 1,400 / 0) | 379 | 62,710 | CSR `int32` + CSR `float32` layer |

Shared semantics: `X` = raw source counts **summed per resolved HGNC id**
(non-negative `int32`, dense for the four CM corpora); one
`layers["log1p_cp10k"] = log1p(1e4·X/row_sum)` `float32`; `var` = that source's
native resolved axis; `uns` empty; `obsm/varm/obsp/varp` empty; no
`pert_idx`/`cell_idx`/`dose_idx`; no Δ; no pharmacology tier filter; one row per
raw well (Tahoe: one pseudobulk profile per plate × well × cell line). The
18,164-gene four-source intersection is **not** a corpus universe any more: it
is the CM projection target owned by `tables/genes.csv` (ticket #16).

### Native gene axis and the 18,164 projection

The axis is the unique resolved HGNC ids of the source, in first-appearance
order over that source's raw gene ids, each represented by its canonical HGNC
symbol (`HgncIndex.ref(hgnc_id).symbol`). Duplicate raw ids that resolve to one
HGNC id are summed. `var` = `symbol` index plus `hgnc_id`, `entrez_id`,
`ensembl_id`, `match_type` (resolution type of the first raw id that introduced
the gene), so a consumer can rebuild the 18,164 projection with an identity gene
map: join `hgnc_id` to `tables/genes.csv` and sum. Rebuilding the five frozen v1
views this way reproduces their `content_sha256` 5/5 (the formal G4 gate in
`--verify`; `execution/build_record/verify.json`) — the corpus stores exactly the
integers the v1 writer produced before it projected.

### `obs` (four CM corpora, 27 columns, index = `well`)

`well, dataset_id, compound_name, drug_key, inchikey, smiles, dose_um,
dose_source, time_h, is_control, vehicle, plate, control_group, cell_type,
donor, replicate, atc, moa, mesh_classes, pharm_atc_moa, pharm_mesh,
annotation_source, lib_size, n_genes_detected, cell_line, n_raw_rows, is_blank`
— ADR 0026's 24 columns in their v1 order plus three additions.

- `dataset_id` uses the v2 short ids (`sctox`, `tenv`); `control_group =
  <dataset_id>::<plate>` uses the same new id. Ticket #17's reconstruction
  remaps `sctox`→`screen_ctox` and `tenv`→`temp_env`.
- `well` is the loader's string **verbatim**, so it keeps the historical id
  prefix (`screen_ctox::…`, `temp_env::…`). This is deliberate provenance: the
  G4 gate remaps `dataset_id`/`control_group`, not the well identity.
- `donor` is unchanged (`NA` for the dtoxs-primary and the GEO arms). The
  recovered raw cell-line label is backfilled into `cell_line`: `A/B/D/E` and
  `PC-A…PC-E` for `dtoxs_primary` (the only place the label exists),
  `MSN02-04R-CM`/`MSN09-04R-CM` for `dtoxs_ipsc`, `NA` for `sctox`/`tenv` (the
  sources name no per-well line).
- `n_raw_rows` records the raw design rows collapsed into that well before the
  loader deduplicates: `dtoxs_primary` 956 → 842 wells (36 keys reused in the
  experiment-5 sub-batches; 26 keys × 5 rows, 10 × 2 rows), 1 for every well of
  every other source.
- `lib_size`/`n_genes_detected` are native-axis values. A view consumer that
  projects onto 18,164 must recompute both from the projected matrix; the v1
  view hash depends on that, not on the corpus values.

Row set: vehicle controls + technical blanks + every SMILES-resolvable
small-molecule compound; no pharmacology tier is applied. The 8
non-small-molecule perturbagen keys (143 wells: four DToxS biologics, four iPSC
antibodies, the `ET-1` peptide) stay excluded with reasons, and the 32 `tenv`
technical blanks are retained rows (`is_blank=True`, never counted as
chemicals) — both recorded in
`execution/build_record/excluded_compounds.csv` (promoted into
`output/corpus/v2/tables/` by ticket #16).

## Building

From the repository root (the repo root is found by walking up to `.git`, never
by counting `../..` — ADR 0030 §4):

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
$PY stage2/02_dataset_corpus/execution/build_corpus.py --check          # raw source counts
$PY stage2/02_dataset_corpus/execution/build_corpus.py --corpus all     # four corpora, one at a time
$PY stage2/02_dataset_corpus/execution/build_corpus.py --corpus tenv    # one corpus
$PY stage2/02_dataset_corpus/execution/build_corpus.py --corpus dtoxs_ipsc --repeat-check
$PY stage2/02_dataset_corpus/execution/build_corpus.py --tables         # seven tables + coverage README
$PY stage2/02_dataset_corpus/execution/build_corpus.py --all            # four corpora, then the tables
```

Each build loads the four raw sources (the compound table is global, so
identity cannot depend on which corpus is being written), checks
`MemAvailable >= 12 GiB`, builds and validates one corpus at a time (G1
structure + row/axis/chemical counts + a raw-source spot check; any failure
aborts), then writes it with `track_times=False`, gzip and an atomic
`.tmp`→rename. Equivalence is content-level (ADR 0033 §7): `content_sha256` is
the v1 recipe minus `uns` (`X` bytes | layer bytes | `obs.to_csv(index=False)` |
`var.to_csv()`); file sha256/bytes are recorded separately.

Records in `execution/build_record/`:

| file | contents |
|---|---|
| `corpus_build.json` | per corpus: rows / perturbed / control / blank / chemicals / axis, source and `n_raw_rows` provenance, dtypes, obs/var columns, `content_sha256`, file sha256, bytes, wall time, peak RSS, gate results; plus the excluded-compound accounting |
| `excluded_compounds.csv` | 8 non-small-molecule keys (143 wells) + the 32 blank wells, with reasons |
| `corpus_reproducibility.json` | gate G5: two builds of one corpus, identical file sha256 |
| `corpus_coverage.json` | ticket #16: the five coverage counts per corpus + the 24 h sub-row, the seven table row counts/hashes, the README hash, the G2/G3 results and the frozen v1 identity-equality checks |

Measured build (2026-11-04, `/usr/bin/time -v`, one corpus per process):

| corpus | rows × axis | chemicals | bytes | wall | peak RSS | file sha256 (first 16) |
|---|---|---|---|---|---:|---:|---|
| `dtoxs_primary` | 779 × 26,328 | 47 | 47,402,239 | 2.4 s | 7.66 GiB | `a29540f403cb407c…` |
| `dtoxs_ipsc` | 312 × 26,328 | 37 | 18,902,113 | 1.6 s | 7.66 GiB | `16e846212c9fc20b…` |
| `sctox` | 3,072 × 31,880 | 42 | 164,875,596 | 11.2 s | 7.66 GiB | `c254d9096c59b7a2…` |
| `tenv` | 4,510 × 19,613 | 460 | 255,148,955 | 11.5 s | 7.66 GiB | `a3d521f0ca0fbc73…` |

The peak RSS is dominated by loading the four raw sources (~56 s per process),
not by the matrix build. G5 PASS: `dtoxs_ipsc` written twice → identical file
sha256 (`16e84621…`). `--verify` runs G1–G6, including the v1 reconstruction gate,
in 76.8 s at 7.66 GiB peak RSS; **all six gates PASS** — G4 reproduces the 5/5
frozen v1 `content_sha256` values (full record:
`execution/build_record/verify.json`).

## Tables and coverage README (ticket #16)

`output/corpus/v2/tables/` holds the seven crosswalk tables and
[`output/corpus/v2/README.md`](output/corpus/v2/README.md) is the artifact
record (corpus coverage five-counts + 24 h sub-row, provenance/hashes, the
rename table, per-corpus `obs`/`var` and the container-asymmetry note, gate
results). Both are generated by `--tables`, which reads the four CM corpora /
`tahoe.h5ad` and re-resolves the raw source gene ids (the corpora keep only
resolved genes), so it pays the same global source load as a corpus build.

| table | rows | content |
|---|---:|---|
| `genes.csv` | 18,164 | the frozen CM projection axis; byte-identical to `_frozen/cm-data-v1/corpus/v1/genes.csv` |
| `gene_crosswalk.csv` | 269,067 | every source native gene id → HGNC symbol: 143,647 four-source rows (identical to the frozen v1 table) + 125,420 Tahoe rows (`token_id` + `ensembl_id`, one row per namespace per gene) |
| `compounds.csv` | 532 | 525 `drug_key` + identity + `atc/moa/mesh_classes` + the three tier booleans + `datasets`/`n_wells`; byte-identical to the frozen v1 table |
| `excluded_compounds.csv` | 9 | 8 non-SMILES keys (143 wells) + the blank row (32 wells), with reasons |
| `source_gene_sets.csv` | 104,149 | per-source resolved gene sets; byte-identical to the frozen v1 table |
| `metadata_cm.csv` | 8,673 | the four CM corpora's `obs` flattened, with `dataset_id` + `legacy_dataset_id` + source file/accession; join-ready on `well` |
| `metadata_tahoe.csv` | 67,018 | the Tahoe `obs` flattened, `dataset_id=tahoe`, `is_blank=False`; join-ready on `profile_id` |

G2 (README coverage counts == `corpus_coverage.json`) and G3 (crosswalk covers
every `var` id of all five corpora; the 18,164 symbols are a subset of each CM
`var`; `compounds.csv` covers every perturbed `drug_key`) are asserted before
the tables are written and again by
`execution/tests/test_tables_contract.py`. Measured `--tables`: 23.7 s table
build after the ~56 s source load, `/usr/bin/time -v` peak RSS 7.78 GiB.

Tests (contract tests skip when the local-only h5ads are absent):

```bash
$PY -m unittest discover -s stage2/02_dataset_corpus/execution/tests -v
```

## Moved in ticket #14

- `02_cm_data/execution/` → `execution/` (code; the five v1 views moved to
  `stage2/_frozen/cm-data-v1/` in ticket #13).
- `04_tahoe_data/execution/build_pseudobulk.py` → `execution/tahoe/`;
  `04_tahoe_data/execution/build_pretrain_view.py` →
  `05_prnet_pretrain/execution/` (the view is consumer-owned).
- `04_tahoe_data/output/pseudobulk/tahoe_pseudobulk.h5ad` →
  `output/corpus/v2/tahoe.h5ad` (renamed only; 25,805,562,874 B, shape
  67,018 × 62,710).
- `04_tahoe_data/output/pseudobulk/{manifest,qa,qa_report,qa_contrasts}` →
  `execution/build_record/tahoe_pseudobulk_*`.
- `04_tahoe_data/output/pseudobulk/{_parts,_state}` → `execution/_parts/`.
- `04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad` and its four
  tracked side files → `05_prnet_pretrain/output/pretrain_view/` (unchanged,
  `content_sha256 = 476812c0…`).
- Deleted: `04_tahoe_data/output/pretrain_view/_view_parts/` (8.26 GB gather
  cache; rebuild is 175 s) and the empty `02_cm_data/` + `04_tahoe_data/` trees.

Move evidence (hashes, bytes, inodes, disk before/after, `_view_parts`
reclaimed, `--check` peak RSS) is in
`execution/build_record/tahoe_corpus_move.json` and
`execution/build_record/tahoe_cm_schema_move.json`.
