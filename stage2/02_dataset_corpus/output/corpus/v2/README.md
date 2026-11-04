# `corpus/v2` — Stage 2 dataset corpora, crosswalk tables and coverage

**Status:** generated artifact. Five h5ad corpora, seven crosswalk tables and this coverage record.
**Governing decision:** [ADR 0033](../../../../../docs/adr/0033-stage2-dataset-corpus.md) (the corpus layer; supersedes the matrix-free ADR 0026 shape and folds `04_tahoe_data` in).
**Built:** `2026-11-04` by `stage2/02_dataset_corpus/execution/build_corpus.py --tables`. The numbers below are also in [`execution/build_record/corpus_coverage.json`](../../../execution/build_record/corpus_coverage.json), so gate G2 can be checked programmatically.

Output is exactly three artifact classes (ADR 0033 decision 6): the five h5ad corpora, the seven tables under `tables/`, and this README. Build QA/manifests live under `execution/build_record/`, never here.

## Corpus coverage

One row per corpus. Five counts, computed from the corpus's own `obs` at its own time point: distinct cell line, chemicals among perturbed wells, distinct doses, time window, and distinct `(cell line, chemical, dose)` combinations. Vehicle controls and technical blanks are listed separately and are **never** counted as chemicals.

Counting rules (the same rules the JSON record carries):

- **cell line** = distinct `cell_line` labels; the primary arm labels its four lines `A/B/D/E` on plates 1/2/5 and `PC-A/PC-B/PC-D/PC-E` on plate 4, so those labels are normalised by removing the `PC-` prefix (recorded: `A, B, D, E, PC-A, PC-B, PC-D, PC-E` → 4 lines). `NA` counts as the one unnamed line of `sctox`/`tenv`.
- **chemical** = distinct `drug_key` over perturbed rows (Tahoe: its native `drug` column).
- **dose** = distinct `dose_um` over perturbed rows, as recorded (not rounded). `sctox` doses are `derived_top_dose` values from a curated per-drug top, so its per-drug 8-point series gives many distinct floats across the 42 drugs (115 values in total) while each drug has 8 doses.
- **time window** = the corpus's own recorded time point(s). `tahoe` has no time column (ADR 0029 update 2), so it is `n.a.`; the atlas is a single 24 h exposure by design (see the sub-row).

| corpus | distinct cell lines | chemicals (perturbed wells) | distinct doses (µM) | time window | (cell line, chemical, dose) | control wells | blank wells |
|---|---:|---:|---:|---|---:|---:|---:|
| `dtoxs_primary` | 4 | 47 | 21 | 48 h | 164 | 102 | 0 |
| `dtoxs_ipsc` | 2 | 37 | 14 | 48 h | 74 | 16 | 0 |
| `sctox` | 1 | 42 | 115 | 24 h | 336 | 552 | 0 |
| `tenv` | 1 | 460 | 3 | 24 h | 1,380 | 272 | 32 |
| `tahoe` | 50 | 379 | 3 | n.a. (no time column) | 56,827 | 1,400 | 0 |

### 24 h comparable subset

A corpus enters this sub-row when its exposure is a single 24 h window. The two `dtoxs` corpora are 48 h, so they are `n.a.`. `tahoe` is a 24 h atlas by design but its `obs` carries no time column (ADR 0029), so its row is the whole corpus and the basis says so — no time field is materialised.

| corpus | distinct cell lines | chemicals | distinct doses | (cell line, chemical, dose) | basis |
|---|---:|---:|---:|---:|---|
| `dtoxs_primary` | n.a. | n.a. | n.a. | n.a. | n.a. (48 h window) |
| `dtoxs_ipsc` | n.a. | n.a. | n.a. | n.a. | n.a. (48 h window) |
| `sctox` | 1 | 42 | 115 | 336 | recorded 24 h |
| `tenv` | 1 | 460 | 3 | 1,380 | recorded 24 h |
| `tahoe` | 50 | 379 | 3 | 56,827 | atlas design, no time column |

Corpus rows, controls and blanks (the contract table):

| corpus | rows | perturbed | control | blank | plates |
|---|---:|---:|---:|---:|---:|
| `dtoxs_primary` | 779 | 677 | 102 | 0 | 4 |
| `dtoxs_ipsc` | 312 | 296 | 16 | 0 | 1 |
| `sctox` | 3,072 | 2,520 | 552 | 0 | 8 |
| `tenv` | 4,510 | 4,206 | 272 | 32 | 16 |
| `tahoe` | 67,018 | 65,618 | 1,400 | 0 | 14 |

## Provenance

The accession/atlas name stays on the raw folder (ADR 0030); the short dataset id is the Stage 2 name. Raw sources are read-only under `pert_data/`; the build makes no network call.

| corpus | accession / atlas | raw item (`pert_data/<item>/source`) | shape (rows × genes) | built | bytes |
|---|---|---|---|---|---:|
| `dtoxs_primary` | LINCS-LDS-1593 — DToxS PromoCell primary CM panel | `pert_data/LINCS-LDS-1593/source` | 779 × 26,328 | 2026-11-04 | 47,402,239 |
| `dtoxs_ipsc` | GSE217423 — DToxS iPSC-CM panel | `pert_data/DToxS-iPSC/source` | 312 × 26,328 | 2026-11-04 | 18,902,113 |
| `sctox` | GSE244740 — ScreenSeq hiPSC-CM cardiotoxicity screen | `pert_data/GSE244740/source` | 3,072 × 31,880 | 2026-11-04 | 164,875,596 |
| `tenv` | GSE262419 — TempO-Seq hiPSC-CM environmental chemical panel | `pert_data/GSE262419/source` | 4,510 × 19,613 | 2026-11-04 | 255,148,955 |
| `tahoe` | Tahoe-100M — Tahoe-100M perturbation atlas (50 cancer cell lines) | `pert_data/Tahoe-100M/source` | 67,018 × 62,710 | 2026-10-30 (built by `04_tahoe_data`), moved unchanged 2026-11-04 (#14) | 25,805,562,874 |

Content and file hashes (equivalence is content-level, never file-byte — ADR 0033 decision 7; the four CM corpora use the ticket #15 digest recipe, `tahoe` the move hash of #14):

| corpus | `content_sha256` | file sha256 |
|---|---|---|
| `dtoxs_primary` | `ad1991974dd997219ddcdb777cd387440ce0730728fc57b308fece78ac6c1735` | `a29540f403cb407c3236e3999eb176a6cc5bc64ce6f85e8089a938eeed9d9d1c` |
| `dtoxs_ipsc` | `6c5e2f2f9806f76ce77401498ff496b01ec2f3aed6178e372be32d0862394f24` | `16e846212c9fc20b556141bbc6b4efbfe4f2355beb53403b4a9e54540b467a1c` |
| `sctox` | `fa8a399021fad60471195d77d8b45851ace67bb87604682c8db0d531da1e203d` | `c254d9096c59b7a25306591c40971066ab65e4cc9362e99ddcc3f517d810fe6b` |
| `tenv` | `81952fe42f1af1107a18892be6a3cbacf1a5b05ceed2feda1bec982b6499957f` | `a3d521f0ca0fbc732f2f6ba5e31ec2f1efb000102c30d64b7c991cc23ed4a886` |
| `tahoe` | `e591e157c3dd65cf23d338040a35d1d1ce2fbd0f3e5a2b3dd445549327a308a9` | `300955e0b7739b87960d1056536eb3c14e10f1376792c695a137eed818ed7f45` |

### Build environment

- Python 3.13.13, anndata 0.12.16, h5py 3.16.0, numpy 2.4.6, pandas 2.3.3
- platform: `Linux-6.18.54.1-microsoft-standard-WSL2-x86_64-with-glibc2.39`
- vendored Stage 1 loaders source commit: `c921f7d`
- pinned annotation assets: `pert_data/Annotation/source/hgnc_complete_set.txt` (HGNC) and `pert_data/Annotation/local/cm_compound_annotation.csv` (compound registry), both pinned by `pert_data/Annotation/local/README.md`
- the same `arc` environment builds the four CM corpora; `tahoe` predates this step (built 2026-10-30 by `04_tahoe_data`, moved unchanged by #14), and its QA numbers live in `execution/build_record/tahoe_pseudobulk_*`.

## Rename table (forward-only, ADR 0033 decision 9)

| historical id / path | current id / path | note |
|---|---|---|
| `screen_ctox` | `sctox` | GSE244740 ScreenSeq hiPSC-CM panel |
| `temp_env` | `tenv` | GSE262419 TempO-Seq hiPSC-CM panel |
| `dtoxs_cm` | `dtoxs_primary` + `dtoxs_ipsc` | the merged v1 view splits back into its two arms |
| `04_tahoe_data` | `02_dataset_corpus` | step folded in; `04` is a gap in the sequence, not a renumber |

Historical records keep the **old** names: a frozen 2026-10 log, manifest or issue that says `screen_ctox`, `temp_env`, `dtoxs_cm` or `04_tahoe_data` is not rewritten. Only new corpus artifacts use the new ids. Inside these tables the `dataset_id` column of `compounds.csv`, `source_gene_sets.csv` and `gene_crosswalk.csv` keeps the four **source** ids (`dtoxs_primary`, `dtoxs_ipsc`, `screen_ctox`, `temp_env`) because those tables are content-identical with the frozen v1 identities; `metadata_cm.csv` carries **both** the new `dataset_id` and the legacy id so a consumer never has to re-map the rename twice.

## Per-corpus `obs` and `var`

### Four CM corpora — `obs` 27 columns, index `well`

```
well, dataset_id, compound_name, drug_key
inchikey, smiles, dose_um, dose_source
time_h, is_control, vehicle, plate
control_group, cell_type, donor, replicate
atc, moa, mesh_classes, pharm_atc_moa
pharm_mesh, annotation_source, lib_size, n_genes_detected
cell_line, n_raw_rows, is_blank
```

`var` = index `symbol`; columns `hgnc_id, entrez_id, ensembl_id, match_type`. `X` is dense `int32` plus a dense `float32` `layers['log1p_cp10k']`; `uns`/`obsm`/`varm`/`obsp`/`varp` are empty; `dataset_id` and `control_group` use the v2 ids while `well` keeps the loader string verbatim (`screen_ctox::…`, `temp_env::…`).

### `tahoe` — `obs` 9 columns, index `profile_id`

```
plate, well, cell_line_id, n_cells, cell_name, drug, dose_um, is_control, control_group
```

`var` = index `symbol` (the atlas gene symbol); columns `ensembl_id, token_id`. `X` is CSR `int32` plus a CSR `float32` layer; one row is a `(plate, well, cell_line)` pseudobulk profile.

### Why the containers are asymmetric

The four CM corpora are small (312–4,510 rows) and share one producer, so they are dense, carry the same 27-column ADR-0026-derived `obs`, and index on the raw well id. `tahoe.h5ad` is a 25.8 GB, 67,018-row CSR artifact whose identity vocabulary is the atlas's: rewriting its `obs` into the CM shape would mutate the artifact without adding information (ADR 0033 decision 10). Its `obs` therefore stays the Tahoe-native nine columns with index `profile_id`, and `is_blank` is supplied as a constant `False` column in `metadata_tahoe.csv` rather than written into the h5ad. Format consistency therefore means the container and value semantics, not an identical column list; joins go through `genes.csv`, `gene_crosswalk.csv`, `compounds.csv`, `metadata_cm.csv` and `metadata_tahoe.csv`.

## Crosswalk tables

| table | rows | sha256 | purpose |
|---|---:|---|---|
| [`tables/genes.csv`](tables/genes.csv) | 18,164 | `19ce44acfd615a601d5bc6617074941e13d51cd6e527c8d519174e5085385f42` | the frozen 18,164-symbol CM projection axis (content-identical to v1) |
| [`tables/gene_crosswalk.csv`](tables/gene_crosswalk.csv) | 269,067 | `bc19c6a84bba8eba7e702b79b715b1032d05c3b47f1f44abd3d64c1ad24d8544` | every source native gene id → HGNC symbol (gate G3) |
| [`tables/compounds.csv`](tables/compounds.csv) | 532 | `dc6ff85322bd8b30a8d35b7e0d540593b881e18b2a76702785cd2a3bc9fc753f` | 525 `drug_key` / 532 key-name rows with `atc/moa/mesh_classes`, the three tier booleans, `datasets`, `n_wells` |
| [`tables/excluded_compounds.csv`](tables/excluded_compounds.csv) | 9 | `0250e0de423c6cd3d2fb0e41f9cc0e83b1bef755592f9c3e82125cfaa7741b24` | 8 non-SMILES perturbagen keys (143 wells) and the 32 technical blanks, with reasons |
| [`tables/source_gene_sets.csv`](tables/source_gene_sets.csv) | 104,149 | `b0152a381cf148747bcf452d64d64d1b304e3d4a0506b2af40cbea39938ff771` | per-source resolved gene sets (drift baseline) |
| [`tables/metadata_cm.csv`](tables/metadata_cm.csv) | 8,673 | `ad350023ace4d1497dbac21088304785cfc1b4586bad3715e7aba67992605886` | the four CM corpora's `obs` flattened, join-ready on `well` |
| [`tables/metadata_tahoe.csv`](tables/metadata_tahoe.csv) | 67,018 | `37e191d64eb0cfa296d48d10b62ea46e4f0839e10eafbf70540c09f5a9a7a1ac` | the Tahoe corpus's `obs` flattened, join-ready on `profile_id` |

`gene_crosswalk.csv` has **269,067 rows**: the 143,647-row four-source table (content-identical to the frozen `_frozen/cm-data-v1/datasets-1/gene_crosswalk.csv`) plus 125,420 Tahoe rows. The four CM sources each expose one native id space (DToxS symbols, `screen_ctox` Ensembl ids, `temp_env` `SYMBOL_probe` ids), so they contribute one row per raw gene id; Tahoe exposes two native namespaces per gene (`token_id` and `ensembl_id`), so each of its 62,710 genes gets one row per namespace. `in_corpus` marks whether the row resolves onto the 18,164 CM axis; `symbol=""`/`hgnc_id=""` marks an unmapped raw id (nothing is guessed).

Row counts: `genes.csv` 18,164; `compounds.csv` 532 key/name rows / 525 `drug_key`; `excluded_compounds.csv` 8 + 1 blank; `source_gene_sets.csv` 104,149 (the sum of the four native axes, 26,328 + 26,328 + 31,880 + 19,613); `metadata_cm.csv` 8,673 (779 + 312 + 3,072 + 4,510); `metadata_tahoe.csv` 67,018.

## Gate results

| gate | result | evidence |
|---|---|---|
| G1 structure + counts | **PASS** (#15) | per-corpus build gates passed: `dtoxs_primary` 12/12; `dtoxs_ipsc` 12/12; `sctox` 12/12; `tenv` 12/12 |
| G2 counts | **PASS** (this build) | this coverage table is rendered from `corpus_coverage.json`; `test_tables_contract` re-parses the table and asserts equality item by item |
| G3 crosswalk coverage | **PASS** (this build) | every `var` id of all five corpora has a crosswalk row (incl. Tahoe `token_id`/`ensembl_id`); the 18,164 `genes.csv` symbols are a subset of each CM corpus `var`; `compounds.csv` covers every non-control/non-blank `drug_key`. Verified before writing; see `test_tables_contract` |
| G4 v1 reconstruction | not run | ticket #17 (needs the frozen v1 view writer; see `.scratch/stage2-corpus-reorg/impl-notes-g4-proof.md`) |
| G5 reproducibility | not run | ticket #17 (ticket #15's `--repeat-check` already passes for `dtoxs_ipsc`) |
| G6 no stale paths | not run | ticket #17 |

## Regenerating the tables and this README

From the repository root (the repo root is found by walking up to `.git`, never by counting `../..` — ADR 0030 §4):

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
# the four CM corpora must already exist (ticket #15)
$PY stage2/02_dataset_corpus/execution/build_corpus.py --tables
$PY -m unittest stage2/02_dataset_corpus/execution/tests/test_tables_contract.py -v
```

The identity half re-reads the four raw sources to resolve every native gene id (the corpora keep only resolved genes), so it pays the same global load as a corpus build (~7.7 GiB peak, ~60 s); the row metadata half only reads each corpus's `obs`/`var` with `h5py` + `read_elem` (no matrices). The tables, this README and `corpus_coverage.json` are all overwritten atomically (`.tmp` → rename).
