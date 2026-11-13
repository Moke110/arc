# Data inventory

> **Measured with `du -sb`** (bytes; decimal units, 1 GB = 10^9 B). The ten top-level trees were measured 2026-10-30T23:01:54+0800; the `stage2/` subtree was re-measured 2026-11-04 after the 02+04 corpus merge, and again on 2026-11-07 after removal of the old benchmark and PRnet steps. The inventory totals in this note use the older measurements. Re-run `du -sb <path>` to refresh; small drifts come from regenerated caches/logs. The ten measured trees sum to **614,789,001,785 B** (614.8 GB); tracked docs plus small top-level files make the whole repository content ≈615 GB (excluding `.git` and the `.scratch/` workspace).
>
> This page is the human index of every large asset deliberately outside version control: where it came from, how big it is, whether it can be rebuilt, and whether it must be backed up. **Machine-readable source of truth** for the external data is each item's `author_log.json` (`sources[]`: origin, version, files, bytes, `tree_sha256`; `local_additions[]`) plus its `manifest.sha256`; `script/verify_all.sh` preserves the 27 legacy 項 and also validates all active items in `pheno_data/registry.json`; the 27/27 result below is the historical 2026-10-30 check. Contract: [ADR 0030](adr/0030-external-data-provenance-and-step-layout.md). A newcomer should read that ADR first, then a category `README.md`, then a step `README.md`.

## What is outside git, and why

Three classes of asset live outside version control:

1. **Downloaded content** — external items' `source/` trees (byte-identical upstream downloads) and declared `local/` additions. Gitignored by ADR 0030, hashed per file, verified offline.
2. **Generated artifacts** — `stage1/**/output` and `stage2/**/output` trees. Expensive (GPU-hours to days) but reconstructible by rerunning the tracked code; only the conclusion tables are force-tracked.
3. **External archives and non-code material** — `archived/`, `reference/`, `wet_data/`, `logs/`.

## Totals

| tree | `du -sb` bytes | ≈ | contents | in git? |
|---|---:|---:|---|---|
| `pert_data/` | 443,111,609,774 | 443.1 GB | 8 項 | source/ + local/ ignored; scripts, author logs, manifests tracked |
| `pert_model/` | 9,843,520,694 | 9.8 GB | 14 項 | source/ + local/ ignored; scripts, author logs, manifests tracked |
| `disease_data/` | 7,423,432 | 7.4 MB | 3 項 | source/ ignored; local/deg/ mixed (JSON tracked, CSV ignored) |
| `drug_emb_model/` | 930,113,535 | 930.1 MB | 2 項 | source/ tracked (code libraries); only `source/models/` + `feature_extract/` ignored |
| `stage1/` | 22,459,911,402 | 22.5 GB | 5 workflows (v0, v1, v1.1, v2, v3) | code + reports tracked (v0 gitignored); fold outputs ignored |
| `stage2/` | 120,327,454,763 | 120.3 GB | 3 steps + non-step dirs | docs/scripts/results tracked; output trees and venvs ignored |
| `archived/` | 3,593,246,074 | 3.6 GB | orphaned HG FDA extension | ignored (provenance only) |
| `reference/` | 1,635,164,057 | 1.6 GB | literature library | ignored |
| `wet_data/` | 296,458 | 296.5 kB | wet-lab raw data | ignored |
| `logs/` | 811,271 | 811.3 kB | runtime logs | ignored |

Notes: `stage2/**/output` carries the bulk of Stage 2 (the tenv benchmark payload and the corpus h5ads); `stage1/**/output` is dominated by fold-level predictions/checkpoints; `reference/PDFs` + `reference/_work` dominate the literature library. The 285 external weight/data symlinks were materialized into the repository in ticket 02 — the archive they came from (`/mnt/f/arc-weights`, 40 GB) remains untouched as an offline backup.

**Hard-link caveat.** `du -sb` counts a hard-linked inode **once per invocation**. There are 15 hard-linked files across `stage1/xpert-fda-miri-v{1,1.1,3}/input/` (notably the 5.57 GB `all_drugs_unimol_arr.npy`), and the venv trees share many inodes internally. Every row below was measured with its **own** `du -sb <path>` invocation, so per-subtree rows match a single-path `du`; sums of rows can exceed a parent total (the `stage1/` total counts the shared UniMol array once, while three workflow rows each count their own copy).

## External data — legacy 27 項 plus registered phenotype sources

Each external item is one directory per download source under ADR 0030: `script/` (downloader + offline verifier), `source/` (untouched download), `local/` (declared local additions), `author_log.json`, `manifest.sha256`. The sizes below are the whole item directory (`du -sb`); the source-only file counts and bytes are in the item's `author_log.json`, and the per-file hashes in its `manifest.sha256`. All 27 verified 27/27 on 2026-10-30.

### pert_data — 8 項 (443.1 GB)

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `pert_data/Tahoe-100M/` | 428,812,221,369 | 428.8 GB | yes — HF `tahoebio/Tahoe-100M` rev `349e2a55…` (CC0-1.0) | re-download 428.8 GB; offline backup of the imported copy is `/mnt/f/arc-weights` (kept outside git) |
| `pert_data/LINCS-LDS-1593/` | 1,669,342,170 | 1.7 GB | yes — clue.io **LDS-1593** release (manual browser step documented in `script/download.md`) | re-download; `local/v1-frozen/` is the only frozen v1-era copy → back up `local/` |
| `pert_data/DToxS-iPSC/` | 2,657,116,588 | 2.7 GB | yes — GEO SuperSeries **GSE217424** (GSE217421 + GSE217423) | re-download from GEO |
| `pert_data/GSE244740/` | 37,483,465 | 37.5 MB | yes — GEO **GSE244740** | re-download from GEO |
| `pert_data/GSE247541/` | 1,002,664 | 1.0 MB | yes — GEO **GSE247541** | re-download from GEO |
| `pert_data/GSE262419/` | 83,449,159 | 83.4 MB | yes — GEO **GSE262419** | re-download from GEO |
| `pert_data/XPert-L1000/` | 9,833,362,269 | 9.8 GB | **NO — author-provided, upstream version unknown** | **back up**: `source/` (9.8 GB) has no recorded upstream URL |
| `pert_data/Annotation/` | 17,628,280 | 17.6 MB | partly — the HGNC complete set is re-downloadable (pinned 2026-10-21); the compound snapshot + `cm_drugs.html` are local | back up `local/` |


### pert_model — 14 項 (9.8 GB)

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `pert_model/PRnet/` | 109,154,121 | 109.2 MB | yes — upstream repo `f19174bd` + authors' checkpoints | re-download; weights also mirrored in `/mnt/f/arc-weights` backup |
| `pert_model/CPA/` | 1,534,463,110 | 1.5 GB | yes — upstream repo `fbd7c025` + Kang/Norman2019 checkpoints | re-download; mirror in `/mnt/f/arc-weights` |
| `pert_model/chemCPA/` | 1,059,126,707 | 1.1 GB | yes — upstream repo `43e830eb` + public ComPert weights | re-download; mirror in `/mnt/f/arc-weights` |
| `pert_model/biolord/` | 180,829,942 | 180.8 MB | yes — upstream repo `b7688790` (no public weights) | re-download |
| `pert_model/PerturbNet/` | 2,206,330,215 | 2.2 GB | yes — upstream repo `e1375b2c` + released weights | re-download; mirror in `/mnt/f/arc-weights` |
| `pert_model/scVIDR/` | 165,201,892 | 165.2 MB | yes — upstream repo `bcb6a96c` + `data.zip` | re-download; mirror in `/mnt/f/arc-weights` |
| `pert_model/Squidiff/` | 1,240,203 | 1.2 MB | yes — upstream repo `abdfc27d` + released checkpoints | re-download |
| `pert_model/CMonge/` | 34,042,646 | 34.0 MB | yes — upstream repo `195f5cd1` + bundled demo data | re-download |
| `pert_model/DECANT/` | 561,086 | 561.1 kB | yes — upstream repo `e50642a9` (no shipped data/weights) | re-download |
| `pert_model/WAVE/` | 1,530,226 | 1.5 MB | yes — upstream repo `e27e453b` + FigShare assets | re-download |
| `pert_model/ExPO/` | 223,097 | 223.1 kB | yes — upstream repo `7e88cbae` | re-download |
| `pert_model/State/` | 4,356,261,848 | 4.4 GB | yes — upstream repo `9bbfe78a` + non-commercial weights | re-download; mirror in `/mnt/f/arc-weights` |
| `pert_model/TranSiGen/` | 67,651,615 | 67.7 MB | yes — upstream repo `e7582e77` (nested `.git` removed; commit recorded) | re-download |
| `pert_model/XPert/` | 126,896,363 | 126.9 MB | **NO — author-provided, version unknown**; `local/` holds the locally authored routing index/logs/metadata/metrics | **back up** `source/` + `local/` |


### disease_data — 3 項 (7.4 MB)

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `disease_data/GSE108940/` | 2,787,159 | 2.8 MB | yes — GEO **GSE108940** (rat, Agilent GE2-v5) | re-download; `local/deg/` DEG lists are derived → back up |
| `disease_data/GSE137293/` | 2,601,105 | 2.6 MB | yes — GEO **GSE137293** (pig, NovaSeq) | re-download; `local/deg/` derived |
| `disease_data/GSE160516/` | 2,030,663 | 2.0 MB | yes — GEO **GSE160516** (mouse, Clariom S) | re-download; `local/deg/` derived |


### disease_data — Stage 2 CM sources (#187, ADR 0066)

Three cardiomyocyte-derived items added for the Stage 2 gene axis. Each is an
ADR 0030 item (tracked `script/`/`author_log.json`/`manifest.sha256`, ignored
`source/` + declared `local/`); the CM-gated per-donor pseudobulk and the
item-local DEG stay in `local/`, and the compact published product is the step's
`output/disease-deg/`. `du -sb`, 2026-11-13.

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `disease_data/GSE270788/` | 951,113,459 | 951.1 MB | yes — GEO **GSE270788** GEX SubSeries: 23 per-GSM `*_filtered_feature_bc_matrix.h5` (931,875,148 B) + cell metadata (1,272,635 B); `script/fetch.sh`. `local/` pseudobulk + DEG are derived | re-download `source/`; tracked summary/manifest in git; `local/` derived |
| `disease_data/GSE255612/` | 1,032,251,016 | 1.03 GB | yes — GEO **GSE255612**: Scanpy mtx (1,003,414,933 B), feature axis (277,568 B) and cell metadata (8,846,174 B); `script/fetch.sh` (parallel ranges). `local/` pseudobulk + DEG derived | re-download `source/`; tracked summary/manifest in git |
| `disease_data/Reichart-CM/` | 3,201,296,579 | 3.2 GB | yes — CELLxGENE object, collection `e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb`, version `27fdbdf9-deae-4aed-bf7f-21964ee62137` (CC BY 4.0), 3,194,863,189 B; `script/fetch.sh` | re-download; the **tracked** axis-aligned per-donor pseudobulk + genotype export (ticket #188 input) is in git |

### disease_data — Stage 2 isogenic CM models (#188, ADR 0066)

Six in-vitro hiPSC-CM bulk items added for the `genetic_cardiomyopathy`
category. Each stores two byte-size-pinned payloads (the deposited expression
matrix + the GEO series matrix) in `source/`; the item-local expression matrix,
design sheet and full DEG table are derived under `local/deg/` with a tracked
`_deg.json` summary, and the compact published product is the step's
`output/disease-deg/genetic_cardiomyopathy/`. `du -sb`, 2026-11-13 (whole item
directory, including the derived `local/` payloads).

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `disease_data/GSE224129/` | 11,849,846 | 11.8 MB | yes — GEO **GSE224129** (GPL24676): RSEM gene matrix 826,308 B + series matrix 3,834 B; `script/fetch.sh` | re-download `source/`; tracked summary/manifest in git; `local/` derived |
| `disease_data/GSE234141/` | 22,693,331 | 22.7 MB | yes — GEO **GSE234141** (GPL28337): counts 1,123,421 B + series matrix 3,460 B; `script/fetch.sh` | re-download; `local/` derived |
| `disease_data/GSE304575/` | 6,430,249 | 6.4 MB | yes — GEO **GSE304575** (GPL24676): raw-count table 817,509 B + series matrix 3,382 B; `script/fetch.sh` | re-download; `local/` derived |
| `disease_data/GSE305648/` | 19,583,646 | 19.6 MB | yes — GEO **GSE305648** (GPL30173): counts 765,978 B + series matrix 5,479 B; `script/fetch.sh` | re-download; `local/` derived |
| `disease_data/GSE253223/` | 6,699,341 | 6.7 MB | yes — GEO **GSE253223** (GPL24676): estimated counts 1,526,983 B + series matrix 4,219 B; `script/fetch.sh` | re-download; `local/` derived |
| `disease_data/GSE194297/` | 4,528,189 | 4.5 MB | yes — GEO **GSE194297** (GPL20301): absolute counts 342,625 B + series matrix 2,746 B; `script/fetch.sh` | re-download; `local/` derived |

`disease_data/Reichart-CM/` gains one additive derivation for ticket #188
(`script/build_genetic_deg.py` → ignored `local/deg/reichart-cm_genetic_deg.csv`
plus a tracked `_deg.json` summary, ~4 MB on disk); its `source/` and #187 bytes
are unchanged.

### drug_emb_model — 2 項 (930.1 MB)

Unlike the other three categories these are **code libraries**: their `source/` trees are tracked; only the downloaded assets (`KPGT/source/models/`, `KPGT/source/feature_extract/`) are ignored.

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `drug_emb_model/HG/` | 141,504 | 141.5 kB | **NO — author-provided** (ADR 0006), version `unknown`; no upstream commit | **back up** `source/` (only copy) |
| `drug_emb_model/KPGT/` | 929,967,615 | 930.0 MB | **NO — vendored, URL/commit `unknown`**; `source/models/` (894 MB weights) + `source/feature_extract/` (34 MB) have no pinned origin | **back up** `source/`; weights also mirrored in `/mnt/f/arc-weights` |

### Not reconstructible — back these up

- `pert_data/XPert-L1000/source/` (9.8 GB), `drug_emb_model/KPGT/source/` (930 MB, mostly model weights) and `drug_emb_model/HG/source/` (author-provided code): all three carry `version: "unknown"` / no upstream commit, so the copy in the repository is the only record.
- Every declared `local/` addition (most are small, but they include the frozen v1-era LINCS package `pert_data/LINCS-LDS-1593/local/v1-frozen/`, the v0-era gene sets, and the annotation snapshot) — see `local_additions[]` in each `author_log.json`.
- `stage1/xpert-fda-miri-v0/` — frozen legacy snapshot whose agent-generated reports cannot be regenerated byte-identically.
- `archived/hg-fda-extension-archive/` and `wet_data/` — no upstream and no synthetic substitute.

## Stage 1 — 22.5 GB

The four finalised screens are rerunnable from the 項 sources with the shared `arc` conda environment; their outputs are ignored. `xpert-fda-miri-v0/` is the frozen legacy snapshot (gitignored; only the canonical reports from its post-processing survive, having been de-duplicated in ticket 10).

| tree | `du -sb` bytes | ≈ | what | reconstructible? | backup |
|---|---:|---:|---|---|---|
| `stage1/xpert-fda-miri-v0/` | 116,820,617 | 116.8 MB | frozen legacy snapshot (v0) | git history `47e1604^` + rerun is approximate only; the on-disk snapshot is the only copy of its agent-generated reports | **back up** (116.8 MB): frozen evidence |
| `stage1/xpert-fda-miri-v1/` | 7,093,994,482 | 7.1 GB | canonical 11-step screen | yes — rerun the tracked workflow on the 項 sources (`arc` conda env); ~hours | no backup needed; rerun |
| `stage1/xpert-fda-miri-v1.1/` | 7,130,609,600 | 7.1 GB | DToxS-LINCS testbed | yes — rerun | no backup needed; rerun |
| `stage1/xpert-fda-miri-v2/` | 9,330,527,774 | 9.3 GB | primary + iPSC-CM | yes — rerun | no backup needed; rerun |
| `stage1/xpert-fda-miri-v3/` | 9,949,375,762 | 9.9 GB | all-CM | yes — rerun | no backup needed; rerun |

Per-workflow rows are isolated `du -sb` values that each include their own hard-linked `input/` copies; the `stage1/` total above counts the shared 5.57 GB UniMol array once, which is why the five rows sum to more than it.

## Stage 2 — the two active steps (2026-11-12, ADR 0062 stable layout)

`stage2/01_bases/` (the foundation) and `stage2/02_models_tenv/` (the tenv model
comparison restored by ADR 0063, tickets #163-#165) are the active steps. Other
Stage 2 directories, including ignored checkpoints and historical predictions,
were deleted at the user's request. Tracked historical code and records remain in
Git; the retired comparison implementation is reachable only through Git history
(`git show 1543781^:stage2/03_modes_tenv/...`) and is untouched by the restored
step. Recovering deleted ignored model fits would require retraining; no current
build depends on them.

Every accepted product now has one stable semantic path (`du -sb`; the retained
rows below were measured 2026-11-12/13, historical sizes for removed payloads are
kept only in git history). #157 removed 15.27 GB of rejected and disposable
payload — the rejected historical v1 Tahoe h5ad (11.20 GB) and two superseded
build-state directories (1.21 GB + 2.83 GB) — after its dependency proof. The
ADR 0066 rebuild (#186-#192) then **deleted the pre-0066 products**
`shared-expression/`, `training-corpora/`, `tenv-raw-counts/`,
`tahoe-raw-counts/`, `tenv-truth/`, `tenv-truth-successor/`,
`historical-training-release/` and the old `gene-axis/` file set, only after the
#192 acceptance replays passed; the old→new mapping and digest evidence are in
`stage2/01_bases/log/acceptance.md` + `inputs.json`. Git retains the compact
files and all code; the ignored multi-GB payloads are re-creatable from the
immutable `source/` trees below.

| stable path under `stage2/01_bases/` | `du -sb` | contents | reconstruction |
|---|---:|---|---|
| `output/gene-identity/` | 16,638,438 | indexed identity (`genes.csv`, `manifest.json`) and the frozen `tables/` axis/crosswalk | compact source gene metadata, pinned HGNC annotation and the historical v1 mapping evidence |
| `output/tenv-hvg/` | 6,243,580 | retained pre-0066 product: descriptive ranking (the tenv-HVG axis source) plus three train-only HVG rankings and deterministic parent-hash fold assignments | `build_hvg.py` — its pre-0066 input `training-corpora/tenv.h5ad` is retired; the same 4,206-well population is available from `output/corpus/tenv.h5ad` under the frozen normalization and the re-point is deferred (README *Retained tenv-HVG ranking*) |
| `output/fold_n_seed/` | 1,159,437 | the accepted cold-drug assignment re-homed (#161) plus the ordered `folds.csv`/`seeds.csv` index tables and manifest | `build_fold_n_seed.py` from the sealed HVG cohort assignment; row-prefix selection only |
| `output/palette/` | 3,204 | the stage palette: six hue families x nine shades as a hex table, its generation recipe and the colour roles (#162) | `output/api/bases/palette.py` `publish()`; the table is a pure function of the recorded recipe |
| `output/api/` | 512,471 | the single public reader/metric implementation imported by build/audit/test code | tracked code, not payload (see `.gitignore` exception and #156 D1) |
| `figure/` | 711,136 | `gene_axis_rank_sum_histogram.png`, its complete plotted-values CSV and the CSV-only renderer (the diagnostic PR figure was deleted by #192) | regeneration reads only `figure/data/*.csv` with `output/` products unreadable |
| `log/` | ~4.1 MB | the post-0066 acceptance record `acceptance.md` + `inputs.json`, the tracked `raw_counts_provenance.json` (the two pre-0066 source-replay records restored from git history) and the per-slice `*_verify.json` replay records (the rest of the pre-0066 tree was deleted, not archived — ADR 0066) | tracked evidence |

### Rebuild foundation products (#186, ADR 0066)

The rebuild's `raw_counts/` -> `corpus/` chain is accepted (#192:
`stage2/01_bases/log/acceptance.md`); the superseded raw, shared-expression and
training products are deleted. Both new areas are gitignored payloads with
tracked manifests, mapping evidence and exclusion table; `du -sb`, 2026-11-13.

| stable path under `stage2/01_bases/` | `du -sb` | contents | reconstruction |
|---|---:|---|---|
| `output/raw_counts/` | 3,037,082,960 | `tenv.h5ad` (4,510 wells x 19,613 native genes, int32), `tenv_probe_mapping.csv`, `tenv_axis_mapping.csv` (frozen `gene_index` -> native column), `tahoe.h5ad` (67,018 profiles x 18,164 axis genes, CSR int32), `token_mapping.csv`, `manifest.json` | `execution/build_raw_counts.py publish`: rehomes the source-accepted payloads and records `source-accepted-replay` provenance against the tracked `log/raw_counts_provenance.json` (the two pre-0066 source-replay acceptance records restored verbatim from git history plus the native tenv build manifest, digest-pinned). `--full-rebuild` rebuilds both from `pert_data/GSE262419/source` + `pert_data/Tahoe-100M/source` and replays from source (recorded ~4 h under the 10 GiB cap) |
| `output/corpus/` | 1,395,893,021 | `tenv.h5ad` (4,206 treated + 272 vehicle), `tahoe_bulk.h5ad` (65,318 treated + 1,400 `DMSO_TF`), `exclusions.csv` (32 blanks + 300 invalid structures), `manifest.json` with the matched-control report; int32, no normalization, no delta | `execution/build_corpus.py publish` from `output/raw_counts/` + `pert_data/Tahoe-100M/source/metadata/drug_metadata.parquet` |
| `output/disease-deg/` | 63,635,080 | 63.6 MB | per category `nhf/`, `icm/`, `af/`, `genetic_cardiomyopathy/`: one rank-sorted compact `<dataset>.deg.csv` + full-axis `<dataset>.rank.csv` per dataset, plus `catalog.json` (canonical record list, `catalog_sha256`) and `manifest.json`; no top-level catalog (ticket #190 aggregates the per-category files). Thirteen datasets: GSE270788 NICM/ICM/AMI, Reichart genotype-negative DCM, Reichart genotype-positive DCM + ARVC, GSE255612 AF, and six isogenic models (GSE224129/234141/304575/305648/253223/194297); `genetic_cardiomyopathy/` is 29,961,789 B of it | `execution/build_disease_deg.py publish --all` from the pinned item-local DEG tables in `disease_data/{GSE270788,GSE255612,Reichart-CM,GSE224129,GSE234141,GSE304575,GSE305648,GSE253223,GSE194297}/local/deg/`; `verify --all` replays and compares every digest. The compact CSVs/catalogs/manifests are force-tracked; the item-local DEG/pseudobulk payloads stay ignored |

Both products are replayed by their drivers' `verify` action (rebuild in
`execution/checkpoints/` and compare every digest). The compact manifests,
mapping evidence and exclusion table are force-tracked; the count payloads are
gitignored.

### Drug-cardiotoxicity class products (#189, ADR 0066)

The `cardiotox_TOP2` / `cardiotox_other` class DEGs are **tracked compact
tables**, not ignored payloads: each per-`(source, class)` DEG CSV is 1.9-2.8 MB
and its full-axis rank CSV 0.38-0.43 MB, 15.7 MB in total across the three
product directories. They are force-added because the repository-wide `*.csv`
ignore rule would otherwise drop them.

| stable path under `stage2/01_bases/` | `du -sb` | contents | reconstruction |
|---|---:|---|---|
| `output/disease-deg/cardiotox/` | 57,716 | pinned crosswalk, classes, class members, excluded codes, per-source eligibility, rejected candidates and manifest | `execution/build_disease_deg_cardiotox.py publish` from the read-only GSE243674 / DToxS-iPSC / GSE244740 sources + `pert_data/Annotation/local/cm_compound_annotation.csv` |
| `output/disease-deg/cardiotox_TOP2/` | 8,826,044 | three pooled class DEG tables + full-axis ranks, `catalog.json`, `manifest.json` | same driver; DToxS pinned edgeR pipeline, GSE243674 DESeq2 and the GSE244740 Welch test over pooled class wells |
| `output/disease-deg/cardiotox_other/` | 6,795,074 | three DToxS non-TOP2 class DEG tables + full-axis ranks, `catalog.json`, `manifest.json` | same driver, DToxS pinned edgeR pipeline over pooled class wells |

There is no gitignored cardiotox payload: the count/R intermediates exist only
inside `execution/checkpoints/` while a candidate is built and are deleted after
the candidate replaces the accepted product. `verify` replays the whole build
from the pinned inputs and compares every digest; the peak RSS of the build
(≈1.5 GiB Python, ≈0.9 GiB R) is recorded in the ticket record, not in a manifest.

### Gene-axis release (#190, ADR 0066)

The two released axes are **tracked compact evidence**, not ignored payloads:
`candidate_rank_sum.csv` is 1.35 MB, the two axis files 23 kB + 58 kB, the
composition table 5 kB and the manifest 25 kB; the figure CSV and PNG are 6 kB +
83 kB. They are force-added because the repository-wide `*.csv` ignore rule would
otherwise drop them. The pre-0066 `gene-axis/` file set (direction composition,
unweighted quota-free fused table, animal half-weight column, tenv non-HVG curve,
ortholog table, `inflection.json`, ~3.54 MB) was **deleted by #192** after
acceptance; nothing read it. `du -sb`, 2026-11-13.

| stable path under `stage2/01_bases/` | `du -sb` | contents | reconstruction |
|---|---:|---|---|
| `output/gene-axis/` | 1,466,000 | released products only (`candidate_rank_sum.csv`, `gene_axis_2000.csv`, `gene_axis_5000.csv`, `axis_composition.csv`, `manifest.json`, ~1.46 MB) | `execution/build_gene_axis.py publish` from the read-only `output/disease-deg/**`, `output/tenv-hvg/descriptive.csv` and `output/gene-identity/genes.csv`; `--verify` replays and compares every digest |
| `figure/gene_axis_rank_sum_histogram.png` | 82,830 | rank-sum histogram over all 18,164 candidates with both axes' `rank_sum` intervals marked (sha256 `49d1dc5d…`) | `figure/script/render_gene_axis_rank_sum_histogram.py` from `figure/data/gene_axis_rank_sum_histogram.csv` alone (`output/` unreadable) |
| `figure/data/gene_axis_rank_sum_histogram.csv` | 6,483 | complete plotted values: histogram bins/shares, both axis intervals/medians and the inflection diagnostic marker | written by the build from the pinned inputs; tracked next to the PNG |

The catalog inputs are the per-category `output/disease-deg/` products; the
manifest pins 33 input digests (frozen axis, tenv-HVG ranking, six catalogs,
six release manifests, 19 used rank frames) and 5 output digests, with no
absolute paths. Build + replay peak RSS ≈140 MiB under the 10 GiB cap.

### Model-facing views `view-2k/` / `view-5k/` (#191, ADR 0066)

The two views are the models' `X`/`control`/`Δ` products: each corpus profile is
normalized once on the **full** 18,164-gene library total (CP10K, natural
`log1p`), then subset to one released axis's ascending `gene_index`; `control`
is the mean of the row's matched vehicle rows normalized the same way and
`Δ = X − control` exactly. The view never re-normalizes on its own subset. The
float32 payloads are gitignored (normalized log values are near-incompressible)
and only the small manifests plus the reader are tracked; `du -sb`, 2026-11-13.

| stable path under `stage2/01_bases/` | `du -sb` | contents | reconstruction |
|---|---:|---|---|
| `output/view-2k/` | 1,371,438,657 | `tenv.h5ad` (4,478 × 2,000 float32: `X`, `layers/control`, `layers/delta`), `tahoe_bulk.h5ad` (66,718 × 2,000), `library_totals.csv` (71,196 profiles' full-axis totals, 3.7 MB), `manifest.json` (5,752 B) | `execution/build_view.py publish` from the read-only `output/corpus/` + `output/gene-axis/`; `verify` replays both views and compares every digest |
| `output/view-5k/` | 3,292,322,079 | `tenv.h5ad` (4,478 × 5,000), `tahoe_bulk.h5ad` (66,718 × 5,000), `library_totals.csv`, `manifest.json` (5,757 B) | same driver; candidates stage in `execution/checkpoints/` and are replaced only after verification |

Measured cost of one publish run: 3 m 35 s, peak RSS 445,512 KiB; one `verify`
replay: 3 m 29 s, peak RSS 443,036 KiB, both streaming 256-row corpus batches
under the 10 GiB cap. The measured full-axis CP10K basis is recorded in each
manifest: `tenv` **12.9299987%** (`G_2000`) / **37.2063189%** (`G_5000`) pooled over all
rows (treated-only 12.9566820% / 37.1590893%), `tahoe_bulk` 10.4600758% /
24.7819987%. ADR 0066's
design-record 12.9% / 36.3% came from an interim pre-#190 axis; the measured
released-axis value is authoritative and the one-normalization decision is
unaffected (≈2.88× denominator ratio). Small-scale reproduction:

```bash
python stage2/01_bases/execution/build_view.py publish   # both views, candidate -> verify -> replace
python stage2/01_bases/execution/build_view.py verify    # replay both, compare every digest
```

The two `manifest.json` files are force-tracked (they declare the payload
digests, input digests and the measured basis); the h5ad payloads and
`library_totals.csv` are ignored and reconstructible by the driver.

Large h5ads, arrays and checkpoint shards stay ignored; compact input tables,
the tracked `output/api` implementation and the records remain tracked. See
`stage2/01_bases/README.md` for the current commands, accepted results,
limitations and historical-versus-current provenance boundary, and
`stage2/01_bases/log/acceptance.md` + `inputs.json` for the post-0066 old→new
binding and the recorded errata (GSE270788 roster; released-axis view basis).

### `stage2/02_models_tenv/` — accepted tenv matrix (2026-11-12, #165/#170)

The step's products are untracked ignored payloads under its `output/`; only the
three figures, their three complete plotted-value CSVs, the CSV-only renderer,
the acceptance records and the raw resource measurements are tracked. The whole
`output/` tree is regenerable from the canonical checkout with
`python execution/run_matrix.py --all`; a completed run set is reused with
`--resume`, and `--accept-only` rescores and replays the published runs without
training or re-inference.

| path under `stage2/02_models_tenv/` | `du -sb` | contents | reconstruction |
|---|---:|---|---|
| `output/inputs/` | 2,041,030,268 | four complete predictor-input releases (`drugs.csv`, `morgan.npy`, per-role baseline/delta arrays, case tables and manifests), one per `(fold, seed)` matrix cell | `execution/script/prepare_inputs.py` from the sealed `01_bases` corpus and fold/seed products |
| `output/runs/` | 967,675,021 | eight run directories: four PRnet fits (`prediction.npy`, `best.pt`, checkpoint/identity/terminal/resource records) and four deterministic `Mean` fold copies | PRnet runs: pinned PRnet environment + `execution/script/prnet_fit_predict.py`; `Mean`: `execution/script/mean_reference.py`; one heavy job at a time under the 10 GiB cap |
| `output/evaluation/` | 57,435,713 | `per_well.csv`, `per_group.csv`, `per_run.csv`, `per_run_curve.csv`, `confusion.csv`, `evaluation.json` and the three plot-value copies from the shared #168 driver | `execution/script/evaluation_driver.py` over the eight complete runs; re-scores from persisted predictions only |
| `output/figure_data/` | 825,190 | copies of the three @d plot-value tables (byte-identical to the tracked `figure/data/`) | `evaluation_driver.py --figure-data` |
| `figure/` | 1,270,546 | `metrics.png` (4 direct panels), `metrics_indirect.png` (8 reliability panels), `learning_curve.png` (PRnet), the three complete plotted-value CSVs and the CSV-only renderer | `python figure/script/render_figures.py --figure-data figure/data --out-dir figure`; `output/` unreadable |
| `log/` | 89,486 | #163 admission, #164 driver, #165 and #166 acceptance records plus the raw `/usr/bin/time -v` measurements (including the three tracked `166-*.time` replay files) and the environment locks | execution records; tracked |

Historical tenv comparison code, figures and records are **not** present in the
working tree; they remain untouched in Git history (`git show
1543781^:stage2/03_modes_tenv/...`) and nothing in this step depends on them.

## Archived, literature, wet-lab, logs

| tree | `du -sb` bytes | ≈ | what | reconstructible? | backup |
|---|---:|---:|---|---|---|
| `archived/` | 3,593,246,074 | 3.6 GB | orphaned HG FDA-extension data/weights (ADR 0006; provenance only) | **NO** — author-given artifact with no public upstream | **back up** (only copy); originally materialized from the external archive |
| `reference/` | 1,635,164,057 | 1.6 GB | background literature library: PDFs (162), markdown transcripts (163), index, explainer pages, `_work/` build scripts | partly — PDFs are re-downloadable from publishers (DOIs in `references.xlsx`), but access can be closed | optional backup; not pipeline-critical |
| `wet_data/` | 296,458 | 296.5 kB | Stage-1 wet-lab raw data (result-figure pptx + animal-experiment plan) | **NO** — experimental raw data | **back up**: record in `docs/wetlab-validation-stage1.md` |
| `logs/` | 811,271 | 811.3 kB | runtime logs (Tahoe fetch/QA history; no pipeline artifact) | n/a — ephemeral | no backup needed |

## Verification and maintenance

- `script/verify_all.sh` — offline `sha256sum -c` over every item manifest; 27/27 pass (idempotent).
- `script/check_paths.sh` — every `pert_data/`/`pert_model/`/`disease_data/`/`drug_emb_model/` path literal in tracked text resolves.
- `script/check_tracked.sh` — no executable file is silently untracked.
- This page's byte counts are `du -sb` snapshots; refresh them when an item is re-downloaded or an output tree is rebuilt. Per-item provenance lives in the items, not here.

## Appendix — orphan dataset table

The acquired CM/screen data sources enter the repository as provenance items
(`pert_data/LINCS-LDS-1593/`, `pert_data/DToxS-iPSC/`, `pert_data/GSE244740/`,
`pert_data/GSE247541/`, `pert_data/GSE262419/`, `pert_data/XPert-L1000/`). The table below
is the orphaned root-level candidate list of external cardiomyocyte hypoxia/IR compound
datasets (compiled 2026-10-14); it is kept here verbatim so the content survives rather
than living as an untracked, one-`git clean`-away file. The original CSV is tracked as
[`tableA_cardiomyocyte_hypoxia_compound_datasets.csv`](../tableA_cardiomyocyte_hypoxia_compound_datasets.csv).

| 数据集 | 造模 | 化合物 | 表型 | 规模 | 获取 |
|---|---|---|---|---|---|
| Heliste 2021 Sci Rep | ✅ HL-1 小鼠心肌细胞(缺氧48 h + 复氧24 h) | ✅ 689 FDA 药物(40/200/1000/5000 nM) | CTG 细胞活力(H/R 下活力降 22%) | 689 化合物 × 4 浓度 × (H/R + 常氧对照) | CC-BY；论文补充材料 S1/S3，未结构化 deposit；https://www.nature.com/articles/s41598-021-96033-z |
| Guo 2012 Circ Res | ✅ H9c2 Seahorse IR 模型 | ✅ 2000 分子 | IR 后细胞死亡评分 | 2000 化合物(单点) | 论文补充材料；Circ Res 110:948-957 |
| Wang 2025 Acta Pharm Sin B | ✅ 3D 心肌组织(CM+CF+BMDM) 缺氧/OGD | ✅ 76 化合物(T. vernicifluum) | 收缩性(可量化) + 抗纤维化 | 76 化合物 | 论文 + GEO GSE275028 |
| Oroxin A 2022 Front Pharmacol | ✅ H9c2 OGD/R | ✅ 2661 天然化合物(10 µM) | CCK-8 活力 | 2661 化合物(单点) | CC-BY；doi 10.3389/fphar.2022.894899 |
| ChEMBL H9c2 H/R 系列 | ✅ H/R | ✅ 每个 AID 1–20 化合物 | 存活/凋亡/坏死、caspase-3、线粒体 ROS、Zn²⁺ | 极小 | PubChem AID 1780097–1780110、2109509–2109554 等(ChEMBL→PubChem) |

### pheno_data — phenotype originals (2026-11-04)

The fifth category follows ADR 0038. Its authoritative inventory is
`pheno_data/registry.json`, which distinguishes active acquired products from
candidate/unavailable products. Active products must have nonempty source-only
manifests and working offline verification; root inventory verification rejects
missing or unregistered items. The 27 legacy items remain mandatory.

`pheno_data/Sirenko-2013/` archives publisher Table 3 XLS, article XML and PMC
metadata from official URLs, with exact byte counts, hashes, media types and
retrieval timestamps in its author log. `source/` and declared `local/` additions
stay ignored; scripts and provenance stay tracked. Copyright is reserved: public
retrieval is not a redistribution license. These originals support the Stage 3
24 h beat-count adapter; numerical table coverage does not establish label or
training eligibility. Run `bash pheno_data/Sirenko-2013/script/verify.sh` offline.

The remaining active phenotype items are `pheno_data/Heliste-2021/`,
`pheno_data/Linn-2024/`, `pheno_data/Liu-2021/`, `pheno_data/Duan-2013/` and
`pheno_data/Huang-2022/`. The phenotype registry determines the mandatory acquired
items; Huang contains public evidence, with its full numerical screen unavailable.
Original images/tables remain ignored; small declared curation snapshots and
download/verification code provide reconstruction evidence. The immutable local
Heliste source-first corpora live in `stage3/01_corpus/output/baseline/` and
`stage3/01_corpus/output/H-R/`; chemical-resolution reports describe admitted and
pending values separately from source availability. Historical endpoint-gated
evidence is retained under `stage3/01_corpus/execution/historical/`. No model view or training result is implied by the archive.

## Shared Stage 2 tenv significance labels

Measured truth is owned by `stage2/01_bases/`, separately from training h5ads.
The current parent-normalized design has 1377 contrasts and a BH family of
exactly 18164 genes, with tested-empty distinct from untestable. Raw well joins,
contrast membership, worker/runtime seals and truth arrays must pass the
independent truth verifier before acceptance. The deleted historical 1380-group,
18165-gene `stage2/evaluation` release is not a current dependency. Tahoe has no
significance-label build. Training readers never load truth labels.
