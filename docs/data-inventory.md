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
| `disease_data/` | 7,449,064,852 | 7.4 GB | 23 項 | source/ ignored; local/deg/ mixed (JSON tracked, CSV ignored) |
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
| `disease_data/GSE255612/` | 1,032,251,016 | 1.03 GB | yes — GEO **GSE255612**: Scanpy mtx (1,003,414,933 B), feature axis (277,568 B) and cell metadata (8,846,174 B); `script/fetch.sh` (parallel ranges). `local/` pseudobulk + DEG derived | re-download `source/`; tracked summary/manifest in git. `01_bases` declares it again as an `AF` tier-1 source under ADR 0071 (it measures 0 v2 DEGs and is excluded from the fusion with that reason) |
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
are unchanged. Measured 2026-11-14 after that derivation, its `du -sb` is
3,207,685,803 B: the 3,201,296,579 B row above predates #188, and the 6,389,224 B
difference is exactly the added `local/deg/` payloads
(`reichart-cm_genetic_deg.csv` 4,336,159 B, `reichart-cm_genetic_matrix.tsv.gz`
2,036,919 B, `reichart-cm_genetic_samples.csv` 1,953 B) plus tracked
`script/build_genetic_deg.py` (8,719 B), `reichart-cm_genetic_deg.json` (3,623 B)
and `author_log.json` growth (1,851 B).

### disease_data — Stage 2 gene-axis sources (#202-#206, ADR 0071)

The eleven items the v2 roster added or re-roled: the tier-2/3 clinical sources
(`GSE35108`, `GSE165303`, `GSE116250`, `GSE141910`, `GSE128188`, `GSE282618`,
`GSE269882`), the AF single-cell/array sources (`GSE138252`, `GSE238242`), the
cardiotox RNA-seq arm (`GSE243674`) and the Kuppe-MI cellxgene object. Each is an
ADR 0030 item (`script/`, ignored `source/`, declared `local/`, tracked
`author_log.json` and `manifest.sha256`); the item-local DEG tables stay in
`local/` and the compact published products are the step's `output/disease-deg/`.
`du -sb` of the whole item directory (derived `local/` included), 2026-11-15.

| 项 (item) | `du -sb` bytes | ≈ | reconstructible? (and from where) | backup expectation |
|---|---:|---:|---|---|
| `disease_data/GSE116250/` | 25,013,473 | 25.0 MB | yes — GEO **GSE116250** (GPL16791 Illumina HiSeq 2500); `script/fetch.sh` | re-download `source/`; tracked summary/manifest in git; `local/` derived |
| `disease_data/GSE128188/` | 4,690,042 | 4.7 MB | yes — GEO **GSE128188**; `script/fetch.sh` | same |
| `disease_data/GSE138252/` | 11,380,739 | 11.4 MB | yes — GEO **GSE138252** (CM-nuclear FANS bulk); `script/fetch.sh` | same |
| `disease_data/GSE141910/` | 92,963,119 | 93.0 MB | yes — GEO **GSE141910**; `script/fetch.sh` | same |
| `disease_data/GSE165303/` | 44,350,465 | 44.4 MB | yes — GEO **GSE165303** (GPL20795 Illumina HiSeq X Ten); the 17.6 GB RAW.tar is deliberately not fetched | same |
| `disease_data/GSE238242/` | 48,138,198 | 48.1 MB | yes — GEO **GSE238242** (multiome snRNA arm only; the snATAC arm is not fetched) | same |
| `disease_data/GSE243674/` | 34,332,098 | 34.3 MB | yes — GEO **GSE243674** (cardiotox RNA-seq counts) | same |
| `disease_data/GSE269882/` | 17,876,737 | 17.9 MB | yes — GEO **GSE269882** (atrial bulk RNA-seq, tier-3 AF) | same |
| `disease_data/GSE282618/` | 5,732,307 | 5.7 MB | yes — GEO **GSE282618** (exploratory HHD bulk; never admitted) | same |
| `disease_data/GSE35108/` | 15,102,745 | 15.1 MB | yes — GEO **GSE35108** (GPL6244 array) | same |
| `disease_data/Kuppe-MI/` | 1,879,197,376 | 1.88 GB | yes — CELLxGENE object, collection `8191c283-0816-424b-9b61-c3e1d6258a77`, version `4780eb7d-4c7c-4630-8705-471f0de8e27a`; `script/fetch.sh` | re-download; the tracked per-donor×zone pseudobulk export is in git |

`disease_data/` totals **7,449,064,852 B over 23 items** on 2026-11-15 (the
`3 項 (7.4 MB)` row in *Totals* is the three legacy MIRI items only, which happen
to sum to 7.42 MB — the coincidence is why the old total looked right).

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

## Stage 2 — the two active steps (2026-11-17, ADR 0062 / 0074)

`stage2/01_bases/` owns the foundation and `stage2/02_models_tenv-5k/` owns
the accepted single-dose native-model comparison. The previous `02_models_tenv`
step is retired after verified successor acceptance; its complete tracked and
ignored tree is retained in the archive listed below. Historical ADRs and sealed
records retain their original paths and scientific identities.

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
| `figure/` | 100,506 | `gene_axis_rank_sum_histogram.png` (two-panel evidence figure), its complete plotted-values CSV and the CSV-only renderer (the diagnostic PR figure was deleted by #192) | regeneration reads only `figure/data/*.csv` with `output/` products unreadable |
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
| `output/disease-deg/` | 138,019,836 | tracked source DEG/rank evidence: five diseases, 41 dataset records (8 nHF, 5 ICM, 5 AF, 11 genetic, 12 TOX); shared studies/control arms explicit, not 41 independent physical studies | `execution/build_disease_deg.py` plus independently replayed `execution/build_tox.py`; complete qualifying tier inventory, no count stop |
| `output/disease-fusion/` | 1,446,730 | full tier/rank/support unions plus five DEG1000 tables; union sizes nHF 1785, ICM 1671, AF 25, genetic 4399, TOX 9589; AF excluded by 100-gene gate | integrated successor driver, independent replay and public readers |

Both products are replayed by their drivers' `verify` action (rebuild in
`execution/checkpoints/` and compare every digest). The compact manifests,
mapping evidence and exclusion table are force-tracked; the count payloads are
gitignored.

### Unified TOX products ADR 0072

The accepted `output/disease-deg/tox/` contains twelve measured DEG/rank pairs,
independent admission and matched-well audits, a catalog and manifest. Compact
products are tracked; source measurements remain immutable external item bytes.
The former mechanism-class products are retired from live output and remain in
git history. Reconstruct TOX with `execution/build_tox.py publish`, verify with
its fresh source replay, and promote with the coherent integrated successor.

### ADR 0073 gene axes, figures and rule document (2026-11-16)

All compact scientific CSVs/catalogs/manifests, figure inputs/PNGs and the human
DOCX are tracked. Axis sizes are 2k/3k/4k/5k, exact and nested; selection uses
DEG1000/f1000 and complete-union rank denominators. Measured fallback is zero
at all four sizes, with fallback still implemented. Current source/provenance
and deterministic replay evidence is in `log/full_source_acceptance.json`;
older ADR 0066–0072 sizes/digests remain sealed historical acceptance records.

| stable path under `stage2/01_bases/` | file bytes | contents / reconstruction |
|---|---:|---|
| `output/gene-axis/` | 1,006,428 | Complete candidate table, four axes, composition/coverage and manifest. `execution/build_gene_axis.py publish` / `--verify` from final disease fusion + frozen identity + tenv-HVG. |
| `gene_axis_fill_rule.docx` | 104,322 | 13-page portrait A4 reading aid, five disease-category subtables (41 records) plus separate tenv-HVG table. `execution/make_gene_axis_docx.py`; input/tool/output pins and every-page QA in `log/gene_axis_document_verify.json`. Old live output/gene-axis copy removed. |
| `figure/gene_axis_rank_sum_histogram.png` | 95,996 | Four-axis (2×2 panel) histogram, rendered only from matching tracked figure/data CSV. |
| `figure/data/gene_axis_rank_sum_histogram.csv` | 18,622 | All candidate/selected histogram bins and axis panel labels. |
| `figure/gene_axis_source_topn_coverage.png` | 75,469 | Four-axis coverage bars using own DEG1000/group denominators; CSV-only renderer. |
| `figure/data/gene_axis_source_topn_coverage.csv` | 1,518 | Explicit complete coverage plot values; excluded AF unavailable with reason. |

### Paired 5k-local expression (ADR 0077)

Only `output/view-5k/` is live. The unchanged 4,478 tenv and 66,718 tahoe_bulk
rows use selected G_5000 int64 totals for CP10K/natural log1p, individually
normalized matched vehicles, and delta = X-control. `stage2-view-v3` binds the
unchanged gene axis, corpus and denominator provenance. All four named axes
remain; 2k/3k/4k expression is retired. Rebuild/replay only this pair with
`execution/build_view.py publish|verify`; read through bounded `bases.view.open_view`.

| product | bytes |
|---|---:|
| `output/view-5k/tenv.h5ad` | 138,696,684 |
| `output/view-5k/tahoe_bulk.h5ad` | 2,942,840,025 |
| `output/view-5k/library_totals.csv` | 9,366,916 |

These three ignored/reconstructible files total **3,090,903,625 bytes**.
Source/corpus/HVG/four-axis and historical downstream bytes remain unchanged.
Old measurements stay in sealed receipts; recovery views outside formal output
are temporary reconstruction evidence, not live model inputs.

### `stage2/02_models_tenv-5k/` — accepted three-fold comparison (#226/#231)

Canonical ignored payloads preserve the accepted complete 54-run matrix and new
18-slot Mean release. Sizes below are regular-file byte sums, excluding directory
metadata. The complete candidate duplicates published runs to retain registrations,
plan, resources and replay evidence; it is required for no-training matrix replay.

| path under `stage2/02_models_tenv-5k/` | files | bytes | contents / reconstruction |
|---|---:|---:|---|
| `output/inputs/` | 12 | 98,265,302 | 1,402 treated wells, 459 parents, 192 unique controls, 5,000 genes, exact three-fold splits/Morgan/manifest; `execution/prepare.py --view ../01_bases/output/view-5k --candidate <empty>` |
| `output/means/` | 62 | 142,776,780 | 18 baseline slots, 15 applicable and three cold-drug Mean_drug N/A; `execution/means.py --inputs output/inputs --candidate <empty>` |
| `output/runs/` | 594 | 2,296,905,455 | 54 newly trained float32 checkpoints/predictions, identities, curves, exact replay and execution receipts; `execution/run_matrix.py --all` with explicit source/input/Mean paths and a new candidate |
| `output/comparison/` | 9 | 13,245,103 | full per-well/parent/seed metrics, paired bootstrap, selection and plot values; `execution/run_matrix.py --accept-only` |
| `execution/checkpoints/matrix-threefold-f32-lr1e4/` | 840 | 2,325,763,750 | full accepted candidate, including all registrations, plan and raw resources; preserved from the verified formal generation |

Tracked report products are `figure/{metrics,paired_differences,learning_curve}.png`,
their complete CSV values under `figure/data/`, the CSV-only renderer and compact
acceptance/resource/environment records under `log/`. Canonical payload SHA256
preservation is recorded in `log/231-canonical-preservation.json`; replay and rebuild
commands are in the step README. Immutable external model sources remain owned by
`pert_model/PRnet` and `pert_model/TranSiGen`. Stage 3 is unchanged.

### Retired `stage2/02_models_tenv/` — preserved old accepted result

The old tree at commit `cf27c49` is archived at
`/home/chang/arc-retired-assets/tenv-02-pre-231-cf27c49/old02.tar`:
214 regular files, 3,069,022,388 source bytes; archive SHA256
`27175ec1d4dca1654500278b79fa3ceb439054f015f9c3b0488fe40b7cf68d08`.
`manifest.json` beside it records every entry; `log/231-retirement.json` in the new
step records renewed source/archive verification and the old→new mapping. This is
same-machine retention, not an offsite backup. Historical source/figures/records
remain accessible through Git history; ignored fits remain recoverable from the
verified archive. Earlier failed five-fold successor candidates retain their
original worktrees and evidence and are not counted as this accepted matrix.

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

## ADR 0072 TOX admission inputs

Seven additional GEO items under `disease_data/` (GSE106297, GSE157282, GSE232331, GSE235470, GSE242692, GSE255505, GSE313746) retain ignored immutable expression and series-matrix bytes; per-item manifests and fetch/verify recipes pin them. GSE243674 declares ignored local regulatory evidence snapshots, rematerialized by `script/fetch_tox_evidence.py`. The GSE244740 author Table 1 PDF is pinned in its manifest. Pre-analysis admission candidates remain in `stage2/01_bases/execution/checkpoints/tox-admission/`, generated by `build_tox_admission.py`; they do not replace the accepted disease release.

### ADR 0073 complete-source reassessment (#221)

GSE183852, GSE248187 and GSE221945 retain ignored complete publisher/author
differential-test workbooks/tables; GSE329393 retains its ignored GEO count
matrix and metadata. Each item's `manifest.sha256` and `author_log.json` pin
the exact bytes and fetch recipe. Their ignored full DEG CSVs are regenerated
by the declared `script/build_deg.py`; tracked summaries pin output digests.
Preserve or re-materialize these bytes in the main checkout before deleting
the implementation worktree.

GSE298023 retains the ignored official 2026 Seurat RDS (4,289,920,683 bytes; SHA256 `78b3667f3715ed847fd2229dcb8000867bc2b90f79f70795551433ebd51e6cb4`). `script/fetch.sh` resumes acquisition and verifies it offline. Declared ignored per-heart CM pseudobulk/count and DEG files are rebuilt by its two source-owned generators; source bytes must survive worktree cleanup.
