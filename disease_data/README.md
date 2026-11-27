# disease_data — public myocardial disease transcriptomes

The disease-side inputs for the Stage 2 gene-axis decision preparation
(ADR 0065) and the Stage 1 MIRI screen, split into one directory per **GEO
series** (项) under the external-data provenance contract
([ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md)).
The routing index and offline verifier enumerate the current items. The initial eight-series inventory below is historical; later source sections record additions. Each item follows:

```
<GSE>/
├── script/          per-series fetch/build/verify recipes + shared map_to_human.R / map_orthologs.R  [tracked]
├── source/          upstream GEO-derived acquisition files (bytes unchanged)                          [gitignored]
├── local/deg/       full-transcriptome <GSE>_deg.csv (ignored, regenerable) + <GSE>_deg.json (tracked)
├── local/expr/      <GSE>_log_expr.tsv.gz + <GSE>_samples.csv, where an expression matrix exists     [ignored]
├── author_log.json  GEO accession, derivation chain, local additions, verification                   [tracked]
├── manifest.sha256  per-file hashes over source/ only                                                [tracked]
└── .gitignore       declares the regenerable local additions (where needed)                          [tracked]
```

`script/verify.sh` checks each `manifest.sha256` offline, and a missing manifest
is a failure. Run `script/verify_all.sh` from the repository root to verify every
external item. Counts/bytes in the initial inventory below describe that historical acquisition, not the current complete-source release.

## Initial eight-series acquisition (historical)

| item | organism / platform | role and comparison | n | `source/` bytes |
|---|---|---|---|---|
| `GSE141910/` | human, Illumina HiSeq 2500 (GPL16791), LV bulk | **HF primary anchor**: HF vs non-failing (NF), limma | 200 HF (166 DCM, 28 HCM, 6 PPCM) / 166 NF | 60,049,299 (`GSE141910_RAW.tar` 60,026,880 + `filelist.txt` 22,419) |
| `GSE116250/` | human, Illumina HiSeq (GPL16791), LV bulk | **HF cheapest anchor + human post-MI (ICM) proxy**: HF vs NF, limma | 50 HF (37 DCM, 13 ICM) / 14 NF | 7,487,245 (`GSE116250_rpkm.txt.gz`) |
| `GSE243674/` | human, hiPSC-CM monolayer, NextSeq 550 (GPL21697) | **chemo-cardiotoxicity, independent**: 10 drug×time vs vehicle contrasts, DESeq2 | 72 wells (6 donors × 6 arms × 2 time points; 12 vehicle) | 2,424,618 (`GSE243674_Counts_RNA_ERMatthews.txt.gz`) |
| `GSE128188/` | human, paired LA+RA appendage, NextSeq 500 (GPL18573) | **atrial fibrillation (>5% prevalence category)**: AF vs SR, DESeq2 | 5 AF / 5 SR donors (20 libraries, donor-pooled) | 647,785 (`GSE128188_CountsV2_EdgeR.txt.gz`) |
| `GSE282618/` | human, post-mortem LV+RV, NextSeq 2000 (GPL30173) | **hypertensive heart disease, exploratory**: HHD vs control, DESeq2 | 3 HHD / 2 control subjects (10 libraries, subject-pooled) | 699,058 (`GSE282618_counts.txt.gz`) |
| `GSE160516/` | mouse (*Mus musculus*), Affymetrix Clariom S | MIRI: pooled I/R (6/24/72 h) vs Sham, limma; Stage 1 primary rank reference | 16: 4 sham + 12 I/R | 2,004,508 (`GSE160516_deg_results.csv`) |
| `GSE108940/` | mouse (*Mus musculus*), Agilent Whole Mouse Genome 4×44K GPL7202 | MIRI: WT I/R vs WT Sham, limma; provenance corrected in #176 | 24 total; 12 WT (6 I/R + 6 sham) | 2,766,695 (`GSE108940_WT_deg_results.csv`) |
| `GSE137293/` | pig (*Sus scrofa*), Illumina NovaSeq | scar vs border zone — an injury-severity gradient, **not disease vs healthy**; supplementary only | 8 pigs × 2 regions | 2,582,175 (`GSE137293_COUNTS.txt.gz` + `GSE137293_deg_results.csv`) |

The three Stage-1 MIRI items (`GSE160516`, `GSE108940`, `GSE137293`) were built
before this round and expose only DEG tables (GSE137293 additionally ships its
counts matrix); the five human series added for the gene-axis decision
preparation ship both a DEG table and a normalized log expression matrix.

### Retained, but **not** axis sources after the ADR 0066 rebuild

ADR 0066 requires a `{disease}-DEG` axis source to be **cardiomyocyte-derived**
(in-vitro CM/organoid/engineered-tissue bulk, or a single-cell CM pseudobulk).
Whole-heart / whole-left-ventricle **tissue bulk** is ineligible even though it
contains cardiomyocytes: its signal is dominated by cell-composition and fibrosis
shifts that do not describe the cardiomyocyte. The six whole-tissue series below
are therefore **retained** (items and their `source/`/`local/` provenance stay
byte-unchanged, and they remain available to Stage 1 / other analyses) but they
are **not** gene-axis sources and must not be re-admitted as one without a new
ADR:

| item | why not an axis source |
|---|---|
| `GSE141910/` | human LV whole-tissue bulk (HF vs NF) — cell-composition/fibrosis signal dominates |
| `GSE116250/` | human LV whole-tissue bulk (HF/ICM vs NF) — same reason |
| `GSE128188/` | human atrial appendage tissue bulk (AF vs SR) — not CM-derived |
| `GSE282618/` | human post-mortem LV/RV tissue bulk (HHD vs control) — not CM-derived |
| `GSE160516/` | mouse whole-LV tissue bulk MIRI — species and tissue-composition signal |
| `GSE108940/` | mouse whole-LV tissue bulk MIRI — same reason |

Their DEG tables still exist for their own purposes; the retired
`mouse_human_orthologs.csv` / `build_mouse_orthologs.R` path went with the old
gene-axis evidence (ADR 0066 *Removals*), so no Stage 2 step consumes a mouse
ranking any more.

## Per-item conventions for the new series

The five human series follow one interface, and it is the only thing the Stage 2
gene-axis step reads from an item (ADR 0065):

- **`local/deg/<GSE>_deg.csv`** — full-transcriptome, gene-level differential
  expression with a signed effect (`log2FC`), test statistic, p-value, BH `padj`,
  expression level and case/control `n`. Contrast direction is disease versus
  control (positive `log2FC` = up in disease); a multi-contrast source keeps one
  long table with a `contrast` column (GSE116250 has three; GSE243674 has ten).
  It is **never L1000-filtered**, is regenerable and gitignored (global `*.csv`
  rule), and carries effect sizes, not just a gene list.
- **`local/deg/<GSE>_deg.json`** — the small tracked summary: accession, species,
  platform, tissue/model, contrast, method, thresholds, `n_genes_tested`,
  `n_up`, `n_down` and digests of the tracked inputs.
- **`local/expr/<GSE>_log_expr.tsv.gz` + `local/expr/<GSE>_samples.csv`** —
  normalized, gene-level, log-scale expression (genes × samples) with a
  sample → case/control sheet. Each item's normalization is documented in its
  `author_log.json` (voom log2 abundance, log2(RPKM+1), log2(CPM+1),
  log2(CP10K+1)). The item's HVG rank is computed from this pair and its DEG rank
  from the DEG table; the axis step reads nothing else.
- **Item-level `.gitignore`** declares the regenerable additions so a fresh
  checkout stays clean: `GSE141910`, `GSE116250`, `GSE128188` and `GSE282618`
  carry a top-level file; `GSE243674` pins its expression outputs with a nested
  `local/expr/.gitignore` (`*` + `!.gitignore`). Every ignored path is declared in
  `local_additions[]`.

## Stage 2 gene-axis CM sources (ticket #187, ADR 0066)

Three items added for the Stage 2 gene axis are **cardiomyocyte-derived**, not
whole tissue: two GEO single-cell deposits and one cellxgene object. Each follows
the same five-member ADR 0030 contract, adds a CM-gate + per-donor pseudobulk
derivation under `local/`, and publishes one compact DEG table per dataset into
the historical [disease DEG release](https://github.com/Moke110/carps/tree/59b188ec0fe6386f2c916e062eaad310c20f00f8/stage2/03_data_train/03.1_hvg-deg-5k/historical/bases/output/disease-deg):

| item | deposit | CM gate | donors | feeds |
|---|---|---|---|---|
| `GSE270788/` | human LV 10x Multiome GEX SubSeries (23 per-GSM filtered h5 + cell metadata; `RAW.tar` never fetched) | `cell.type.l1 == "Cardiomyocyte"` (13,803 CM, 52–1,576/donor, all qualify) | 23 | `nHF` NICM; `ICM` ICM + AMI |
| `GSE255612/` | human left-atrium 10x snRNA (mtx + genes axis + cell metadata) | `cell_type == "Cardiomyocytes"` (32,438 CM) | 34 (18 AF / 16 control) | `AF` |
| `Reichart-CM/` | cellxgene CM-only h5ad (collection `e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb`, version `27fdbdf9-deae-4aed-bf7f-21964ee62137`; 311,418 cells) | CM-only object (`cell_type == "cardiac muscle cell"`) | 79 (52 DCM / 18 normal / 8 ARVC / 1 NCCM) | `nHF` genotype-negative DCM; ticket #188 export for `genetic_cardiomyopathy` |

Conventions specific to these three items:

- **`local/pseudobulk/`** — the per-donor pseudobulk, summed over each donor's CM
  cells **before** any test. GSE270788 and Reichart-CM store integer counts;
  GSE255612's deposit is log1p(CP10K) normalised expression (verified on every
  cell: `sum(expm1(value)) == 1e4`), so its pseudobulk is the per-donor mean
  CP10K. The matrices are declared additions and stay out of `source/`.
- **`local/deg/<item>_deg.csv` + tracked `_deg.json`** — the item-local
  full DEG table (long, one row per gene × contrast) and its tracked summary.
  The DEG is the single live implementation in
  `stage2_CarPS-VAE/util/output/api/stage2_util/derivation/disease_deg.R` (limma: voom for counts,
  `lmFit` for the log-scale deposit); the item scripts only bind inputs and run
  it. GSE270788's three arms share one control set, so its table comes from
  **one joint `~ 0 + arm` fit** over all 23 donors (three contrasts), not three
  case-control fits.
- **Tracked compact export** — `Reichart-CM/local/pseudobulk/…csv.gz` (axis genes
  × 78 donors, integer counts) and `Reichart-CM/local/genotype/…csv` (per-donor
  genotype split) are force-tracked because they are the ticket #188
  `genetic_cardiomyopathy` input. The single NCCM donor is excluded from that
  export; genotype-positive DCM and the 8 ARVC donors are exported, not used for
  `nHF`.
- **Erratum carried from #187.** GSE270788's deposited arm sizes are Donor 7 /
  NICM 6 / ICM 5 / AMI 5 (both the cell metadata grouped by `sample` and the
  series-matrix `!Sample_characteristics_ch1`), not the survey's `6 / 7`; the
  deposit is used unchanged and ADR 0066's roster counts are superseded.

`script/verify_all.sh` includes all three items; `bash script/verify.sh` in each
item checks `source/` offline against its `manifest.sha256`.

## Stage 2 isogenic CM models (ticket #188, ADR 0066)

Six **isogenic hiPSC-CM bulk** models join the `genetic_cardiomyopathy`
category, where the two Reichart genotype-split exports (genotype-positive DCM
and ARVC donors) are the tier-1 human CM-pseudobulk datasets and the six
in-vitro models enter at tier 2 under ADR 0066's **Isogenic exception**. Each
item follows the same five-member ADR 0030 contract, and one shared adapter
(`stage2_CarPS-VAE/util/output/api/stage2_util/derivation/bulk_deg.py`) binds the deposited matrix +
GEO series matrix and runs the one live DEG implementation
(`stage2_util/derivation/disease_deg.R`); the item script only declares its design.

| item | deposit | design | published dataset |
|---|---|---|---|
| `GSE224129/` | MYBPC3 (cMyBP-C) HCM engineered cardiac tissue, isogenic +/+ (4), +/- (4), -/- (4); RSEM gene matrix | one joint `~0 + group`; `+/-` and `-/-` vs `+/+`, collapsed by maximum effect | `gse224129__mybpc3_hcm` |
| `GSE234141/` | TTN truncating-variant DCM hiPSC-CM vs isogenic WT/WT, 2D monolayer and EHT, atrial- and ventricular-like, 3 replicates per cell (24 libraries) | one joint fit over 8 composite strata; four paired TTNtv-vs-WT contrasts, collapsed by maximum effect | `gse234141__ttntv_dcm` |
| `GSE304575/` | LMNA p.His222Pro DCM hiPSC-CM vs CRISPR-corrected isogenic control (6/6) | `~0 + group`; one contrast | `gse304575__lmna_h222p_dcm` |
| `GSE305648/` | RBM20 DCM hiPSC cardioids, WTC11: KO / PL / RQ / wildtype (3 each) | `~0 + group + batch` (batch rotated A/B/C); three mutant-vs-WT contrasts, collapsed by maximum effect | `gse305648__rbm20_dcm` |
| `GSE253223/` | ARVC-model wild-type hiPSC-CM, siPKP2 (4/day) vs siNeg (1/day), days 2/4/6/8 | `~0 + group + day`; four harvest days pooled with day blocked; one knockdown contrast | `gse253223__pkp2_knockdown_arvc` |
| `GSE194297/` | DMD point-mutation cardiac organoids day 56 vs CRISPR-corrected isogenic line (3/3) | `~0 + group`; one contrast | `gse194297__dmd_cardiomyopathy` |

Conventions specific to these six items:

- **`script/fetch.sh`** stores exactly two payloads per item: the deposited
  expression matrix and the GEO series matrix (the sample titles/characteristics
  that carry genotype, cell type, format, batch and harvest day). No other
  payload is fetched; both are byte-size pinned and verified offline by
  `script/verify.sh` against `manifest.sha256` over `source/`.
- **Sample binding is by name, never by position.** The matrix column tokens are
  matched to the series-matrix titles/descriptions with a recorded basis and a
  fail-closed bijection check, so no library can be silently dropped or twice
  matched.
- **`local/deg/<item>_expr.tsv.gz` + `_samples.csv` + `_deg.csv`** are the
  item-local regenerable derivation (matrix, design sheet, full long DEG table);
  the **tracked** `local/deg/<item>_deg.json` summary records the deposit
  digests, the sample-mapping basis, the isogenic-control arm, the
  strata/pooling decision, the design and its contrasts with `n_case`/
  `n_control`, the integerisation statistics when the deposit is not integer
  (rule, maximum absolute deviation, library-size scale, and how many genes'
  rounded rows are all-zero so the count model drops them), and the R/limma
  versions.
- **Integerisation is one shared, tested helper** (`bases.bulk_deg.integerise_counts`,
  round to nearest with ties to even). It is used for `GSE224129` (RSEM
  expected counts) and `GSE253223` (estimated counts); the other four deposits
  are integer and go straight to the strict-count voom path. The shared
  `disease_deg.R` count contract is **not** relaxed: #186/#187/#189 recorded
  products depend on it. A deposit that is already normalised takes the
  `lmFit` non-count path instead of being rounded.
- **`Reichart-CM/`** keeps #187's tracked pseudobulk + genotype export
  untouched; ticket #188 adds `script/build_genetic_deg.py`, which runs the same
  live implementation over those **tracked** inputs (read-only, the pseudobulk is
  not re-derived) and publishes the `dcm_genotype_positive-normal` (44 vs 18) and
  `arvc-normal` (8 vs 18) contrasts into `local/deg/reichart-cm_genetic_deg.csv`
  + a tracked `.json` summary.
- **Registration**: all six items are in `script/verify_all.sh`'s
  `expected_items` (additive), so a missing item or manifest turns the health
  check red.

## Build recipes and pinned R

`script/fetch.sh` documents and re-runs the GEO FTP download into `source/`
(never edited afterwards); `script/build_deg.R` writes the DEG table and its
tracked JSON; `script/build_expression.R` writes the expression pair where the
series supports one. Each item carries its own copies of the shared
`script/{map_to_human.R,map_orthologs.R}` helpers (no shared-helper dependency
between items) and adds a series-specific `script/lib_<GSE>.R` when the label map
needs pinning (GSE141910, GSE116250). One deliberate cross-category read exists:
`GSE243674/script/build_deg.R` maps its ENTREZID keys with the pinned
`pert_data/Annotation/` HGNC table (read-only, digest recorded in the author log),
because that mapping table is a repository asset rather than an item download.

Pinned R interpreters — do not substitute the default `Rscript` on `PATH`:

```bash
# Human DEG + expression builds (R 4.5.3; limma, DESeq2, edgeR, org.Hs.eg.db, jsonlite)
cd disease_data/<GSE>
/home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_deg.R
/home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_expression.R

# GSE108940 mouse->human ortholog remap (R 4.4.1; homologene, org.Mm.eg.db, org.Hs.eg.db).
# The `arc` env's R/Rscript wrappers are mis-patched; invoke exec/R directly.
cd disease_data/GSE108940
R_HOME=/home/chang/miniconda3/envs/arc/lib/R \
  /home/chang/miniconda3/envs/arc/lib/R/bin/exec/R \
  --no-save --no-restore --file=script/build_deg.R
```

Re-running a build reproduces its DEG CSV digest exactly (the tracked JSON
records the counts and input digests). `source/` bytes are never modified, so
`manifest.sha256` keeps verifying after a reprocess.

## GSE108940 is mouse, not rat

The item was originally labelled *Rattus norvegicus* / Agilent 8×60K and mapped
rat→human (`inTax=10116`). NCBI GEO says **GSE108940 is *Mus musculus*, GPL7202
(Agilent-014868 Whole Mouse Genome 4×44K G4122F), 24 samples**, and the local DEG
table's 28,037 `Gene` strings are the platform's mouse `GeneName` set (4,463
`*Rik`, 908 `Olfr*`, 3,561 `AK*`, 606 `ENSMUST*`; zero rat `Rn#`/`ENSRNO`/`RGD`
identifiers). The item is corrected to mouse with `inTax=10090`; the only L1000
list content change is **TP53 (newly UP)**. `source/` bytes and the manifest are
unchanged, and `stage1/` is not edited — the Stage 1 union impact (+TP53) is
recorded in
[`docs/research/stage2-stage3-gene-axis-gse108940-provenance.md`](../docs/research/stage2-stage3-gene-axis-gse108940-provenance.md).

## Stage-1 L1000 lists stay Stage-1-only

The three original MIRI items keep their `local/deg/*.json` human-symbol
**L1000-978 landmark** lists as Stage 1 inputs; that role and format are
unchanged. `GSE160516` and `GSE137293` are byte-unchanged; `GSE108940`'s list was
regenerated under the corrected mouse mapping (5 UP + 6 DOWN), the one content
change from #176. None of the three is a gene-axis source: an L1000-978 list
cannot rank 18,164 genes. The five new DEG tables and the two mouse
full-transcriptome `source/*_deg_results.csv` tables are used instead, and no new
artefact is filtered to the 978 landmarks. `GSE137293` (pig scar vs border) is
supplementary only — a severity gradient, never a disease-vs-healthy fusion
source.

## The `>5%` prevalence criterion

The prevalence taxonomy behind the new series is
[`docs/research/stage2-stage3-gene-axis-datasets.md`](../docs/research/stage2-stage3-gene-axis-datasets.md)
(GBD 2019 arithmetic over "diseases of the heart": IHD ≈56%, AF/AFL ≈17%,
RHD ≈11% globally, non-rheumatic valve disease ≥6.8%, HHD ≈5.3%). `GSE128188`
covers the AF category and `GSE282618` the borderline HHD category (exploratory,
n=5 subjects). No RHD or valve myocardium series was usable, and no human
acute-MI bulk case–control series exists.

## Downstream consumers

- **Local Stage 1 MIRI screen** — `stage1/08_miri_phenotype/miri_phenotype.py` reads
  `disease_data/<GSE>/local/deg/*.json` (and `source/*_deg_results.csv`) for the
  Stage-1 L1000 reference lists. That interface is unchanged.
- **Stage 2 gene axis** — the [HVG-DEG step](../stage2_CarPS-VAE/02_gene_axis/02.1_hvg-deg/README.md)
  owns published selection and evidence. Shared DEG workers and rank helpers
  are documented in [Stage 2 util](../stage2_CarPS-VAE/util/README.md).

## How a new GEO series would be added

1. Create `disease_data/<GSE>/` with `script/`, `source/`, `local/deg/` and, if
   an expression matrix is derived, `local/expr/`.
2. Put the upstream acquisition files in `source/` via `script/fetch.sh`; do not
   edit downloaded bytes.
3. Copy both shared helpers into `script/` and add `script/build_deg.R` (full
   transcriptome, not L1000) plus `script/build_expression.R` when applicable.
4. Write `author_log.json` (accession, organism, platform, derivation chain,
   `local_additions[]`), `script/verify.sh` and the item `.gitignore`; declare
   every generated path.
5. Add one line to `AGENTS.md`, one row to the table above, and the item to the
   `script/verify_all.sh` inventory; hash `source/` with
   `python3 .scratch/data-provenance-restructure/ops/hash_items.py --only disease_data/<GSE>`.
6. Run `script/verify_all.sh`, `script/check_paths.sh` and
   `script/check_tracked.sh` from the repository root.

## Unified TOX source admission (ADR 0072)

The independent pre-DEG ledger and its evidence are retained in the
[historical Bases records](https://github.com/Moke110/carps/blob/59b188ec0fe6386f2c916e062eaad310c20f00f8/stage2/03_data_train/03.1_hvg-deg-5k/historical/bases/README.md). These additions have immutable GEO matrices and sample
metadata, fetch recipes, offline manifests and per-item design records:

| items | deposited measurements | eligible comparison |
|---|---|---|
| GSE242692, GSE313746 | raw counts, symbols / Entrez respectively | DOX versus matched donor/time vehicle; 6R in GSE313746 is donor 6 technical replication |
| GSE255505 | fractional estimated raw counts | WT PSC-CM DOX versus DMSO, 3/3 |
| GSE232331 | raw counts, engineered heart tissue | DOX versus DMSO with experiment dates retained; ACLA/AMR identity/evidence unresolved |
| GSE235470 | Salmon gene counts, cardiac spheroids | unformulated DOX versus DMSO, 3/3; protein formulation SPEDOX is distinct |
| GSE157282, GSE106297 | FPKM | DOX versus control; log2(FPKM+1) non-count analysis, with co-exposures and replicate families explicit |

All seven are tier 2. Admission does not assert a DEG yield or a toxicity phenotype
at every dose/time. Per-source stratification and existing statistical tests remain
analysis requirements. Under ADR 0073 the source ladder processes every qualifying declared tier without a DEG-count stop; the former stopping rule is historical. `GSE243674/local/tox_evidence/` contains
independent regulatory evidence snapshots, reproducible with its
`script/fetch_tox_evidence.py`; these are declared local additions, separate from
its unchanged GEO source bytes.

## Complete-source reassessment (ADR 0073, #221)

The current registry collects all qualifying tiers without a DEG-count stop.
The earlier ADR 0066 whole-tissue exclusions and ADR 0071 stopping rationale
above describe the historical builds; current human atrial/ventricular bulk
qualifies at tier 3, including AF GSE128188 and the GSE116250 ICM arm.
The dated source-review records in the Stage 2 source registry retain finite
primary-source admission, rejection, uncertainty and reentry evidence.

| added item | input and analysis | disease/tier |
|---|---|---|
| GSE183852 | Full author Cardiomyocytes DCM/Donor DESeq2 table; missing padj retained, no local refit | nHF / 1 |
| GSE248187 | Full author DSP-knockout EHT/isogenic edgeR test table | genetic / 2 |
| GSE221945 | Full author TBX5 G145R Brugada/corrected DESeq2 table, effect sign reversed | genetic / 2 |
| GSE329393 | Deposited RBM20 counts, matched-family joint limma-voom fit, exact approved HGNC Entrez mapping | genetic / 2 |

Each item has immutable source bytes, offline verification, declared local
derivations and digest-pinned fetch/build recipes. Author-table adapters retain
full tested rows and author statistics; they do not assert a local raw-count
refit or infer independent donors from library counts.

GSE298023 is also reconstructed from the official 2026 Seurat object (4,289,920,683 bytes), using author `ident == Cardiomyocytes` and `heart_number` identities. Sparse integer counts are summed per donor across repeated scar/nonscar regions before one joint run-pattern-adjusted DCM/ICM-vs-NF fit. Actual CM-positive cohorts are 6 DCM, 7 ICM and 4 NF; controls are shared across the two disease records. The raw 142,024-cell matrix stays sparse; normalised/scaled assays are released before aggregation. The unchanged 10 GiB cap is required; rebuilding uses `script/build_pseudobulk.R` then `script/build_deg.py`.
