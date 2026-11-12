# disease_data — public myocardial disease transcriptomes

The disease-side inputs for the Stage 2 gene-axis decision preparation
(ADR 0065) and the Stage 1 MIRI screen, split into one directory per **GEO
series** (项) under the external-data provenance contract
([ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md)).
Eight items are present. Each is:

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
external item; the `disease_data/` category has eight items and
78,661,383 B (75.0 MiB) of `source/` bytes.

## The eight 项

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

## Build recipes and pinned R

`script/fetch.sh` documents and re-runs the GEO FTP download into `source/`
(never edited afterwards); `script/build_deg.R` writes the DEG table and its
tracked JSON; `script/build_expression.R` writes the expression pair where the
series supports one. Each item carries its own copies of the shared
`script/{map_to_human.R,map_orthologs.R}` helpers (no cross-item dependency) and
adds a series-specific `script/lib_<GSE>.R` when the label map needs pinning
(GSE141910, GSE116250).

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

- **Stage 1 MIRI screen** — `stage1/*/08_miri_phenotype/miri_phenotype.py` reads
  `disease_data/<GSE>/local/deg/*.json` (and `source/*_deg_results.csv`) for the
  Stage-1 L1000 reference lists. That interface is unchanged.
- **Stage 2 gene axis** — `stage2/01_bases/execution/build_gene_axis.py`
  (ADR 0065, #183) reads the five human items' `local/deg/<GSE>_deg.csv` and
  `local/expr/` pair, plus the two mouse items' full-transcriptome
  `source/*_deg_results.csv` (DEG-only; no expression matrix is shipped). It
  never reads the L1000 lists. The evidence product and its six reported
  quantities are described in
  [`docs/research/stage2-stage3-gene-axis-rank-aggregation.md`](../docs/research/stage2-stage3-gene-axis-rank-aggregation.md).

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
