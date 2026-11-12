# Stage 2 gene-axis — GSE108940 provenance resolution

**Status:** resolution note (2026-11-13). Implements issue **#182** (which implements **#176**).
Not an ADR; it records a provenance correction and its Stage 1 impact.

**Scope.** The provenance of `disease_data/GSE108940/` is corrected in place: species,
platform, ortholog taxid and the derived L1000 lists. `source/` bytes are unchanged (the
manifest still verifies), and `stage1/` is not edited — the Stage 1 impact is recorded below
for a future ticket.

## Question

`disease_data/GSE108940/source/GSE108940_WT_deg_results.csv` (28,037 DEG rows) was
labelled *Rattus norvegicus* / Agilent SurePrint G3 8×60K and mapped rat→human
(`inTax=10116`), while NCBI GEO describes GSE108940 as *Mus musculus* on GPL7202. Which
species/platform actually produced the local table, and is the source still usable?

## Bottom line

**(a) It is the mouse GPL7202 series — no rat accession is involved.** The local CSV is
exactly the non-control `GeneName` set of the GSE108940 Agilent FEATURES files, in mouse
annotation namespaces; the rat label came from a mis-recording, not from different raw
data. The item is corrected to *Mus musculus* / GPL7202 / `inTax=10090`. The remapped L1000
list is **5 UP + 6 DOWN** (was 4 + 6): the only gene that changes is **TP53 (newly UP)**.
`source/` is byte-identical, so `manifest.sha256` and the item's `retrieved_at` /
`verified_at` / `hashed_at` are unchanged.

## Evidence

### 1. Independent GEO evidence (E-utilities / GEO FTP, fetched 2026-11-13)

| # | source | observation | verdict |
|---|---|---|---|
| 1 | `esummary` of GSE108940 (`db=gds`, uid 200108940) | `accession=GSE108940`, `taxon=Mus musculus`, `gpl=7202`, `n_samples=24`, `gdstype=Expression profiling by array`, `suppfile=TXT`, title "Cardiac Muscle Ring Finger-1 (MuRF1) and Transcriptional Regulation of Ischemia-Reperfusion Injury In Vivo" | mouse, 24-sample array series |
| 2 | `esummary` of GPL7202 (`db=gds`, uid 100007202) | "Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Probe Name version)", `taxon=Mus musculus`, in situ oligonucleotide; 41,534 60-mer probes / >41,000 mouse genes | 4×44K mouse array, not an 8×60K SurePrint G3 array |
| 3 | `GSE108940_series_matrix.txt.gz` header | `!Series_platform_id "GPL7202"`, `!Series_platform_taxid "10090"`, `!Series_sample_taxid "10090"`, 24 GSMs (GSM2916264–GSM2916287); design = WT/MuRF1-Tg × sham/I-R, 3 biological replicates × above/below LAD = 6 samples per group | taxid 10090; the WT sham / WT I-R design named in `build_deg.R` |
| 4 | GEO FTP `suppl/filelist.txt` | `GSE108940_RAW.tar` 380,979,200 B; 24 members `GSM*_..._GE2-v5_95_Feb07_*.txt.gz`, 15,763,484–15,957,459 B each (380,956,088 B member total); 6 WT_SHAM + 6 WT_IR + 6 TG_SHAM + 6 TG_IR | Agilent 4×44K FEATURES extractions, 6 per group; the 12 WT files are exactly the WT contrast |
| 5 | GEO FTP `matrix/` | `GSE108940_series_matrix.txt.gz` 4,956,665 B (processed matrix only) | no expression matrix is local; only the DEG table |

### 2. Row-by-row local evidence

The local DEG table's `Gene` column equals the platform annotation of the GSE108940 FEATURES
files, not an annotation from any other array:

| # | observation | verdict |
|---|---|---|
| 6 | One downloaded GSE108940 FEATURES file (GSM2916264, 45,018 feature rows) has **28,052 unique non-empty `GeneName`** values; removing the 15 control probes (`GE_BrightCorner`, `DarkCorner`, `ETG*` ×13) leaves **28,037 — exactly the unique `Gene` set of the local DEG CSV**. All 28,037 CSV genes are present in that file's `GeneName` column (100%); a second FEATURES file (GSM2916276, WT I-R) has the identical 28,052-name set (symmetric difference 0). | the DEG CSV is derived from these GPL7202 FEATURES files |
| 7 | The FEATURES rows carry mouse annotation (`Mus musculus ...` descriptions, `Mus musculus` RIKEN clones); the same strings appear in the CSV. | mouse origin |
| 8 | The FEATURES file has 45,018 rows — 4×44K scale (GPL7202), not the ~63,000-feature 8×60K SurePrint G3 layout claimed before. | not an 8×60K rat array |

### 3. Annotation-namespace census of the 28,037 DEG genes

| bucket | count | note |
|---|---|---|
| `*Rik` (RIKEN cDNA, mouse nomenclature) | 4,463 | mouse-only suffix |
| `Olfr*` (olfactory receptor, mouse nomenclature) | 908 | mouse-only prefix (rat orthologs are `Olr*`) |
| `AK*` (GenBank, RIKEN mouse cDNA clones) | 3,561 | mouse clone accessions |
| `ENSMUST*` (mouse Ensembl transcripts) | 606 | mouse-only accession namespace |
| `[A-Z]{2}[0-9]{6}` GenBank clone accessions (includes the AK set) | 4,501 | 940 are non-AK mouse EST/cDNA accessions |
| `TC*` (TIGR consensus) | 505 | legacy accession style |
| `Gm[0-9]+` (predicted mouse genes) | 105 | mouse-only |
| `Klra[0-9]+` (Ly49, mouse) | 19 | mouse |
| `RefSeq` `NM_`/`XM_`/`NR_`/`XR_` | 9 | |
| accession/transcript-ID strings, total | 5,621 (20.05%) | |
| in `org.Mm.eg.db` SYMBOL | 14,520 (51.79%) | |
| in `org.Rn.eg.db` SYMBOL | 12,113 (43.20%) | shared orthologous names |
| mouse-only / rat-only / shared / neither | 2,477 / 70 / 12,043 / 13,447 | |
| `ENSRNO*`, `RGD*`, rat `Rn.<UniGene>` IDs | 0 | no rat-specific identifiers at all |
| `Rn[0-9]+` (rat gene symbols) | 0 | — |

The 70 "rat-only" strings are annotation-era names that also exist in rat (`Igh-6`,
`Ly6c`, `Hist1h2ac`, …), not rat-origin genes. The earlier candidates in the ticket fact
table were prefix artefacts: `^Gm` = 119 counts canonical symbols such as `Gmppb`,
`Gmfg`, `Gmps`; `^Rn` = 86 counts mouse symbols such as `Rnf14`, `Rnps1`, `Rnaseh1`,
`Rnu33` — not a single rat `Rn#` symbol.

### Conclusion

(a) **Mouse GPL7202.** The source is usable; there is no mis-recorded rat accession to
identify, and the item does not need to be dropped. The rat label was a provenance error in
`author_log.json` / `local/deg/*.json` / `build_deg.R` only.

## Correction applied (all inside `disease_data/GSE108940/`)

| file | change |
|---|---|
| `author_log.json` | species `Mus musculus (mouse)`; method = GSE108940 RAW-tar Agilent GPL7202 FEATURES + limma; local additions mouse→human; `notes` rewritten with the mouse evidence. `source/` entries, `bytes`, `tree_sha256`, `manifest`, `retrieved_at`, `verified_at`, `hashed_at` unchanged because the source bytes are unchanged. |
| `local/deg/GSE108940_MIRI.json` | species `Mus musculus`; platform `Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (GPL7202)`; description now mouse MIRI + homologene `inTax=10090, outTax=9606`; up/down lists regenerated (5 UP + 6 DOWN). |
| `local/deg/GSE108940_human_l1000.csv` | regenerated with the mouse mapper; same schema (`human_gene, direction, logFC, P.Value`) and same thresholds (`P<0.05`, `|logFC|>0.3`). Never tracked (global `*.csv` rule). |
| `script/build_deg.R` | header/flow rewritten to mouse GPL7202; `library(org.Mm.eg.db)`; `map_deg_to_l1000(species="mouse")`; JSON metadata corrected. The script now also runs when the raw FEATURES files are absent (it then maps the existing `source/GSE108940_WT_deg_results.csv`), so `local/deg/` is reproducible without re-downloading 363 MB of raw data. |
| `source/`, `manifest.sha256`, `script/verify.sh`, `script/map_to_human.R` | **unchanged**; `sha256sum -c manifest.sha256` still passes (`4d7505e41b952177d3e983cfc1e7a065e1315b6adaad062c976c1a2087e09098`). |

### Old vs new L1000 list (the disease input Stage 1 consumes)

| | old (rat mapper, `inTax=10116`) | new (mouse mapper, `inTax=10090`) |
|---|---|---|
| UP | ACLY, CCNA1, OXCT1, PPARD (4) | ACLY, CCNA1, OXCT1, PPARD, **TP53** (5) |
| DOWN | CXCL2, DDX10, DECR1, IER3, IL13RA1, PRKAG2 (6) | unchanged (6) |
| differing genes | — | **TP53 added to UP** (from mouse-only symbol `Trp53`; logFC 0.364, P 0.0245) |

All 10 old rows are reproduced exactly by the rat mapper with the same helper and thresholds
(verified: set-equal to the old `GSE108940_human_l1000.csv`), so the only change is the one
gene the rat mapper could not resolve.

## Stage 1 impact (recorded here; Stage 1 is not edited)

**Consumers.** `stage1/xpert-fda-miri-{v1,v1.1,v2,v3}/08_miri_phenotype/miri_phenotype.py`
read every `disease_data/<GSE>/local/deg/*.json` (including `GSE108940_MIRI.json`) and write
`08_miri_phenotype/output/MIRI_up-down.json` (majority-vote union) and `MIRI_rank.json`.
`09_miri_efficacy_scores/miri_efficacy_scores.py` then consumes both (WCS/RGES, Zhang/ssCMap,
pathway RGES), and steps 10–11 build on it. Those four step-08 scripts label GSE108940
`"rat MIRI"` — a stale label left untouched with the rest of Stage 1.

**Effect of the correction.** Re-running the step-08 union logic with the corrected list:

| | old union | new union | difference |
|---|---|---|---|
| UP | 78 | 79 | **+TP53** |
| DOWN | 36 | 36 | none |
| conflicts | IER3, IL13RA1 (GSE160516 UP vs GSE108940 DOWN → tie-removed) | same | none |

`MIRI_rank.json` is built from the primary mouse series GSE160516 alone and is **unaffected**.
So a Stage 1 re-run of `08 → 09 → 10 → 11` would add TP53 to the MIRI reference and change the
efficacy scores and the tracked step-09 figures; that re-run is a **new ticket**, not part of
#182. Nothing under `stage1/` was edited.

## What a full-transcriptome recompute would need

The item has **no expression matrix** locally — only the limma DEG table
(`source/GSE108940_WT_deg_results.csv`, 2,766,695 B). A full-transcriptome recompute
(required to replace the 978-gene L1000 lists with an 18,164-axis-compatible DEG table)
would need:

| file | size | why |
|---|---|---|
| `GSE108940_RAW.tar` (`suppl/`) | **380,979,200 B (363 MB)** | 24 Agilent FEATURES `.txt.gz` (15.76–15.96 MB each) — the inputs `script/build_deg.R` steps 1–3 parse |
| `filelist.txt` | 2.6 KB | member names → GSM/group map |
| `GSE108940_series_matrix.txt.gz` (`matrix/`) | 4,956,665 B | processed matrix; sufficient for a matrix-level reanalysis, not for the FEATURES/limma path the script implements |

The recompute would re-run `build_deg.R` steps 1–3 over the 24 FEATURES files (12 WT samples
in the contrast), then the same mouse→human mapping at `P<0.05`, `|logFC|>0.3`; no download
was performed for this ticket.

## Reproduction (exact commands and environment)

Pinned R (the `arc` env's `R`/`Rscript` wrappers are mis-patched; `exec/R` works), with the
worktree's gitignored landmark CSV supplied via `ARC_L1000_CSV` (a full checkout resolves it by
default). Run from `disease_data/GSE108940/`:

```bash
cd disease_data/GSE108940
R_HOME=/home/chang/miniconda3/envs/arc/lib/R \
ARC_L1000_CSV=/home/chang/arc/pert_model/XPert/local/metadata/l1000_gene_info_978.csv \
/usr/bin/time -v /home/chang/miniconda3/envs/arc/lib/R/bin/exec/R \
  --no-save --no-restore --file=script/build_deg.R
# → "No raw Agilent FEATURES files in source/ — mapping the existing DEG table only"
# → "GSE108940 done: 5 UP + 6 DOWN (human L1000)"
# wall 0:02.84, peak RSS 447.6 MiB (before the run: free -g → 47 GB total, 42 GB available)
```

`script/build_deg.R` writes `local/deg/GSE108940_human_l1000.csv` and
`local/deg/GSE108940_MIRI.json`; it does not touch `source/` in this mode. Namespace census and
GEO evidence were produced with `org.Mm.eg.db`/`org.Rn.eg.db` (`keys(..., "SYMBOL")` + Entrez
`mapIds`) and E-utilities (`esearch`/`esummary`, `db=gds`) as listed above.

Verification (offline):

```bash
bash disease_data/GSE108940/script/verify.sh          # sha256sum -c ../manifest.sha256 → OK
script/check_paths.sh                                 # green
```

## Residual notes and out-of-scope flags

1. **`source/` is not a pure download.** Per `author_log.json`, `source/GSE108940_WT_deg_results.csv`
   is the limma output written by the historical pipeline, not a byte set downloaded from GEO
   (GEO's series has only `GSE108940_RAW.tar` + `filelist.txt`). ADR 0030's "source/ is the
   untouched download" therefore does not literally hold for this item. #182 requires the
   source bytes to stay identical, so this is recorded, not changed.
2. **Shared helper copy.** `script/map_orthologs.R` (byte-identical across the three MIRI items)
   still documents GSE108940 as rat and accepts only `rat`/`pig`. Editing it here would break the
   copies-identical invariant and the other items are outside this ticket's write scope; the
   helper does not participate in the corrected `build_deg.R` path (which sources
   `map_to_human.R:mouse2human` directly).
3. **Category-level docs still say rat.** `disease_data/README.md` (table row + "Rat vs human"
   section) and `disease_data/AGENTS.md` still describe GSE108940 as rat; they were outside the
   write scope of this ticket and need a follow-up edit.
4. **Stage 1 label.** `stage1/*/08_miri_phenotype/miri_phenotype.py` still prints/labels
   `"GSE108940": "rat MIRI"`; Stage 1 is frozen for this ticket (impact recorded above).
5. **L1000 only.** The corrected list remains a 978-landmark set; it is not usable as the
   Stage 2 18,164-gene axis input until the full-transcriptome recompute described above is
   done (a separate ticket).

## Appendix — evidence artefact digests (2026-11-13)

```
0350953c02b83d81b16ecb9bca9b6ba1d6f8ac6b5758a46f996377dae947ad4a  GSM2916264_FEATURES.txt.gz   (WT sham, above LAD)
e0332c155dd5b51cd617efd1fe5c2dd82fcae56c381f03c7c76c8aca6cbaa3d4  GSM2916276_FEATURES.txt.gz   (WT I/R, above LAD)
e20a8bd5dc4cdf4d91b0cea6bfaa4960b683a34bd8292d3b013dfac28859d125  GSE108940_series_matrix.txt.gz
bd2f3199bea338b87d12b41f451ea899ac3e52bc5fb62d8e935514504e45c62d  GPL7202.annot.gz
c4b8b30b8463ef065eac0517efcc286ec096652abf91cfe366024be03527bdd1  filelist.txt
```
