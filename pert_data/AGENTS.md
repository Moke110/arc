# AGENTS.md — pert_data (cardiomyocyte drug-perturbation data lake)

Data home for the new-arc cardiomyocyte (CM) perturbation pipeline. Holds raw/
processed **transcriptomic** perturbation data + the LINCS L1000 pretrain corpus.
Papers for each dataset live in `reference/pert_data/` (PDF + extracted `.md`).

**Status 2026-08-18:** the GSE146095/GSE146096 raw GEO deposit was **deleted**
from `dtoxs_primary/` (its data is a strict subset of the LINCS LDS-1593
package); the LINCS LDS-1593 well-level package was **released flat into
`dtoxs_primary/`** (the `lincs-lds-1593/` wrapper removed — `Data/`,
`Metadata/`, `ProcessingPipeline/`, `ReadMe.txt`, `LDS-1593.tar.gz`,
`DToxS_summary.md` now live directly under `dtoxs_primary/`); and the
**v1-era package + v1's frozen finetune data were frozen byte-identically at
`dtoxs_primary/v1-frozen/`** (what xpert-fda-miri-v1 used). Earlier: 2026-08-17
`dtoxs-xpert/` deleted / re-homed, `cm/` grouping removed (ADR 0011);
GSE221595 removed 2026-08-12 (user decision). All accessions verified against
GEO on 2026-08-12. See `reference/pert_data/publications.md`.

## Dataset overview

| Dataset dir | GEO accession(s) | Cell type | Platform | Drugs/chemicals | Samples | Exposure | Citation |
|---|---|---|---|---|---|---|---|
| `dtoxs_primary/` | LINCS LDS-1593 well-level package (Level 3 of LDG-1452; formerly also GSE146096/GSE146095, deleted 2026-08-18 as a strict subset) | primary human CM-like (PromoCell, 4 donors) | 3′-DGE RNA-seq | 47 small molecules + TRS antibody + IGF-1/TNF-α proteins (+ non-KI perturbations) | 956 well records | 48 h, single dose | van Hasselt *Nat Commun* 2020;11:4809 |
| `dtoxs_ipsc/` | GSE217421 (Conv, 1171) + GSE217423 (DGE, 344), SuperSeries GSE217424 | iPSC-CM (6 donors) | bulk RNA-seq | 54 FDA drugs (23 KIs, 4 mAbs, 4 anthracyclines, 7 cardiac-acting, 16 non-cardiac-acting) | 1515 | 48 h, single dose | Hansen *Nat Commun* 2024;15:7968 |
| `gse244740/` | GSE244740 | hiPSC-CM | ScreenSeq (GPL24676) | 42 (33 cardiotoxicants + 9 non-cardiotoxicants) | 3072 | 24 h, 8 concentrations × 3 reps | Rosell-Hidalgo *Expert Opin Drug Metab Toxicol* 2024;20(7):685–707 |
| `gse262419/` | GSE262419 (SuperSeries; GSE244109 is a contained SubSeries — do not download separately) | iPSC-CM | TempO-Seq (GPL16791, ~20k targeted genes) | 464 chemicals, 12 classes | 4558 | 24 h transcriptomic; 90 min functional | Tsai *Chem Res Toxicol* 2024;37(8):1428–1444 |
| `gse247541/` | GSE247541 | hiPSC-CM (hCM) | bulk RNA-seq | 4 antivirals (apilimod, remdesivir, ritonavir, lopinavir) | 15 (3 DMSO + 4×3) | 6 d, multi-dose | Xu *Adv Sci* 2022;9(30):e2203388 |
| `l1000-xpert/` | GSE92742 / GSE70138 (datalake) | ~70 cancer cell lines (no CM) | L1000 | ~20k compounds | 3 h5ad (pretrain) | — | XPert pretrain-only, never for CM eval or new-model dev |

## Per-dataset contents

### `dtoxs_primary/` — DToxS primary CM (van Hasselt 2020): released LINCS LDS-1593 + v1-frozen
- **`Data/` `Metadata/` `ProcessingPipeline/` `ReadMe.txt` `LDS-1593.tar.gz`** —
  the **LINCS LDS-1593 well-level package** (Level 3 of LINCS deposit
  LDG-1452, DToxS center), released flat into `dtoxs_primary/` on 2026-08-18
  (the former `lincs-lds-1593/` wrapper was removed): per-batch well-level
  read-counts (`Data/Level3/Sequencing Experiment Number {1,2,4,5}/...`, 956
  well records, uncompressed), the LINCS metadata suite (`Metadata/` incl.
  `Small_Molecule_Metadata.txt` with SMILES + batch codes), the DCIC design
  xlsx, processing specs, and `LDS-1593.tar.gz` (the raw portal tarball).
  Byte-identical to the former `pert_data/dtoxs-xpert/` deposit (deleted
  2026-08-17). It carries the **non-KI perturbations** (estradiol, amiodarone,
  cyclosporin-A, milrinone, anthracyclines, ...), antibodies (TRS), proteins
  (IGF/TNF), and the CTRL-DMSO / CTRL-EtOH vehicle wells.
- **`v1-frozen/`** — the **v1-era data**, frozen byte-identically on 2026-08-18
  so `workflow/xpert-fda-miri-v1` stays reproducible after the GSE deletion +
  release: `dtoxs.h5ad` (a copy of v1's frozen `input/dtoxs.h5ad`, 521 samples
  / 36 drugs) + a byte-identical copy of the well-level package
  (`Data/` `Metadata/` `ProcessingPipeline/` `ReadMe.txt`). v1's step-02
  baseline reads from here (`pert_data/dtoxs_primary/v1-frozen/Data/Level3`).
- **Deleted 2026-08-18:** the **raw GEO deposit** (GSE146096 DGE + GSE146095
  Conv files) was removed — its content (23 KIs + CTRL, and the CER/IBR/LEN/NIN
  Conv wells) is a strict subset of the LINCS package. SuperSeries
  **GSE146097** was never downloaded (sub-series convention, same as
  `dtoxs_ipsc/`).

### `dtoxs_ipsc/` — DToxS iPSC-CM (Hansen 2024)
- `GSE217421_*` (Conv arm) and `GSE217423_*` (DGE arm): RNA-seq read counts, experiment design, configs, series matrices, source code, README. `GSE217421_Reference-Library.tar.gz` is ~1 GB (reference genome/lib — do not unpack casually).
- SuperSeries `GSE217424`. Transcriptomics only in GEO; WGS of the 6 cell lines is in dbGAP (phs002088.v1.p1), donor-level not drug-level.

### `gse244740/` — Rosell-Hidalgo (ScreenSeq)
- `GSE244740_processed_data_counts.txt.gz` (ScreenSeq counts) + series matrix. GPL24676.
- 42 compounds × 8 concentrations × 3 replicates, 24 h; DMSO + sunitinib controls.
- GEO holds only ScreenSeq. The paper's Ca²⁺ transient + high-content imaging + ATP data are in the paper/supplementary (see `reference/pert_data/gse244740.md`).

### `gse262419/` — Tsai (TempO-Seq)
- **16 per-plate CSVs** `GSE262419_Plate01..16.csv.gz` + `GSE262419_hash.csv.gz` + series matrix. (The old `GSE262419_iCM_normcounts.csv.gz` name from the design doc does **not** exist — the 16 plate CSVs replace it.)
- 464 chemicals × 8 concentrations; 16 plates × 29 test chemicals + vehicle/media/positive controls per plate.
- GEO holds expression only. Functional + viability data are in the paper's Supporting Information.

### `gse247541/` — Xu (COVID-19 antivirals)
- `GSE247541_raw_counts.csv.gz` + series matrix. Small (15 samples).
- 4 cardiotoxic antivirals; paper also screened 21 repurposed drugs and validated rescue compounds. Functional/imaging data in the paper only.

### `l1000-xpert/`
- `data/l1000_mdmt_full_336852.h5ad`, `l1000_mdmt_268022_pretrain.h5ad`, `l1000_mdmt_example.h5ad`; `L1000_summary.md`.
- Pretrain corpus for the XPert drug-response latent space / drug encoder; never enters CM evaluation, and **not used in new perturbation-model development** (kept only as the XPert pretrain unit).

## Removed datasets (provenance only)

- **GSE221595 (8 TKIs, hiPSC-CM)** — removed 2026-08-12 (user decision); directory deleted, planning docs updated. Its GEO series matrix was a single pooled pseudo-sample (per-drug data lived in `GSE221595_RAW.tar`). Do not reintroduce without user request.

## Per-drug non-transcriptomic metrics (papers + supplementary)

Beyond transcriptomics, per-drug readouts measured **on the CMs**:

| Dataset | Non-transcriptomic per-drug metrics | Coverage |
|---|---|---|
| DToxS primary | cell viability (Hoechst nuclear count), mitochondrial stress (MitoTracker Red) | dose-response, 48 h — but only a **subset of selected KIs**, not all 26; poorly correlated with clinical risk (paper Fig 1e/Supp, not in deposit) |
| DToxS iPSC-CM | **none** on CMs (transcriptomics only; WGS is donor-level) | — |
| GSE244740 | Ca²⁺ transient: amplitude, frequency, peak width, FWHM, rise/decay times, peak count (acute + 24 h); HCI: cell count, nuclear size, DNA structure, mitochondrial mass, MMP (TMRE), cellular Ca²⁺; cellular ATP (CellTiter-Glo) | all 42 compounds, 8-point dose-response |
| GSE262419 | Ca²⁺ flux functional: beat frequency, QT prolongation, asystole; cell viability; mitochondrial health (intact-mitochondria area, mitochondrial integrated intensity) | all 464 chemicals, concentration-response, 90 min exposure |
| GSE247541 | viability (calcein-AM/PI); apoptosis (TUNEL); MMP (JC-1); sarcomere organization (α-actinin); Ca²⁺ transients (Fluo-4, incl. 1–3 Hz pacing); hEHT contractile function | 4 cardiotoxic drugs (+ rescue screen) |

**No non-transcriptomic metric is measured across all 5 datasets.** Cell viability is
the closest to universal (4/5 — missing in DToxS iPSC-CM); Ca²⁺/beat-rate functional
readouts are in 3/5 (GSE244740, GSE262419, GSE247541), absent from both DToxS arms.

## Gotchas

- All CM datasets are log2/read-count expression space — z-score the predicted DEG in the ranking pipeline (workflow step 06); do not z-score at training.
- TempO-Seq (GSE262419) is probe-based; probe→gene mapping needed in step_04; platform imbalance risk (~70% of training samples).
- Concentration-response datasets (GSE244740, GSE262419) need per-plate matched controls for Δ computation.
- Only `Small_Molecule_Metadata.txt` (dtoxs_primary/Metadata) and resolved SMILES carry structures; antibodies/proteins have none. Chemicals without SMILES are excluded from training (user decision).
- `dtoxs_primary/DToxS_summary.md` is the well-level field-coverage summary (956 wells; 206 control; 750 treated; 672 wells SMILES-resolvable); `State` codes like `END`/`SUN+`/`DIC` have inconsistent metadata aliases — match to canonical IDs carefully.
- `dtoxs_primary/v1-frozen/` is a byte-identical frozen copy of the v1-era well-level package + v1's frozen `dtoxs.h5ad` — do not modify; the live package lives at the `dtoxs_primary/` root.
