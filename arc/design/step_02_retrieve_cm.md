# Step 2 — Data Retrieval: Cardiomyocyte Datasets

## Goal

Download all 5 cardiomyocyte drug-perturbation datasets from GEO and extract processed expression matrices plus experiment-design metadata. GSE221595 is removed from the pipeline (user decision 2026-08-12) — see step_01 §1.5b.

## Dependencies

- **Step 1 (Data Inventory):** Provides the validated GEO accessions, expected sample counts, and dataset metadata needed for retrieval.

## Inputs

- Validated GEO accessions from step_01:
  - DToxS primary CM: GSE146096, GSE146097
  - DToxS iPSC-CM: GSE217421, GSE217423, GSE217424
  - Rosell-Hidalgo 2023: GSE244740 (resolved)
  - GSE262419 (464 chemicals, TempO-Seq, ~4,558 samples)
  - GSE247541 (4 antiviral drugs, ~15 samples)
- GEO database (accessible via geoquery R package or direct FTP)
- DToxS data portal (dtoxs.org) and LINCS Data Portal as alternative sources for curated versions

## Outputs

- For each dataset: raw or processed expression matrix (read counts or normalized expression)
- For each dataset: sample-level metadata including drug name, dose, timepoint, donor, cell type, platform, batch
- For each dataset: gene-level identifiers (probe IDs, gene symbols, or Ensembl IDs)
- All outputs saved to shared checkpoint storage with dataset-specific subdirectories

## Behavior

### 2.1 Retrieval strategy

For each dataset, attempt retrieval in this order:
1. **geoquery (R):** Query GEO via `GEOquery::getGEO()` for the series accession. Download the series matrix file which contains processed expression data + sample metadata.
2. **Direct FTP:** If geoquery fails or the series matrix is too large, download supplementary files directly from `ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnnnn/GSEXXXXX/suppl/`.
3. **DToxS portal:** For DToxS datasets, check dtoxs.org for curated normalized fold-change matrices (the van Hasselt 2020 paper mentions these are deposited there).
4. **LINCS Data Portal:** For DToxS datasets, check the LINCS Data Portal as an alternative source.

### 2.2 Per-dataset retrieval

**DToxS primary CM (GSE146096, GSE146097):**
- GSE146096 is the DGE (differential gene expression) series; GSE146097 is the full transcriptome series.
- Retrieve processed read-count matrices + metadata.
- Expected: ~502 samples, ~26 kinase inhibitors, 4 donors, 3'-DGE RNA-seq platform, 48h treatment.
- Extract metadata fields: drug name, dose, donor ID, treatment time, batch.

**DToxS iPSC-CM (GSE217421, GSE217423, GSE217424):**
- Three related series covering the same study (Hansen 2024).
- Retrieve bulk RNA-seq processed data + metadata.
- Expected: ~1171 samples, 54 FDA drugs (23 KIs, 4 mAbs, 4 anthracyclines, 7 cardiac-acting, 16 non-cardiac-acting), 6 donors, 48h.
- Extract metadata fields: drug name, drug class, dose, donor ID, treatment time.

**GSE221595: REMOVED.** Not retrieved. Excluded from the pipeline per user decision (2026-08-12).

**Rosell-Hidalgo 2023 (GSE244740):**
- Accession **GSE244740** (resolved & verified 2026-08-12). 3072 samples, GPL24676, ScreenSeq.
- hiPSC-CM treated with 33 cardiotoxicants + 9 non-cardiotoxic compounds (42 total), 8 concentrations × 3 replicates, 24 h.
- Processed counts (`GSE244740_processed_data_counts.txt.gz`) + series matrix already downloaded to `pert_data/gse244740/`.
- Extract metadata fields: drug name, cardiotoxicity class (cardiotoxicant vs non-cardiotoxic), concentration, timepoint (24 h).

**GSE262419 (Tsai/Rusyn 2024, EPA-funded):**
- SuperSeries with ~4,558 samples — the largest single download in this project. Budget extra time.
- Platform: GPL16791 (Illumina HiSeq 2500), TempO-Seq (~20,000 targeted genes).
- Supplementary files: **16 per-plate CSVs** (`GSE262419_Plate01..16.csv.gz`) + `GSE262419_hash.csv.gz` (verified 2026-08-12; the design doc's earlier `GSE262419_iCM_normcounts.csv.gz` does not exist). All 17 downloaded to `pert_data/gse262419/`.
- BioProject: PRJNA1091686.
- Expected: 464 chemicals across 12 chemical classes, concentration-response design (multiple concentrations per chemical), iPSC-CM.
- Extract metadata fields: chemical name, chemical class (12 categories), concentration, timepoint, donor/batch/plate identifier, control_flag.
- **Do NOT separately download GSE244109** — it is a SubSeries of GSE262419 and its PFAS data is fully contained within the SuperSeries. Downloading both would create duplicate samples.
- Note: The TempO-Seq platform produces probe-level expression data (not transcript-level). Gene identifiers will be probe IDs that need mapping to gene symbols in step_04.

**GSE247541 (COVID-19 antiviral cardiotoxicity):**
- Small dataset: ~15 samples (3 DMSO controls + 4 drugs × 3 replicates).
- Platform: RNA-seq.
- Expected: 4 antiviral drugs (apilimod, remdesivir, ritonavir, lopinavir), iPSC-derived cardiomyocytes (hCM).
- Extract metadata fields: drug name, dose, timepoint, control_flag.
- Quick download; parse metadata from sample titles/characteristics.

### 2.3 Metadata extraction

For every retrieved dataset, extract and standardize these metadata fields:
- **drug_name:** raw drug/chemical name as recorded in the GEO sample title/characteristics
- **dose:** concentration value + unit (if available) — the treatment level for single-dose designs
- **concentration:** numeric concentration value for concentration-response designs (distinct from dose — concentration is the specific value in a concentration-response series; for single-dose datasets, concentration = dose). For GSE262419, this is a key variable with multiple values per chemical.
- **timepoint:** treatment duration
- **donor_id:** donor or cell line identifier
- **cell_type:** primary cardiomyocyte or iPSC-derived cardiomyocyte
- **platform:** sequencing platform / library prep method
- **study:** dataset/study identifier (DToxS-primary, DToxS-iPSC, Rosell-Hidalgo, GSE262419, GSE247541)
- **batch:** batch or plate identifier (if available)
- **control_flag:** whether the sample is a vehicle/DMSO control

### 2.4 Control sample identification

For each dataset, identify which samples are vehicle/DMSO controls. These are needed for Δ computation in step_04. Controls may be labeled as "DMSO", "vehicle", "untreated", or dose=0. Document the control identification logic per dataset.

For GSE262419: the concentration-response design uses DMSO/vehicle controls per concentration plate. Each concentration plate should have its own matched control. Verify that controls are present for each plate/batch and document the plate-control matching logic — this is critical for Δ computation in step_04, since each treated sample must be matched to the correct concurrent control.

For GSE247541: 3 DMSO controls are explicitly present. Match each drug replicate to the appropriate control (if plate-matched, use plate control; otherwise use the mean of the 3 DMSO controls).

## Unit tests

1. **Download completeness:** For each dataset, verify the number of downloaded samples matches the expected count from step_01 (within tolerance). Flag any dataset with >10% sample count discrepancy. For GSE262419, verify ~4,558 samples are present. For GSE247541, verify ~15 samples.
2. **Metadata field presence:** For each dataset, verify all 10 metadata fields (drug_name, dose, concentration, timepoint, donor_id, cell_type, platform, study, batch, control_flag) are extractable. Fields that are unavailable in the source data should be documented as "not available" rather than missing. For single-dose datasets, concentration may equal dose.
3. **Control identification:** For each dataset, verify at least one control sample is identified. If no controls are found, flag the dataset — Δ computation will fail without controls. For GSE262419, verify controls are present per plate/batch.
4. **Drug name extraction:** For each dataset, verify the number of unique drug/chemical names extracted matches the expected count from step_01. Flag discrepancies. For GSE262419, verify ~464 unique chemicals are extracted.
5. **Expression matrix integrity:** For each dataset, verify the expression matrix has no all-zero rows (genes) or all-zero columns (samples). Verify matrix dimensions match sample count × gene count. For GSE262419 (TempO-Seq), verify probe-level data is present and probe IDs are extractable for gene mapping in step_04.
6. **Sub-series dedup:** Verify GSE244109 was NOT separately downloaded. Verify no duplicate samples exist between GSE262419 and any GSE244109 data that may have been inadvertently retrieved.

## Integration test

Validate that all 5 cardiomyocyte datasets are retrieved and metadata is consistent:
- All 5 datasets have expression matrices with non-zero dimensions.
- All 5 datasets have metadata tables with the required fields.
- The union of unique drug/chemical names across all 5 datasets is reported and cross-referenced with step_01's expected drug/chemical counts.
- Control samples are identified in all 5 datasets.
- No dataset has a sample count discrepancy >10% from step_01 expectations.
- The total combined sample count across all 5 datasets is reported.
- GSE244109 was not separately downloaded (no duplicate samples from the GSE262419 SuperSeries).
- GSE262419's TempO-Seq probe IDs are documented for gene-symbol mapping in step_04.

## Acceptance criteria

- [ ] All 5 cardiomyocyte datasets downloaded from GEO
- [ ] GSE221595 absent (removed per user decision 2026-08-12)
- [ ] Expression matrices extracted for each dataset
- [ ] Metadata tables with 10 standardized fields for each dataset
- [ ] Control samples identified in each dataset (per-plate for GSE262419)
- [ ] Drug/chemical names extracted and unique count per dataset matches step_01
- [ ] Total combined sample count reported
- [ ] GSE244109 not separately downloaded (sub-series dedup verified)

## Risks & notes

- GEO supplementary files can be large and slow to download. Budget ~2–4h for the original datasets, plus extra time for GSE262419 (~4,558 samples, the largest single download). GSE247541 is small and quick.
- Some datasets may store expression data as supplementary files (not in the series matrix), requiring FTP download and manual parsing.
- DToxS datasets span multiple GSE accessions that must be combined — verify they don't have overlapping samples.
- Rosell-Hidalgo accession is resolved (GSE244740) — no author request needed.
- Metadata formats vary wildly across GEO series — drug names may be embedded in sample titles, characteristics columns, or supplementary metadata files. Expect per-dataset parsing logic.
- GSE262419's TempO-Seq data is probe-based, not transcript-based. Probe-to-gene mapping must be performed in step_04. Some probes may map to multiple genes or no genes — document mapping ambiguity.
- GSE262419's concentration-response metadata may be complex (chemical × concentration × plate). Ensure the concentration field is correctly parsed for each sample, as it is a key M3 input.
- GSE247541's small sample size (~15) means limited statistical power, but it adds a unique drug class (antivirals).
