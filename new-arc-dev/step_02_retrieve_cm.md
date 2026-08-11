# Step 2 — Data Retrieval: Cardiomyocyte Datasets

## Goal

Download all 4 cardiomyocyte drug-perturbation datasets from GEO and extract processed expression matrices plus experiment-design metadata.

## Dependencies

- **Step 1 (Data Inventory):** Provides the validated GEO accessions, expected sample counts, and dataset metadata needed for retrieval.

## Inputs

- Validated GEO accessions from step_01:
  - DToxS primary CM: GSE146096, GSE146097
  - DToxS iPSC-CM: GSE217421, GSE217423, GSE217424
  - GSE221595 (8 TKIs)
  - Rosell-Hidalgo 2023 (accession from step_01)
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

**GSE221595 (8 TKIs):**
- High-throughput RNA-seq of hiPSC-CM treated with 8 FDA-approved TKIs.
- Multi-dose, multi-timepoint design (day 1, 3, 5 for crizotinib and ponatinib; other TKIs at varying doses).
- Retrieve processed counts + metadata.
- Extract metadata fields: drug name, dose, timepoint, batch.

**Rosell-Hidalgo 2023:**
- Use the accession resolved in step_01.
- ScreenSeq RNA-seq of hiPSC-CM treated with 33 cardiotoxicants + 9 non-cardiotoxic compounds.
- Retrieve processed data + metadata.
- Extract metadata fields: drug name, cardiotoxicity class (cardiotoxicant vs non-cardiotoxic), dose, timepoint.

### 2.3 Metadata extraction

For every retrieved dataset, extract and standardize these metadata fields:
- **drug_name:** raw drug name as recorded in the GEO sample title/characteristics
- **dose:** concentration value + unit (if available)
- **timepoint:** treatment duration
- **donor_id:** donor or cell line identifier
- **cell_type:** primary cardiomyocyte or iPSC-derived cardiomyocyte
- **platform:** sequencing platform / library prep method
- **study:** dataset/study identifier (DToxS-primary, DToxS-iPSC, GSE221595, Rosell-Hidalgo)
- **batch:** batch or plate identifier (if available)
- **control_flag:** whether the sample is a vehicle/DMSO control

### 2.4 Control sample identification

For each dataset, identify which samples are vehicle/DMSO controls. These are needed for Δ computation in step_04. Controls may be labeled as "DMSO", "vehicle", "untreated", or dose=0. Document the control identification logic per dataset.

## Unit tests

1. **Download completeness:** For each dataset, verify the number of downloaded samples matches the expected count from step_01 (within tolerance). Flag any dataset with >10% sample count discrepancy.
2. **Metadata field presence:** For each dataset, verify all 9 metadata fields (drug_name, dose, timepoint, donor_id, cell_type, platform, study, batch, control_flag) are extractable. Fields that are unavailable in the source data should be documented as "not available" rather than missing.
3. **Control identification:** For each dataset, verify at least one control sample is identified. If no controls are found, flag the dataset — Δ computation will fail without controls.
4. **Drug name extraction:** For each dataset, verify the number of unique drug names extracted matches the expected drug count from step_01. Flag discrepancies.
5. **Expression matrix integrity:** For each dataset, verify the expression matrix has no all-zero rows (genes) or all-zero columns (samples). Verify matrix dimensions match sample count × gene count.

## Integration test

Validate that all 4 cardiomyocyte datasets are retrieved and metadata is consistent:
- All 4 datasets have expression matrices with non-zero dimensions.
- All 4 datasets have metadata tables with the required fields.
- The union of unique drug names across all 4 datasets is reported and cross-referenced with step_01's expected drug counts.
- Control samples are identified in all 4 datasets.
- No dataset has a sample count discrepancy >10% from step_01 expectations.
- The total combined sample count across all 4 datasets is reported.

## Acceptance criteria

- [ ] All 4 cardiomyocyte datasets downloaded from GEO
- [ ] Expression matrices extracted for each dataset
- [ ] Metadata tables with 9 standardized fields for each dataset
- [ ] Control samples identified in each dataset
- [ ] Drug names extracted and unique drug count per dataset matches step_01
- [ ] Total combined sample count reported

## Risks & notes

- GEO supplementary files can be large and slow to download. Budget ~2–4h for all 4 datasets.
- Some datasets may store expression data as supplementary files (not in the series matrix), requiring FTP download and manual parsing.
- DToxS datasets span multiple GSE accessions that must be combined — verify they don't have overlapping samples.
- Rosell-Hidalgo accession may not be in GEO at all — if so, this dataset must be requested from authors or excluded, reducing the non-cardiotoxic drug pool.
- Metadata formats vary wildly across GEO series — drug names may be embedded in sample titles, characteristics columns, or supplementary metadata files. Expect per-dataset parsing logic.
