# Step 1 — Data Inventory & Accession Validation

## Goal

Catalog all drug-perturbed transcriptome datasets relevant to cardiomyocytes (with >3 distinct drugs) and the L1000 pretrain corpus, validate GEO accessions, and produce a structured inventory with dataset metadata for downstream retrieval.

## Dependencies

- **Step 0 (Review):** The architecture decision determines what data is needed (cardiomyocyte training data + L1000 pretrain corpus). The >3-drug filter and dataset tiering follow from the review's data requirements.

## Inputs

- GEO database (queryable via geoquery / direct FTP)
- rummageo index at the datalake (`human-gse-processed-meta.json`, 16,920 human GEO series)
- LINCS L1000 datalake mount
- Published literature references for dataset discovery (van Hasselt 2020, Hansen 2024, Rosell-Hidalgo 2023)

## Outputs

- A validated dataset inventory table with: GEO accession(s), cell type, platform, drug count, drug list, sample count, timepoint, donor count, tier (1/2/3), and role (train/eval/pretrain)
- A non-cardiotoxic drug pool list (~25 drugs across datasets)
- Accession resolution for any datasets with uncertain GEO IDs (especially Rosell-Hidalgo 2023)

## Behavior

### 1.1 Cardiomyocyte drug-perturbation datasets (>3-drug filter)

Apply the user's filter: only include cardiomyocyte datasets with MORE THAN 3 distinct drugs. After filtering, 4 datasets remain:

| Dataset | GEO Accession(s) | Cell type | Platform | Drugs | ~Samples | Timepoint | Reference | Tier |
|---|---|---|---|---|---|---|---|---|
| DToxS primary CM | GSE146096, GSE146097 | Primary human cardiomyocytes (adult heart ventricle, PromoCell, 4 donors) | 3'-DGE RNA-seq | ~26 kinase inhibitors | ~502 | 48h | van Hasselt 2020 Nature Comms | 1 |
| DToxS iPSC-CM | GSE217421, GSE217423, GSE217424 | iPSC-derived cardiomyocytes (6 donors) | Bulk RNA-seq | 54 FDA drugs (23 KIs, 4 mAbs, 4 anthracyclines, 7 cardiac-acting, 16 non-cardiac-acting) | ~1171 | 48h | Hansen 2024 Nature Comms | 1 |
| GSE221595 | GSE221595 | hiPSC-CM | High-throughput RNA-seq | 8 TKIs (multi-dose, multi-time: day 1/3/5) | verify | day 1/3/5 | Rosell-Hidalgo 2023 (verify) | 2 |
| Rosell-Hidalgo 2023 | verify GEO accession | hiPSC-CM | RNA-seq (ScreenSeq) | 33 cardiotoxicants + 9 non-cardiotoxic compounds (42 total) | verify | verify | Rosell-Hidalgo 2023 Expert Opin Drug Metab Toxicol | 2 |

**Estimated total cardiomyocyte drug-perturbed samples:** ~1,700–2,000+ across ~90+ unique drugs (with overlap), 3 platforms, 2 cell types (primary CM, iPSC-CM), 10+ donors/batches.

### 1.2 Non-cardiotoxic drug pool

The non-cardiotoxic drug pool (~25 drugs) comes from:
- 16 non-cardiac-acting drugs in DToxS-iPSC (GSE217421/217423/217424)
- 9 non-cardiotoxic compounds in Rosell-Hidalgo 2023

This is the key data for the non-cardiotoxic drug prediction requirement.

### 1.3 Non-cardiomyocyte drug-perturbation data (pretrain only)

| Dataset | Source | Cell types | Drugs | Profiles | Role |
|---|---|---|---|---|---|
| LINCS L1000 | Datalake mount + GEO (GSE92742, GSE70138) | ~70 cancer cell lines (A375, A549, MCF7, PC3, HEPG2, HT29, HA1E, HCC515...) — no cardiomyocytes | ~20,000 compounds | ~1.3M profiles (Level 3/5) | Pretrain drug molecular encoder + drug-response latent space; NOT used for cardiomyocyte evaluation |
| LINCS L1000 signatures | Datalake (`RNAseq_transcriptomics_signatures/`) | Same | 875 drug up/down signatures (GMT) | 875 | Quick-look drug signature reference |

### 1.4 Accession validation

For each dataset in the inventory:
- Verify the GEO accession resolves to the expected series (correct title, organism, platform, sample count).
- For Rosell-Hidalgo 2023: search GEO for the accession (not yet confirmed). Search terms: "Rosell-Hidalgo cardiotoxicity ScreenSeq hiPSC cardiomyocyte".
- For GSE221595: verify sample count and drug list from GEO metadata.
- For DToxS accessions: verify both arms (GSE146096/146097 for primary CM, GSE217421/217423/217424 for iPSC-CM) resolve correctly.

### 1.5 Datasets excluded by the >3-drug filter

Document datasets that were found but excluded because they have ≤3 distinct drugs, for transparency:
- GSE213311 (doxorubicin + indisulam rescue, ~2 drugs)
- GSE232331 (doxorubicin only, 1 drug)
- GSE230638 (doxorubicin + metformin, 2 drugs)
- Burridge 2016 (doxorubicin dose-response, 1 drug)
- GSE113546 (statins, verify drug count)
- GSE129793 (calcium channel blockers, verify drug count)

## Unit tests

1. **>3-drug filter:** For each dataset in the final inventory, verify the drug count is >3. Any dataset with ≤3 drugs must be in the excluded list, not the inventory.
2. **Accession format:** All GEO accessions match the pattern GSE\d+. No malformed accession strings.
3. **Cell type validation:** Every cardiomyocyte dataset has cell type = "primary cardiomyocyte" or "iPSC-derived cardiomyocyte" (or variant). No non-cardiomyocyte datasets in the cardiomyocyte inventory.
4. **Non-cardiotoxic drug pool:** Verify the non-cardiotoxic drug pool list has entries from at least 2 datasets. Each drug has a source dataset annotation.
5. **L1000 exclusion:** Verify LINCS L1000 is NOT in the cardiomyocyte dataset inventory — it is in the pretrain-only section.
6. **Tier assignment:** Tier 1 datasets have verified accessions and >10 drugs. Tier 2 datasets have verified or pending accessions and >3 drugs.

## Integration test

Validate the complete inventory is internally consistent and retrieval-ready:
- Every dataset in the inventory has a resolvable GEO accession (or a documented pending-verification status).
- The total drug count across cardiomyocyte datasets (counting unique drugs) is reported.
- The non-cardiotoxic drug pool size is reported.
- No dataset appears in both the cardiomyocyte inventory and the pretrain-only section.
- The inventory table has no missing required fields (accession, cell type, platform, drug count, tier).

## Acceptance criteria

- [ ] 4 cardiomyocyte datasets cataloged with >3 drugs each
- [ ] Rosell-Hidalgo GEO accession resolved or documented as pending
- [ ] GSE221595 sample count verified from GEO
- [ ] Non-cardiotoxic drug pool listed (~25 drugs with source datasets)
- [ ] LINCS L1000 cataloged as pretrain-only with datalake path
- [ ] Excluded datasets (≤3 drugs) documented for transparency

## Risks & notes

- Rosell-Hidalgo 2023 GEO accession is not yet confirmed — this is the primary accession-resolution risk. If it cannot be found in GEO, the dataset may need to be requested from authors or excluded.
- GSE221595 sample count needs verification from GEO metadata — the exact number of samples affects the total training set size estimate.
- The non-cardiotoxic drug pool (~25 drugs) is small. This is the fundamental data gap for non-cardiotoxic drug prediction. The unseen-drug evaluation (step_12) will reveal whether molecular-structure-based generalization works for this class.
- Drug name → SMILES mapping (done in step_04) may reveal that some "non-cardiotoxic" drugs overlap with cardiotoxic drugs in other datasets, reducing the effective non-cardiotoxic pool.
