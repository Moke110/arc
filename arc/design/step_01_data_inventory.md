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

Apply the user's filter: only include cardiomyocyte datasets with MORE THAN 3 distinct drugs. GSE221595 was removed per user decision (2026-08-12) — see §1.5b. After filtering, 5 datasets remain:

| Dataset | GEO Accession(s) | Cell type | Platform | Drugs/Chemicals | ~Samples | Timepoint | Reference | Tier |
|---|---|---|---|---|---|---|---|---|
| DToxS primary CM | GSE146096, GSE146097 | Primary human cardiomyocytes (adult heart ventricle, PromoCell, 4 donors) | 3'-DGE RNA-seq | ~26 kinase inhibitors | ~502 | 48h | van Hasselt 2020 Nature Comms | 1 |
| DToxS iPSC-CM | GSE217421, GSE217423, GSE217424 | iPSC-derived cardiomyocytes (6 donors) | Bulk RNA-seq | 54 FDA drugs (23 KIs, 4 mAbs, 4 anthracyclines, 7 cardiac-acting, 16 non-cardiac-acting) | ~1171 | 48h | Hansen 2024 Nature Comms | 1 |
| Rosell-Hidalgo 2023 | GSE244740 | hiPSC-CM | RNA-seq (ScreenSeq) | 33 cardiotoxicants + 9 non-cardiotoxic compounds (42 total) | 3072 | 24h | Rosell-Hidalgo 2024 Expert Opin Drug Metab Toxicol | 2 |
| GSE262419 | GSE262419 | iPSC-derived cardiomyocytes | TempO-Seq (~20,000 targeted genes) | 464 chemicals (12 classes: drugs + environmental + food constituents); only SMILES-resolvable chemicals enter training | 4,558 (concentration-response design) | concentration-response (verify exact exposure time from paper) | Tsai 2024 Chem Res Toxicol (PMID 39046974) | 1 |
| GSE247541 | GSE247541 | iPSC-derived cardiomyocytes (hCM) | RNA-seq | 4 antiviral drugs (apilimod, remdesivir, ritonavir, lopinavir) | ~15 (3 DMSO + 4 drugs × 3 replicates) | verify | COVID-19 antiviral cardiotoxicity study | 2 |

**Estimated total cardiomyocyte drug-perturbed samples:** ~6,500+ across ~500+ unique chemicals/drugs (with overlap, subject to SMILES resolution), 5 platforms, 2 cell types (primary CM, iPSC-CM), 10+ donors/batches.

**Note on GSE244109:** GSE244109 (26 PFAS, 197 samples) is a SubSeries of the GSE262419 SuperSeries. Its data is fully contained within GSE262419. It is NOT added separately to avoid sample duplication.

### 1.2 Non-cardiotoxic drug pool

The non-cardiotoxic drug pool comes from:
- 16 non-cardiac-acting drugs in DToxS-iPSC (GSE217421/217423/217424)
- 9 non-cardiotoxic compounds in Rosell-Hidalgo 2023
- GSE262419 chemicals without significant transcriptomic effect: the Tsai 2024 paper reports 69 of 464 chemicals had significant transcriptomic effects, leaving ~395 chemicals as potential non-cardiotoxic/low-toxicity controls (subject to SMILES resolution). These are environmental chemicals and food constituents, not validated non-cardiotoxic drugs — the assumption "no significant transcriptomic effect" ≈ "non-cardiotoxic" is an approximation that must be validated.

This expanded pool (~25 + ~395 = ~420 potential non-cardiotoxic chemicals) is the key data for the non-cardiotoxic drug prediction requirement. However, the environmental chemicals differ fundamentally from clinical drugs in structure and exposure context.

### 1.3 Non-cardiomyocyte drug-perturbation data (pretrain only)

| Dataset | Source | Cell types | Drugs | Profiles | Role |
|---|---|---|---|---|---|
| LINCS L1000 | Datalake mount + GEO (GSE92742, GSE70138) | ~70 cancer cell lines (A375, A549, MCF7, PC3, HEPG2, HT29, HA1E, HCC515...) — no cardiomyocytes | ~20,000 compounds | ~1.3M profiles (Level 3/5) | Pretrain drug molecular encoder + drug-response latent space; NOT used for cardiomyocyte evaluation |
| LINCS L1000 signatures | Datalake (`RNAseq_transcriptomics_signatures/`) | Same | 875 drug up/down signatures (GMT) | 875 | Quick-look drug signature reference |

### 1.4 Accession validation

For each dataset in the inventory:
- Verify the GEO accession resolves to the expected series (correct title, organism, platform, sample count).
- For Rosell-Hidalgo 2023: accession resolved to **GSE244740** (3072 samples, GPL24676, ScreenSeq, 42 compounds × 8 concentrations × 3 replicates, 24 h) — verified 2026-08-12.
- For DToxS accessions: verify both arms (GSE146096/146097 for primary CM, GSE217421/217423/217424 for iPSC-CM) resolve correctly.
- For GSE262419: verify SuperSeries status, sample count (~4,558), platform (GPL16791 / TempO-Seq), and that GSE244109 is listed as a SubSeries. Supplementary files are 16 per-plate CSVs (`GSE262419_Plate01..16.csv.gz`) + `GSE262419_hash.csv.gz` (verified 2026-08-12).
- For GSE247541: verify sample count (~15), platform (RNA-seq), and drug list (apilimod, remdesivir, ritonavir, lopinavir).

### 1.5 Datasets excluded by the >3-drug filter

Document datasets that were found but excluded because they have ≤3 distinct drugs, for transparency:
- GSE213311 (doxorubicin + indisulam rescue, ~2 drugs)
- GSE232331 (doxorubicin only, 1 drug)
- GSE230638 (doxorubicin + metformin, 2 drugs)
- Burridge 2016 (doxorubicin dose-response, 1 drug)
- GSE113546 (statins, verify drug count)
- GSE129793 (calcium channel blockers, verify drug count)
- GSE244109 (26 PFAS, 197 samples) — NOT excluded by the >3-drug filter, but NOT added separately because it is a SubSeries of GSE262419. Its PFAS data is fully contained within the GSE262419 SuperSeries. Adding both would create duplicate samples.

### 1.5b Dataset removed by user decision

- **GSE221595 (8 TKIs, hiPSC-CM)** — removed from the pipeline per user decision on 2026-08-12. The data (`pert_data/gse221595/`) was deleted and all planning docs updated. Note: the GEO series matrix holds a single pooled pseudo-sample (1 sample; per-drug data lives in `GSE221595_RAW.tar`), so it contributed no per-drug replicate-level samples. Retained here only for provenance.

### 1.6 Concentration-time variable ranges

Record the concentration and timepoint ranges per dataset. This metadata is needed for step_04 (Δ computation with matched controls) and step_06 (M4 covariate embedder). The model design is unchanged — dose/concentration is already an M3 input, and each concentration point is a separate training sample with its own Δ.

| Dataset | Concentration range | # Concentrations per chemical | Timepoint(s) | Design type |
|---|---|---|---|---|
| DToxS primary CM | single dose (verify exact value) | 1 | 48h | single-dose |
| DToxS iPSC-CM | single dose (verify exact value) | 1 | 48h | single-dose |
| Rosell-Hidalgo 2023 | 8 concentrations (verify range from paper) | 8 | 24h | concentration-response |
| GSE262419 | concentration-response (verify range from paper) | multiple per chemical (drives 4,558 samples from ~464 chemicals) | verify from paper | concentration-response |
| GSE247541 | single dose (verify) | 1 | verify | single-dose |

**Key implication for GSE262419:** The concentration-response design means each chemical contributes multiple samples (one per concentration point). Each concentration point gets its own Δ computed against the matched DMSO/vehicle control for that concentration plate. The concentration value is passed to M3 as the dose input. This inflates sample count but not chemical diversity — the unseen-drug split (step_05) is per-chemical, so the number of unique chemicals is what matters for generalization assessment.

## Unit tests

1. **>3-drug filter:** For each dataset in the final inventory, verify the drug/chemical count is >3. Any dataset with ≤3 drugs must be in the excluded list, not the inventory.
2. **Accession format:** All GEO accessions match the pattern GSE\d+. No malformed accession strings.
3. **Cell type validation:** Every cardiomyocyte dataset has cell type = "primary cardiomyocyte" or "iPSC-derived cardiomyocyte" (or variant). No non-cardiomyocyte datasets in the cardiomyocyte inventory.
4. **Non-cardiotoxic drug pool:** Verify the non-cardiotoxic drug pool list has entries from at least 2 datasets. Each drug/chemical has a source dataset annotation. For GSE262419 entries, annotate as "environmental chemical — no significant transcriptomic effect" to distinguish from validated non-cardiotoxic drugs.
5. **L1000 exclusion:** Verify LINCS L1000 is NOT in the cardiomyocyte dataset inventory — it is in the pretrain-only section.
6. **Tier assignment:** Tier 1 datasets have verified accessions and >10 drugs/chemicals. Tier 2 datasets have verified or pending accessions and >3 drugs/chemicals.
7. **Sub-series dedup:** Verify GSE244109 does NOT appear as a separate dataset in the inventory. Verify GSE262419 is present and annotated as a SuperSeries containing GSE244109.
8. **Concentration-time ranges:** Verify section 1.6 has an entry for every dataset in the inventory. Verify GSE262419's entry notes concentration-response design with multiple concentrations per chemical.

## Integration test

Validate the complete inventory is internally consistent and retrieval-ready:
- Every dataset in the inventory has a resolvable GEO accession (or a documented pending-verification status).
- The total drug/chemical count across cardiomyocyte datasets (counting unique chemicals by InChIKey after SMILES canonicalization in step_04) is reported.
- The non-cardiotoxic drug pool size is reported, with a breakdown of validated non-cardiotoxic drugs vs. GSE262419 environmental chemicals without significant transcriptomic effect.
- No dataset appears in both the cardiomyocyte inventory and the pretrain-only section.
- The inventory table has no missing required fields (accession, cell type, platform, drug/chemical count, tier).
- GSE244109 is documented as a sub-series of GSE262419, not as a separate inventory entry.
- GSE221595 appears nowhere in the inventory (removed per user decision).
- The concentration-time variable ranges table (section 1.6) has entries for all 5 datasets.

## Acceptance criteria

- [ ] 5 cardiomyocyte datasets cataloged with >3 drugs/chemicals each
- [ ] Rosell-Hidalgo GEO accession resolved (GSE244740)
- [ ] GSE221595 removed from inventory and all planning docs
- [ ] GSE262419 verified as SuperSeries with ~4,558 samples, TempO-Seq platform, GSE244109 documented as sub-series
- [ ] GSE247541 verified with ~15 samples, 4 antiviral drugs
- [ ] Non-cardiotoxic drug pool listed (~25 validated drugs + ~395 GSE262419 environmental chemicals, subject to SMILES resolution)
- [ ] Concentration-time variable ranges recorded for all 5 datasets (section 1.6)
- [ ] LINCS L1000 cataloged as pretrain-only with datalake path
- [ ] Excluded datasets (≤3 drugs) documented for transparency

## Risks & notes

- Rosell-Hidalgo 2023 GEO accession is resolved (GSE244740, verified 2026-08-12); no author request needed.
- GSE221595 removed per user decision (2026-08-12) — no longer affects the total training set size estimate.
- The non-cardiotoxic drug pool is now large (~420 potential chemicals) but the GSE262419 contribution (~395 environmental chemicals) is fundamentally different from validated non-cardiotoxic drugs. The assumption "no significant transcriptomic effect" ≈ "non-cardiotoxic" is an approximation. These environmental chemicals may have cardiac effects not captured by transcriptomics, or may lack clinical relevance.
- GSE262419 uses TempO-Seq (~20,000 targeted genes), a probe-based platform fundamentally different from the RNA-seq platforms used by other datasets. Gene-space intersection in step_04 may shrink the shared gene space. TempO-Seq's predetermined probe coverage may miss genes captured by full-transcript RNA-seq.
- GSE262419's 464 chemicals include environmental chemicals and food constituents that may not have resolvable SMILES. Per user decision, chemicals without SMILES are excluded from training. The effective chemical count may be smaller than 464. SMILES resolution rate must be reported in step_04.
- GSE262419's concentration-response design inflates sample count (~4,558 samples from ~464 chemicals × multiple concentrations). This inflates per-chemical sample count but not chemical diversity. The unseen-drug split (step_05) is per-chemical, so the number of unique chemicals is what matters for generalization assessment.
- GSE262419's ~4,558 samples represent ~70% of total training data. Platform-balanced sampling in step_11 is critical to prevent TempO-Seq from dominating training.
- GSE247541 has only ~15 samples and 4 drugs — it is small but adds a unique drug class (antivirals) not present in other datasets.
- Drug name → SMILES mapping (done in step_04) may reveal that some "non-cardiotoxic" drugs overlap with cardiotoxic drugs in other datasets, reducing the effective non-cardiotoxic pool.
