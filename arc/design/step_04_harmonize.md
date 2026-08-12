# Step 4 — Data Harmonization & Preprocessing

## Goal

Harmonize all retrieved datasets (6 cardiomyocyte + L1000) into a shared gene space, normalize expression, compute Δ profiles, and assign covariate labels — producing training-ready data matrices for both pretraining and cardiomyocyte training.

## Dependencies

- **Step 2 (Retrieve CM):** Provides raw cardiomyocyte expression matrices and metadata.
- **Step 3 (Retrieve L1000):** Provides L1000 Δ-profiles, metadata, and SMILES mappings.

## Inputs

- 6 cardiomyocyte dataset expression matrices with metadata (from step_02)
- L1000 Δ-profile matrix with metadata and SMILES (from step_03)
- L1000 landmark gene list (978 genes)
- Reference gene mapping (HGNC symbols / Ensembl IDs)
- TempO-Seq probe-to-gene mapping (for GSE262419)

## Outputs

- **L1000 pretrain set:** normalized Δ-profile matrix on 978 landmark genes + covariate labels (cell line, batch) + drug SMILES + dose
- **Cardiomyocyte training set:** normalized expression matrices on shared gene space (~10,000 genes) + Δ-profiles + covariate labels (cell_type, platform, study, donor) + drug SMILES + dose
- **Shared gene space definition:** the gene list used for cardiomyocyte training (L1000 landmark genes ∪ top highly-variable genes across cardiomyocyte datasets)
- **Drug name → SMILES mapping** for all cardiomyocyte drugs (via PubChem/ChEMBL)

## Behavior

### 4.1 Gene harmonization

Map all datasets to a shared gene space:
1. **Gene ID standardization:** Convert all gene identifiers (probe IDs, Ensembl IDs, gene symbols) to HGNC gene symbols. Use a reference mapping (e.g., org.Hs.eg.db or mygene.info). For GSE262419 (TempO-Seq), map probe IDs to gene symbols using the platform's probe annotation. Document probes that map to multiple genes (ambiguous) or no genes (unmapped) — exclude both from the shared gene space.
2. **L1000 pretrain gene space:** Use the 978 L1000 landmark genes. The L1000 Δ-profiles are already on this space.
3. **Cardiomyocyte training gene space:** Define as the union of:
   - 978 L1000 landmark genes (for compatibility with pretrained decoder)
    - Top N highly-variable genes across the 5 cardiomyocyte datasets (N chosen so total ~10,000)
4. **Gene intersection:** For cardiomyocyte datasets, intersect each dataset's genes with the shared gene space. Document genes lost in the intersection.
5. **TempO-Seq (GSE262419) special handling:** TempO-Seq targets ~20,000 genes via predetermined probes. Intersect the TempO-Seq gene set with the shared gene space defined above. Document:
   - Genes in TempO-Seq but NOT in the shared gene space (TempO-Seq-unique genes lost)
   - Genes in the shared gene space but NOT in TempO-Seq (genes missing from TempO-Seq — these will have missing values for GSE262419 samples)
   - The final shared gene space size after intersection (may be smaller than ~10,000 if TempO-Seq probe coverage doesn't overlap well with RNA-seq-discovered HVGs)
    - If the shared gene space shrinks significantly (>20% reduction), consider re-computing HVGs using only the gene set common to all 5 datasets as a fallback.

### 4.2 Normalization

Within each dataset (do NOT cross-dataset normalize):
1. **Library-size normalize:** Divide each sample's expression by total counts, scale to a common library size (e.g., 1e6).
2. **Log1p transform:** Apply log1p to normalized counts.
3. **No cross-dataset batch correction at this stage** — batch/platform effects are handled by the model's MMD regularization (M1) and covariate embeddings (M4) during training, not by preprocessing.

### 4.3 Δ-profile computation (cardiomyocyte)

For each cardiomyocyte drug-treated sample, compute Δ = treated − control:
1. **Matched control:** For each treated sample, find the matched control (same donor, same batch, same timepoint, vehicle/DMSO). If the dataset is paired (control and treated from same donor/batch), use the paired control.
2. **Unpaired fallback:** If no matched control exists, use the donor-mean or batch-mean control expression.
3. **Δ = log1p(treated) − log1p(control).** This is the prediction target for the model.
4. Store both the control expression (input to M1) and the Δ (target for M5) for each sample.
5. **Concentration-response handling (GSE262419):** Each concentration point is a separate sample with its own Δ. Match each treated sample to the concurrent DMSO/vehicle control from the same concentration plate (per-plate control matching from step_02). The concentration value is stored as the dose input for M3 — the model design is unchanged. One chemical with N concentrations contributes N separate (control, Δ, concentration) training samples. This inflates sample count but not chemical diversity.

### 4.4 Drug name → SMILES mapping (cardiomyocyte)

For every unique drug/chemical across the 6 cardiomyocyte datasets:
- Look up canonical SMILES via PubChem (by compound name or CID).
- For drugs not found in PubChem, try ChEMBL cross-reference.
- For drugs with no resolvable SMILES, flag and exclude from training (the drug encoder requires molecular structure).
- Report the mapping success rate per dataset.
- Cross-reference with L1000 SMILES mappings to identify drugs present in both L1000 and cardiomyocyte datasets.
- **GSE262419 special handling:** The 464 chemicals include environmental chemicals, industrial compounds, and food constituents that may have lower SMILES resolution rates than FDA drugs. Per user decision, only chemicals with resolvable SMILES enter training. Report:
  - SMILES resolution rate for GSE262419 (expected to be lower than drug-only datasets)
  - Number of chemicals excluded due to missing SMILES
  - Chemical class breakdown of resolved vs. unresolved chemicals (e.g., environmental chemicals may have lower coverage than pharmaceutical drugs)
  - Some environmental chemicals (e.g., PFAS, industrial compounds) may have PubChem entries but non-standard SMILES — verify that RDKit can parse the resolved SMILES.

### 4.5 Covariate label assignment

For each cardiomyocyte sample, assign the covariate tuple: {cell_type, platform, study, donor}:
- **cell_type:** primary cardiomyocyte or iPSC-derived cardiomyocyte
- **platform:** 3'-DGE RNA-seq, bulk RNA-seq, ScreenSeq, TempO-Seq, RNA-seq (5 categories)
- **study:** DToxS-primary, DToxS-iPSC, Rosell-Hidalgo, GSE262419, GSE247541 (5 categories)
- **donor:** donor or cell line identifier

For each L1000 sample, assign: {cell_line, batch} (these are the pretrain covariates, discarded in Phase 2).

### 4.6 Data quality checks

For each dataset:
- Report the number of samples with matched controls vs. fallback controls.
- Report the number of genes lost in the shared gene space intersection.
- Report the number of drugs/chemicals with/without SMILES.
- Flag any dataset with >20% missing SMILES or >50% fallback controls.
- For GSE262419: additionally report (a) TempO-Seq probe-to-gene mapping success rate, (b) the number of shared gene space genes missing from TempO-Seq, (c) the final shared gene space size after TempO-Seq intersection, (d) SMILES resolution rate broken down by chemical class (12 categories).

## Unit tests

1. **Gene space consistency:** Verify the cardiomyocyte training gene space contains all 978 L1000 landmark genes. Verify the total gene count is ~10,000 (±2,000 — wider tolerance due to TempO-Seq intersection). Document the final gene count after TempO-Seq intersection.
2. **Normalization correctness:** For a random sample, verify post-normalization expression values are log1p-scaled and non-negative. Verify library sizes are approximately equal across samples within a dataset.
3. **Δ computation correctness:** For a paired sample (treated + matched control), verify Δ = treated − control element-wise. Verify Δ values can be positive or negative (not all zero or all positive). For GSE262419 concentration-response samples, verify each concentration point has its own Δ computed against the correct per-plate control.
4. **SMILES mapping:** For each cardiomyocyte dataset, verify the fraction of drugs/chemicals with SMILES is reported. Flag any dataset with <80% SMILES coverage. For GSE262419, report SMILES coverage separately and do NOT flag if <80% (environmental chemicals are expected to have lower coverage — per user decision, chemicals without SMILES are simply excluded).
5. **Covariate completeness:** For each cardiomyocyte sample, verify all 4 covariate fields (cell_type, platform, study, donor) are non-null. Verify platform has 5 categories and study has 5 categories. For each L1000 sample, verify cell_line and batch are non-null.
6. **No cross-dataset normalization:** Verify that normalization statistics (library size, gene means) are computed within-dataset, not across datasets. Check that a sample's normalized expression depends only on its own dataset's statistics.
7. **TempO-Seq probe mapping:** For GSE262419, verify that probe IDs have been mapped to gene symbols. Verify the fraction of probes with successful gene mapping is reported. Verify ambiguous and unmapped probes are excluded.

## Integration test

Validate the harmonized data is training-ready:
- All 5 cardiomyocyte datasets are on the same shared gene space (final count documented after TempO-Seq intersection).
- L1000 pretrain data is on 978 landmark genes.
- Every cardiomyocyte training sample has: control expression, Δ profile, drug/chemical SMILES (or is excluded for missing SMILES), dose/concentration, and 4 covariate labels.
- Every L1000 pretrain sample has: Δ profile, drug SMILES, dose, and 2 covariate labels.
- The drug SMILES sets for L1000 and cardiomyocyte datasets are cross-referenced (overlap reported).
- No dataset has >20% missing SMILES or >50% fallback controls (or these are documented as accepted risks). GSE262419 may have >20% missing SMILES — this is expected and accepted per user decision (chemicals without SMILES are excluded).
- GSE262419 concentration-response samples each have a concentration value and a per-plate matched control Δ.

## Acceptance criteria

- [ ] All 5 datasets mapped to shared gene space (final count documented after TempO-Seq intersection)
- [ ] Within-dataset normalization applied (library-size + log1p)
- [ ] Δ profiles computed for all cardiomyocyte drug/chemical-treated samples (including GSE262419 concentration-response points)
- [ ] Drug/chemical name → SMILES mapping completed; GSE262419 SMILES resolution rate reported by chemical class
- [ ] Chemicals without SMILES excluded from training (per user decision)
- [ ] Covariate labels assigned to all samples (5 platforms, 5 studies)
- [ ] Data quality report generated (matched vs. fallback controls, gene loss, SMILES coverage, TempO-Seq probe mapping)

## Risks & notes

- Cross-platform gene-space intersection may lose genes — especially for 3'-DGE (which captures only 3' ends) vs. bulk RNA-seq (full transcript). The shared gene space may be smaller than expected if 3'-DGE gene coverage is limited.
- TempO-Seq (GSE262419) adds a TempO-Seq platform with ~20,000 targeted genes. The intersection of TempO-Seq probes with the RNA-seq-discovered HVG gene space may shrink the shared gene space. If the reduction is >20%, consider re-computing HVGs using only genes common to all 5 datasets.
- Some cardiomyocyte drugs may not have resolvable SMILES (especially monoclonal antibodies in DToxS-iPSC — mAbs don't have small-molecule SMILES). These drugs cannot be encoded by the molecular drug encoder and must be excluded or handled separately.
- GSE262419 environmental chemicals may have lower SMILES resolution rates than FDA drugs. Per user decision, chemicals without SMILES are excluded — the effective training set from GSE262419 may be smaller than 464 chemicals.
- Δ computation depends on quality control matching. Poor control matching (e.g., using batch-mean controls when paired controls exist but weren't labeled) will introduce noise into the Δ target. For GSE262419, per-plate control matching is critical — mismatched controls across concentration plates would introduce plate-level batch effects into the Δ.
- The decision NOT to apply cross-dataset batch correction at preprocessing is deliberate — batch effects are handled by the model (MMD + covariate embeddings). This preserves biological signal that might be removed by aggressive batch correction.
- GSE262419's concentration-response design means one chemical contributes multiple Δ profiles (one per concentration). This is not a data quality issue — it is additional training signal. The concentration value is passed to M3 as the dose input, so the model can learn concentration-dependent effects.
