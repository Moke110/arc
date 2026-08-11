# Step 3 — Data Retrieval: L1000 Pretrain Corpus

## Goal

Retrieve and prepare the LINCS L1000 drug-perturbation corpus for use as the drug-encoder pretraining source (Phase 1).

## Dependencies

- **Step 1 (Data Inventory):** Provides the L1000 source locations (datalake mount + GEO accessions GSE92742, GSE70138) and expected profile counts.

## Inputs

- LINCS L1000 datalake mount (`/mnt/datalake/LINCS1000/`)
  - `RNAseq_transcriptomics_genesets/single_drug_perturbations-v1.0.gmt` — 875 drug signatures
  - `RNAseq_transcriptomics_signatures/human_geo_sigs.tsv` — 1.6GB signature file
- GEO accessions: GSE92742 (Level 3 expression), GSE70138 (Level 3 expression, supplementary)
- LINCS L1000 landmark gene list (978 genes)

## Outputs

- L1000 Δ-profile matrix: drug-perturbed expression change (Δ) per sample, restricted to drug perturbations (exclude genetic perturbations)
- L1000 metadata table: drug name, dose, timepoint, cell line, batch for each profile
- Drug name → SMILES mapping for all L1000 drugs (via PubChem/ChEMBL cross-reference)
- L1000 cell-line covariate labels for the covariate embedder (M4) during pretraining

## Behavior

### 3.1 Source selection

Two retrieval paths, use whichever is more complete:

**Path A — Datalake signatures (quick reference):**
- Load the 875 drug signatures from the GMT file. These are pre-computed up/down gene sets, not full profiles.
- Load the human_geo_sigs.tsv for additional signature context.
- Use these for quick drug-coverage assessment, NOT for pretraining (they are aggregated signatures, not per-sample profiles).

**Path B — GEO Level 3/5 data (pretraining source):**
- Download GSE92742 (main L1000 batch) and GSE70138 (supplementary batch) from GEO.
- Level 3 data contains per-sample expression on 978 landmark genes (after L1000 deconvolution).
- Level 5 data contains pre-computed differential expression signatures (perturbed vs. control).
- For pretraining, use Level 3 data to compute per-sample Δ = treated − matched_control.

### 3.2 Drug perturbation filtering

L1000 contains both drug (chemical) and genetic (knockdown/overexpression) perturbations. Filter to drug perturbations only:
- Exclude perturbation type = "trt_sh" (shRNA knockdown) or "trt_oe" (overexpression).
- Keep perturbation type = "trt_cp" (chemical compound).
- Expected: ~20,000 compounds, ~1.3M drug-perturbation profiles across ~70 cell lines.

### 3.3 Δ-profile construction

For each drug-perturbed sample, compute Δ = treated_expression − control_expression:
- Match each treated sample to its appropriate control (same cell line, same batch, same timepoint, DMSO/vehicle).
- If no matched control exists, use the cell-line mean DMSO profile as the control baseline.
- Δ is computed on the 978 landmark gene space.
- Store as a matrix: samples × 978 genes.

### 3.4 Drug name → SMILES mapping

For every unique drug in the L1000 drug-perturbation set:
- Look up canonical SMILES via PubChem (by compound name or CID).
- For drugs not found in PubChem, try ChEMBL cross-reference.
- For drugs with no resolvable SMILES, flag and exclude from pretraining (the drug encoder requires molecular structure).
- Report the fraction of L1000 drugs with successfully mapped SMILES.

### 3.5 Cell-line covariate labels

Extract cell-line identifiers for each profile. These become the covariate labels for M4 (covariate embedder) during L1000 pretraining:
- Core cell lines: A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3
- Additional cell lines in the full L1000 set (~70 total)
- These covariate embeddings are discarded in Phase 2 — only the M4 architecture carries over.

## Unit tests

1. **Perturbation type filter:** Verify all profiles in the output have perturbation type = "trt_cp" (chemical). No genetic perturbations present.
2. **Landmark gene count:** Verify the Δ-profile matrix has exactly 978 columns (landmark genes). No missing or extra genes.
3. **SMILES coverage:** Verify the fraction of drugs with mapped SMILES is reported. Flag if <80% of drugs have SMILES (would significantly reduce pretraining coverage).
4. **Control matching:** For a random sample of 100 treated profiles, verify each has a matched control (same cell line, same batch). Report the fraction with matched controls vs. cell-line-mean fallback.
5. **Δ non-triviality:** Verify the Δ-profile matrix is not all zeros. Report the fraction of non-zero entries. Flag if <5% non-zero (would indicate a computation error).
6. **Cell-line label presence:** Verify every profile has a non-null cell-line label. Report the number of unique cell lines.

## Integration test

Validate the L1000 corpus is pretraining-ready:
- The Δ-profile matrix has the expected dimensions (~1.3M samples × 978 genes, or a documented subset if downsampling is needed for compute).
- The metadata table has drug name, dose, timepoint, cell line, and batch for every profile.
- Every drug in the metadata table either has a SMILES string or is flagged for exclusion.
- The cell-line covariate labels cover all profiles.
- No genetic perturbations are present.
- The total drug count with SMILES is reported and is >10,000 (L1000's value is broad chemical-space coverage).

## Acceptance criteria

- [ ] L1000 drug-perturbation profiles retrieved (Level 3 from GEO or datalake)
- [ ] Genetic perturbations filtered out
- [ ] Δ-profile matrix constructed (treated − matched control on 978 landmark genes)
- [ ] Drug name → SMILES mapping completed for >80% of drugs
- [ ] Cell-line covariate labels extracted for all profiles
- [ ] Metadata table complete (drug, dose, timepoint, cell line, batch)

## Risks & notes

- L1000 Level 3 data is large (~50–100 GB). Download may take several hours. Consider using the datalake-mounted version if it contains per-sample profiles.
- L1000 never enters cardiomyocyte evaluation — it is pretrain-only. This is a hard constraint from the user's decision.
- Some L1000 drugs may not have resolvable SMILES (especially older or withdrawn compounds). The drug encoder can only encode drugs with molecular structure.
- The ~70 L1000 cell lines are all cancer cell lines — drug mechanisms may differ from cardiomyocytes. The pretrain gives chemical-structure-aware embeddings, not cardiomyocyte-specific drug responses. Phase 2 fine-tuning (step_11) adapts the mapping.
- If GEO download is infeasible due to size, consider using the datalake signatures (875 drugs) as a reduced pretraining set, with the tradeoff of less chemical-space coverage.
