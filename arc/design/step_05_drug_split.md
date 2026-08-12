# Step 5 — Drug Overlap Analysis & Unseen-Drug Split

## Goal

Compute the drug/chemical overlap matrix across all 5 cardiomyocyte datasets, canonicalize drug SMILES, and construct the unseen-drug train/test split that the model will be evaluated on.

## Dependencies

- **Step 1 (Data Inventory):** Provides the expected drug lists and drug class annotations per dataset.
- **Step 4 (Harmonize):** Provides the canonicalized drug SMILES and the harmonized cardiomyocyte data with covariate labels.

## Inputs

- Harmonized cardiomyocyte data with drug SMILES and covariate labels (from step_04)
- Drug class annotations from step_01 (KIs, anthracyclines, cardiac-acting, non-cardiac-acting)
- L1000 drug SMILES set (from step_04, for cross-referencing pretrain vs. eval drugs)

## Outputs

- **Drug overlap matrix:** which drugs appear in which cardiomyocyte datasets
- **Canonical SMILES registry:** every unique drug across all datasets with its canonical SMILES and dataset sources
- **Unseen-drug split:** train_drugs list and test_drugs list, with class balance and scaffold diversity annotations
- **Split report:** per-class drug counts in train vs. test, scaffold diversity metrics, leakage check results

## Behavior

### 5.1 SMILES canonicalization

For every drug/chemical across all 5 cardiomyocyte datasets:
1. Canonicalize SMILES using RDKit (`MolToSmiles(MolFromSmiles(smiles))`).
2. Compute InChIKey for unambiguous cross-dataset drug identity.
3. Merge drugs that share the same InChIKey (same drug, different name spellings across datasets).
4. Report the number of drugs/chemicals before and after canonicalization merge.
5. **SMILES-only filter (per user decision):** Only chemicals with resolvable SMILES (from step_04) enter the split. Chemicals without SMILES are excluded entirely — document the exclusion count, broken down by dataset and chemical class. For GSE262419, this may exclude a significant fraction of environmental chemicals.

### 5.2 Drug overlap matrix

Compute which drugs/chemicals appear in which cardiomyocyte datasets:
- For each unique drug/chemical (by InChIKey), record which datasets contain it.
- Drugs present in multiple datasets = shared training signal.
- Drugs unique to one dataset = candidates for held-out evaluation.
- Cross-reference with L1000: which cardiomyocyte drugs also appear in L1000 (these get pretrained drug embeddings; held-out drugs that are in L1000 still get pretrained embeddings, which is the intended transfer mechanism).
- The overlap matrix now spans 5 datasets. GSE262419's ~464 chemicals (subject to SMILES filtering) will have minimal overlap with the drug-focused datasets (DToxS, Rosell-Hidalgo) — most GSE262419 chemicals are unique to that dataset. GSE247541's 4 antivirals are also likely unique (no antivirals in other datasets). Document the overlap structure: which chemicals are shared vs. dataset-unique.

### 5.3 Unseen-drug split construction

**Split strategy:** Hold out entire drugs from the cardiomyocyte training set. The model must predict the cardiomyocyte transcriptomic response to drugs it has never seen (only their molecular structure).

**Split construction rules:**
1. **Stratify by drug class:** Ensure the held-out set contains drugs from each class: KIs, anthracyclines, cardiac-acting, non-cardiac-acting, antiviral (GSE247541), environmental chemical (GSE262419 chemicals without significant transcriptomic effect).
2. **Non-cardiac-acting priority:** The non-cardiac-acting class is the key test for non-cardiotoxic drug prediction. With GSE262419's ~395 chemicals without significant transcriptomic effect (subject to SMILES resolution), the non-cardiotoxic pool is now large enough for a robust hold-out (e.g., 20–30% of ~420 = ~84–126 held-out chemicals). Hold out a meaningful fraction.
3. **Scaffold diversity:** Compute Bemis-Murcko scaffolds for all drugs/chemicals. Ensure held-out drugs have scaffold diversity (not just scaffold variants of training drugs). If a held-out drug shares its scaffold with a training drug, the unseen-drug test is easier (scaffold leakage). GSE262419 environmental chemicals may introduce novel scaffolds not seen in drug datasets — this is a benefit for scaffold diversity but may make the unseen-drug test harder.
4. **No drug leakage:** Verify that no held-out drug/chemical (by InChIKey) appears in the training set. This is a hard constraint.
5. **L1000 overlap note:** Held-out drugs may appear in L1000 (pretrain). This is NOT leakage — L1000 is pretrain-only and never enters cardiomyocyte evaluation. The pretrained drug encoder is the transfer mechanism for unseen drugs.
6. **Environmental chemical class:** GSE262419 chemicals span 12 chemical classes. For stratification purposes, group them as "environmental chemical" (all GSE262419 chemicals without significant transcriptomic effect) vs. "drug" (GSE262419 chemicals that are FDA drugs or have drug-like properties). The split should hold out chemicals from both sub-groups.

### 5.4 Split report

Generate a report with:
- Per-class drug/chemical counts in train vs. test (classes: KIs, anthracyclines, cardiac-acting, non-cardiac-acting, antiviral, environmental chemical)
- Scaffold diversity: number of unique Bemis-Murcko scaffolds in train vs. test, and the number of test drugs whose scaffold also appears in train
- L1000 overlap: fraction of held-out drugs that have L1000 pretrain embeddings
- Sample counts: total samples in train vs. test (per dataset and overall). Note that GSE262419 concentration-response samples inflate per-chemical sample counts — report both per-chemical and per-sample counts.
- SMILES exclusion summary: number of chemicals excluded due to missing SMILES, broken down by dataset and class

## Unit tests

1. **SMILES canonicalization:** For a set of known drug name variants (e.g., "imatinib" vs "Gleevec" vs "STI-571"), verify they all canonicalize to the same InChIKey.
2. **No drug leakage:** Verify the intersection of train_drugs and test_drugs (by InChIKey) is empty. This is a hard-fail test.
3. **Class balance:** Verify the test set contains at least one drug/chemical from each class present in the full dataset. Report per-class train/test counts. Classes include: KIs, anthracyclines, cardiac-acting, non-cardiac-acting, antiviral, environmental chemical.
4. **Scaffold diversity:** Verify the number of test drugs whose Bemis-Murcko scaffold also appears in train is reported. Flag if >50% of test drugs share scaffolds with train (scaffold leakage risk).
5. **Non-cardiac-acting coverage:** Verify the test set contains at least 10 non-cardiac-acting chemicals (the expanded pool from GSE262419 now supports this). If the non-cardiac-acting pool is still too small after SMILES filtering, document this as a limitation.
6. **Sample count sanity:** Verify the train/test sample split ratio is approximately 80/20 (or the chosen ratio). Flag if any dataset is entirely in train or entirely in test. Note: GSE262419's concentration-response design means per-chemical sample counts are inflated — verify the split is per-chemical (not per-sample).
7. **SMILES-only filter:** Verify no chemical without SMILES appears in either train or test set. Verify the exclusion count is reported.

## Integration test

Validate the split is evaluation-ready:
- The split is deterministic (same input → same split) — verify by running twice and checking identical output.
- No drug/chemical leakage (hard constraint).
- Every drug/chemical in train and test has a canonical SMILES (from step_04). No SMILES-missing chemicals in either set.
- The split report contains all required fields (per-class counts including environmental and antiviral, scaffold diversity, L1000 overlap, sample counts, SMILES exclusion summary).
- The non-cardiac-acting class has robust representation in the test set (≥10 chemicals, or a documented limitation if SMILES filtering reduced the pool).
- The split is per-chemical (not per-sample) — verify that all concentration-response samples for a held-out chemical are in the test set, and all samples for a training chemical are in the train set.

## Acceptance criteria

- [ ] All drugs/chemicals canonicalized to InChIKey; duplicate names merged
- [ ] Drug/chemical overlap matrix computed across all 5 cardiomyocyte datasets
- [ ] Unseen-drug split constructed with class stratification (6 classes including environmental and antiviral)
- [ ] No drug/chemical leakage verified (hard constraint)
- [ ] No SMILES-missing chemicals in train or test (per user decision)
- [ ] Scaffold diversity reported; scaffold leakage risk flagged if >50%
- [ ] Non-cardiac-acting chemicals robustly represented in test set (≥10, or limitation documented)
- [ ] Split report generated with all required fields including SMILES exclusion summary

## Risks & notes

- The non-cardiotoxic pool is now large (~420 potential chemicals) but dominated by GSE262419 environmental chemicals. These differ fundamentally from clinical drugs in structure and exposure context. The unseen-drug test on environmental chemicals tests structural generalization, not clinical non-cardiotoxicity prediction.
- Scaffold leakage is a concern: many kinase inhibitors share scaffolds (e.g., quinazoline, pyrimidine). If held-out KIs share scaffolds with training KIs, the unseen-drug test is easier than truly novel scaffolds. GSE262419 environmental chemicals may introduce novel scaffolds, partially mitigating this.
- Drugs unique to one dataset (e.g., GSE247541's 4 antivirals, GSE262419's environmental chemicals) are natural held-out candidates, but holding out an entire dataset's unique drugs means that dataset contributes only training samples for shared drugs, reducing its training value.
- GSE262419's concentration-response design means each chemical has multiple samples. The split is per-chemical — all concentration points for a held-out chemical go to test. This means the test set may have many samples per held-out chemical, while the train set has many samples per training chemical. Report both per-chemical and per-sample split ratios.
- The split must be deterministic and saved to a file so that all downstream steps (training, evaluation) use the same split.
- SMILES filtering may exclude a significant fraction of GSE262419 chemicals. The effective chemical count for splitting may be much smaller than 464. Report the post-filtering chemical count before constructing the split.
