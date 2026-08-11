# Step 5 — Drug Overlap Analysis & Unseen-Drug Split

## Goal

Compute the drug overlap matrix across all cardiomyocyte datasets, canonicalize drug SMILES, and construct the unseen-drug train/test split that the model will be evaluated on.

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

For every drug across all cardiomyocyte datasets:
1. Canonicalize SMILES using RDKit (`MolToSmiles(MolFromSmiles(smiles))`).
2. Compute InChIKey for unambiguous cross-dataset drug identity.
3. Merge drugs that share the same InChIKey (same drug, different name spellings across datasets).
4. Report the number of drugs before and after canonicalization merge.

### 5.2 Drug overlap matrix

Compute which drugs appear in which cardiomyocyte datasets:
- For each unique drug (by InChIKey), record which datasets contain it.
- Drugs present in multiple datasets = shared training signal.
- Drugs unique to one dataset = candidates for held-out evaluation.
- Cross-reference with L1000: which cardiomyocyte drugs also appear in L1000 (these get pretrained drug embeddings; held-out drugs that are in L1000 still get pretrained embeddings, which is the intended transfer mechanism).

### 5.3 Unseen-drug split construction

**Split strategy:** Hold out entire drugs from the cardiomyocyte training set. The model must predict the cardiomyocyte transcriptomic response to drugs it has never seen (only their molecular structure).

**Split construction rules:**
1. **Stratify by drug class:** Ensure the held-out set contains drugs from each class: KIs, anthracyclines, cardiac-acting, non-cardiac-acting.
2. **Non-cardiac-acting priority:** The non-cardiac-acting class is the key test for non-cardiotoxic drug prediction. Hold out a meaningful fraction (e.g., 20–30%) of non-cardiac-acting drugs.
3. **Scaffold diversity:** Compute Bemis-Murcko scaffolds for all drugs. Ensure held-out drugs have scaffold diversity (not just scaffold variants of training drugs). If a held-out drug shares its scaffold with a training drug, the unseen-drug test is easier (scaffold leakage).
4. **No drug leakage:** Verify that no held-out drug (by InChIKey) appears in the training set. This is a hard constraint.
5. **L1000 overlap note:** Held-out drugs may appear in L1000 (pretrain). This is NOT leakage — L1000 is pretrain-only and never enters cardiomyocyte evaluation. The pretrained drug encoder is the transfer mechanism for unseen drugs.

### 5.4 Split report

Generate a report with:
- Per-class drug counts in train vs. test
- Scaffold diversity: number of unique Bemis-Murcko scaffolds in train vs. test, and the number of test drugs whose scaffold also appears in train
- L1000 overlap: fraction of held-out drugs that have L1000 pretrain embeddings
- Sample counts: total samples in train vs. test (per dataset and overall)

## Unit tests

1. **SMILES canonicalization:** For a set of known drug name variants (e.g., "imatinib" vs "Gleevec" vs "STI-571"), verify they all canonicalize to the same InChIKey.
2. **No drug leakage:** Verify the intersection of train_drugs and test_drugs (by InChIKey) is empty. This is a hard-fail test.
3. **Class balance:** Verify the test set contains at least one drug from each class present in the full dataset. Report per-class train/test counts.
4. **Scaffold diversity:** Verify the number of test drugs whose Bemis-Murcko scaffold also appears in train is reported. Flag if >50% of test drugs share scaffolds with train (scaffold leakage risk).
5. **Non-cardiac-acting coverage:** Verify the test set contains at least 2 non-cardiac-acting drugs (the key test class). If the non-cardiac-acting pool is too small to split, document this as a limitation.
6. **Sample count sanity:** Verify the train/test sample split ratio is approximately 80/20 (or the chosen ratio). Flag if any dataset is entirely in train or entirely in test.

## Integration test

Validate the split is evaluation-ready:
- The split is deterministic (same input → same split) — verify by running twice and checking identical output.
- No drug leakage (hard constraint).
- Every drug in train and test has a canonical SMILES (from step_04).
- The split report contains all required fields (per-class counts, scaffold diversity, L1000 overlap, sample counts).
- The non-cardiac-acting class has representation in the test set (or a documented limitation if the pool is too small).

## Acceptance criteria

- [ ] All drugs canonicalized to InChIKey; duplicate names merged
- [ ] Drug overlap matrix computed across all 4 cardiomyocyte datasets
- [ ] Unseen-drug split constructed with class stratification
- [ ] No drug leakage verified (hard constraint)
- [ ] Scaffold diversity reported; scaffold leakage risk flagged if >50%
- [ ] Non-cardiac-acting drugs represented in test set (or limitation documented)
- [ ] Split report generated with all required fields

## Risks & notes

- The non-cardiotoxic drug pool is small (~25 drugs). Holding out 20–30% means only ~5–8 non-cardiotoxic drugs in the test set. This limits statistical power for non-cardiotoxic drug prediction evaluation.
- Scaffold leakage is a concern: many kinase inhibitors share scaffolds (e.g., quinazoline, pyrimidine). If held-out KIs share scaffolds with training KIs, the unseen-drug test is easier than truly novel scaffolds.
- Drugs unique to one dataset (e.g., Rosell-Hidalgo's 9 non-cardiotoxic compounds) are natural held-out candidates, but holding out an entire dataset's unique drugs means that dataset contributes only training samples for shared drugs, reducing its training value.
- The split must be deterministic and saved to a file so that all downstream steps (training, evaluation) use the same split.
