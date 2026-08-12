# Step 12 — Evaluation: Metrics & Baselines

## Goal

Evaluate the trained cardiomyocyte drug-perturbation model on the unseen-drug test split using Δ-metrics, rank metrics, and distributional metrics, with mandatory comparisons against simple baselines — producing an honest assessment of whether the model adds value.

## Dependencies

- **Step 11 (CM Training):** Provides the trained model checkpoint.

## Inputs

- Trained model checkpoint (M1–M5 weights)
- Unseen-drug test split (from step_05) — held-out drugs and their cardiomyocyte samples
- Harmonized cardiomyocyte test data: control expression, true Δ, drug SMILES, dose, covariates

## Outputs

- Evaluation report with per-metric, per-class, per-regime results
- Baseline comparison table (model vs. each baseline on each metric)
- Per-drug prediction quality summary
- Honest assessment: does the model beat baselines?

## Behavior

### 12.1 Evaluation metrics

**Δ-metrics (primary):**
- **Pearson_delta:** Correlation of predicted vs. true Δ across genes, averaged over samples. Measures whether the model captures the direction and magnitude of expression change.
- **MSE_delta:** Mean squared error between predicted and true Δ. Measures absolute prediction accuracy.
- **frac_correct_direction_{20,50,100,200}:** Fraction of top-k differentially expressed genes (DEGs) where the predicted direction (up/down) matches the true direction. Measures whether the model identifies the right genes and their direction of change.

**Absolute metrics (secondary, expected to be misleadingly high):**
- **Pearson_absolute:** Correlation of predicted_treated vs. true_treated absolute expression. Expected to be high because most genes don't change — this metric is dominated by the baseline expression, not the drug effect.
- **MSE_absolute:** MSE on absolute expression. Reported for completeness but not the primary metric.

**Rank metrics:**
- **rankRMSE:** Root mean squared error of ranked drug effects. Assesses whether the model correctly orders drug effects.
- **cosine-LogFC-rank:** Cosine similarity of log-fold-change rank profiles (PerturBench convention). Assesses whether the model correctly orders genes by perturbation magnitude.

**Distributional metrics:**
- **MMD in PCA space:** Maximum mean discrepancy between predicted and true Δ distributions in PCA space. Measures whether the model captures the distribution of responses, not just the mean.

### 12.2 Baselines (mandatory)

The model MUST beat these baselines to claim success. Per Ahlmann-Eltze and PerturBench, simple baselines are hard to beat:

1. **No-change baseline:** Predict Δ = 0 for all drugs. The model must beat this on Δ-metrics (Pearson_delta, MSE_delta, frac_correct_direction).
2. **Mean-drug-response baseline:** Predict Δ = mean Δ across all training drugs. The model must beat this to show it learns drug-specific effects.
3. **Linear PCA model:** A linear model on PCA-reduced expression (Ahlmann-Eltze's strong baseline). The model must beat this to justify the deep-learning complexity.
4. **Per-class mean baseline:** Predict Δ = mean Δ for the drug's class (KI, anthracycline, cardiac-acting, non-cardiac-acting, antiviral, environmental chemical). The model must beat this to show it learns within-class drug-specific effects.

### 12.3 Generalization regimes

Report performance separately for each regime:
- **Unseen drug, seen drug class:** Held-out drug whose class is represented in training (e.g., held-out KI, other KIs in training). Expected to be easier.
- **Unseen drug, unseen drug class:** Held-out drug whose class is not in training (e.g., held-out non-cardiac drug, no non-cardiac drugs in training — hardest regime).
- **Cross-platform:** Train on iPSC-CM datasets, test on primary-CM DToxS (and vice versa). Tests cross-cell-type generalization.
- **Cross-platform (TempO-Seq → RNA-seq):** Train on TempO-Seq (GSE262419) samples, test on RNA-seq platform samples (or vice versa). This is the hardest cross-platform test, as TempO-Seq is a targeted probe-based platform fundamentally different from RNA-seq. If the model succeeds here, MMD regularization effectively integrated the 5th platform. If it fails, TempO-Seq-specific effects are leaking into predictions.

### 12.4 Per-class reporting

Report all metrics separately for each drug/chemical class:
- Kinase inhibitors (KIs)
- Anthracyclines
- Cardiac-acting drugs
- Non-cardiac-acting drugs (the key test class for non-cardiotoxic drug prediction)
- Antiviral drugs (GSE247541: apilimod, remdesivir, ritonavir, lopinavir)
- Environmental chemicals (GSE262419 chemicals without significant transcriptomic effect — the expanded non-cardiotoxic pool)

### 12.5 Per-drug analysis

For each held-out drug/chemical, report:
- Predicted vs. true Δ Pearson correlation
- Top-20 DEG direction match rate
- Whether the model beats the no-change baseline for this specific drug/chemical
- Drug/chemical structural similarity to nearest training drug/chemical (Tanimoto similarity) — to contextualize prediction quality
- Chemical class annotation (drug vs. environmental chemical) — to assess whether the model generalizes differently to drug-like vs. non-drug-like structures

## Unit tests

1. **Pearson_delta correctness:** Given predicted_Δ and true_Δ, verify Pearson_delta is in [-1, 1]. Verify it equals 1 when predicted = true. Verify it equals 0 when predicted is random.
2. **frac_correct_direction correctness:** Given predicted_Δ and true_Δ, verify frac_correct_direction_k is in [0, 1]. Verify it equals 1 when all top-k directions match. Verify it equals ~0.5 for random predictions.
3. **No-change baseline:** Verify the no-change baseline (Δ=0) produces Pearson_delta = 0 (or undefined — document behavior) and frac_correct_direction = 0.5 (random direction).
4. **Mean-drug-response baseline:** Verify the mean-drug-response baseline produces a non-trivial Pearson_delta (captures the average drug effect direction).
5. **MMD computation:** Given two sets of Δ vectors, verify MMD is non-negative, zero when the sets are identical, and positive when they differ.
6. **Per-class reporting:** Verify metrics are computed and reported separately for each drug/chemical class (KIs, anthracyclines, cardiac-acting, non-cardiac-acting, antiviral, environmental chemical). Verify no class is missing from the report.
7. **Regime separation:** Verify the four generalization regimes (unseen drug/seen class, unseen drug/unseen class, cross-platform, cross-platform TempO-Seq→RNA-seq) are reported separately.

## Integration test

Validate the full evaluation pipeline:
- The trained model (from step_11) loads correctly and produces predictions for all held-out drugs/chemicals.
- All 4 baselines are computed and reported.
- All metrics (Δ-metrics, absolute, rank, distributional) are computed for the model and all baselines.
- The evaluation report contains per-class (6 classes) and per-regime (4 regimes) breakdowns.
- The honest assessment is stated: does the model beat each baseline on each primary metric?
- Per-drug analysis is generated for all held-out drugs/chemicals, including chemical class annotation.
- No held-out drug/chemical was used in training (verify against step_05 split).

## Acceptance criteria

- [ ] All primary Δ-metrics computed (Pearson_delta, MSE_delta, frac_correct_direction_{20,50,100,200})
- [ ] All 4 baselines computed and compared
- [ ] Rank metrics computed (rankRMSE, cosine-LogFC-rank)
- [ ] Distributional metric computed (MMD in PCA space)
- [ ] Per-class results reported (KIs, anthracyclines, cardiac-acting, non-cardiac-acting, antiviral, environmental chemical)
- [ ] Per-regime results reported (unseen/seen class, unseen/unseen class, cross-platform, cross-platform TempO-Seq→RNA-seq)
- [ ] Per-drug analysis generated for all held-out drugs/chemicals
- [ ] Honest assessment stated: model vs. baselines on each primary metric

## Risks & notes

- Given Ahlmann-Eltze and PerturBench evidence, there is a real risk the model does not beat linear baselines. If it doesn't, report this honestly and identify which component (if any) shows promise (e.g., the model may beat baselines on non-cardiac-acting drugs but not on KIs).
- The non-cardiac-acting test class now has robust representation (~84–126 held-out chemicals from the expanded pool). However, most of these are GSE262419 environmental chemicals, not validated non-cardiotoxic drugs. Report results separately for validated non-cardiotoxic drugs vs. environmental chemicals to assess whether the model generalizes differently to these sub-groups.
- The absolute metrics (Pearson_absolute, MSE_absolute) are expected to be misleadingly high because most genes don't change. Do not use these as primary evidence of model quality — they are reported for completeness and to demonstrate why Δ-metrics are the primary metrics.
- The cross-platform regime (train on iPSC-CM, test on primary-CM) is the hardest test of the MMD regularization. If the model fails here but succeeds on same-platform unseen-drug tests, it indicates the cross-platform integration is incomplete.
- The cross-platform TempO-Seq→RNA-seq regime is the newest and potentially hardest test. TempO-Seq's probe-based measurement differs fundamentally from RNA-seq. If the model fails specifically here, it indicates that MMD regularization did not fully integrate TempO-Seq with RNA-seq platforms — consider whether TempO-Seq needs a separate harmonization track (a design decision that was deferred in step_04).
- If the model beats the no-change baseline but not the linear PCA baseline, it means the model captures some drug-specific signal but not enough to justify deep-learning complexity. This is a partial success worth reporting.
- Environmental chemicals may have very different structures from the L1000 pretrain compounds. The pretrained M2 encoder may produce poor embeddings for these chemicals. If per-class results show environmental chemicals performing much worse than drug classes, this indicates the drug encoder doesn't generalize to non-drug chemical space.
