# Step 13 — Figures & Report

## Goal

Generate publication-quality figures and a comprehensive PDF report summarizing the data inventory, model architecture, training results, and evaluation findings.

## Dependencies

- **Step 12 (Evaluation):** Provides the evaluation results (metrics, baselines, per-class/per-regime/per-drug breakdowns).
- **Step 11 (CM Training):** Provides training logs (loss curves, MMD curves, per-platform loss).
- **Step 5 (Drug Split):** Provides the drug overlap matrix and split report.
- **Step 1 (Data Inventory):** Provides the dataset inventory table.

## Inputs

- Evaluation report from step_12 (all metrics, baselines, per-class/per-regime/per-drug results)
- Training logs from step_11 (loss curves, MMD loss, per-platform loss)
- Drug overlap matrix and split report from step_05
- Dataset inventory from step_01

## Outputs

- A set of publication-quality figures (PNG format)
- A comprehensive PDF report summarizing the full project

## Behavior

### 13.1 Figures

Generate the following figures (data plots using seaborn + matplotlib, saved as PNG):

1. **Data inventory overview:** Bar chart of sample counts per dataset, grouped by cell type and platform. Drug/chemical count annotations per dataset. GSE262419's 4,558 samples will dominate the chart — use a log scale or split into two panels (small datasets: DToxS, Rosell-Hidalgo, GSE247541; large dataset: GSE262419) to keep all datasets visible.
2. **Drug overlap heatmap:** Heatmap showing which drugs/chemicals appear in which datasets. Drugs/chemicals on y-axis, datasets on x-axis, filled cells indicating presence. With ~500+ chemicals across 5 datasets, the full heatmap may be too dense — filter to show only chemicals present in ≥2 datasets (shared drugs) plus a summary count of dataset-unique chemicals.
3. **Unseen-drug split visualization:** Scatter plot of held-out vs. training drugs/chemicals in chemical space (e.g., t-SNE of Morgan fingerprints), colored by train/test split and shaped by drug/chemical class. The chemical space now includes environmental chemicals from GSE262419 — the t-SNE will show broader coverage and potentially distinct clusters for drug-like vs. environmental chemical space.
4. **Training loss curves:** Line plot of training and validation loss over epochs, with MMD loss and reconstruction loss on separate panels. Per-platform loss as additional lines (5 platforms).
5. **MMD convergence:** Line plot of MMD loss over epochs, showing cross-platform latent alignment improving during training.
6. **Evaluation: model vs. baselines:** Grouped bar chart comparing the model and all 4 baselines on primary Δ-metrics (Pearson_delta, MSE_delta, frac_correct_direction_20). Faceted by drug/chemical class (6 classes: KIs, anthracyclines, cardiac-acting, non-cardiac-acting, antiviral, environmental chemical).
7. **Evaluation: per-regime comparison:** Bar chart of model performance (Pearson_delta) across the four generalization regimes (unseen/seen class, unseen/unseen class, cross-platform, cross-platform TempO-Seq→RNA-seq).
8. **Per-drug prediction quality:** Heatmap or dot plot of per-drug/chemical Pearson_delta for all held-out drugs/chemicals, annotated with drug/chemical class and Tanimoto similarity to nearest training drug/chemical. With ~100+ held-out chemicals, this figure may need to show only the top/bottom performers or use a compact dot plot format.
9. **Predicted vs. true Δ scatter:** For a selection of held-out drugs, scatter plot of predicted Δ vs. true Δ across genes, with Pearson correlation annotated.

### 13.2 PDF report

Generate a comprehensive PDF report (using the pdf-report-generation skill) with the following sections:

1. **Executive summary:** Project goal, approach, key findings, honest assessment of model vs. baselines.
2. **Data inventory:** Table of all 5 datasets, sample counts, drug/chemical counts, platforms. Non-cardiotoxic drug pool description (validated drugs + GSE262419 environmental chemicals). TempO-Seq platform description. Concentration-response design note for GSE262419. GSE244109 sub-series relationship note.
3. **Architecture:** Model diagram (M1–M5 composition), evidence basis for each component, key design decisions.
4. **Training:** Phase 1 (L1000 pretrain) and Phase 2 (cardiomyocyte training) summaries, loss curves, MMD convergence.
5. **Evaluation:** Primary results table (model vs. baselines on all metrics), per-class breakdown, per-regime breakdown, per-drug analysis.
6. **Honest assessment:** Does the model beat baselines? Which components show promise? What are the limitations?
7. **Limitations & future work:** Small sample size for non-TempO-Seq datasets, non-cardiotoxic drug scarcity (mitigated by GSE262419 environmental chemicals but with caveats), platform heterogeneity challenges (5 platforms including TempO-Seq), concentration-response sample inflation, environmental chemical SMILES coverage, potential improvements.

### 13.3 Figure quality checks

For each figure:
- Verify axes are labeled with units.
- Verify color palettes are colorblind-friendly.
- Verify text is readable (font size, no overlapping labels).
- Verify figure is saved as PNG at publication resolution (≥300 DPI).
- Run media output check on each figure.

## Unit tests

1. **Figure generation:** Verify each of the 9 figures is generated as a PNG file with non-zero file size. Verify no figure is blank or corrupted.
2. **Figure content:** For each figure, verify it contains the expected data (e.g., the training loss curve figure has lines for training and validation loss, with 5 platform lines). Verify axes are labeled.
3. **Report generation:** Verify the PDF report is generated with non-zero file size. Verify it contains all 7 sections.
4. **Report completeness:** Verify the report references all figures. Verify the evaluation section contains the model-vs-baseline comparison table with 6 drug/chemical classes and 4 regimes.
5. **Media output check:** Run visual QA on each figure — verify no blank, unreadable, clipped, or low-quality figures.
6. **Data inventory accuracy:** Verify the data inventory figure and report section reflect 5 datasets, 5 platforms, and the GSE262419/GSE247541 additions. Verify GSE244109 is noted as a sub-series, not a separate dataset.

## Integration test

Validate the figures and report are complete and accurate:
- All 9 figures are generated, non-corrupted, and pass media output check.
- The PDF report contains all 7 sections and references all figures.
- The evaluation results in the report match the evaluation results from step_12 (no data fabrication).
- The honest assessment in the report is consistent with the baseline comparison (if the model doesn't beat baselines, the report says so).
- The data inventory in the report matches step_01 (5 datasets, 5 platforms, GSE244109 sub-series noted).
- The architecture description in the report matches step_00.

## Acceptance criteria

- [ ] 9 figures generated as PNG, all pass media output check
- [ ] Data inventory figure reflects 5 datasets (with log scale or split panels for GSE262419 dominance)
- [ ] Drug overlap heatmap shows 5 datasets (filtered for readability)
- [ ] PDF report generated with all 7 sections
- [ ] Report references all figures
- [ ] Evaluation results in report match step_12 output (6 classes, 4 regimes)
- [ ] Honest assessment consistent with baseline comparison
- [ ] Data inventory and architecture descriptions match earlier steps

## Risks & notes

- If the model does not beat baselines (a real possibility per Ahlmann-Eltze and PerturBench), the report must state this honestly. The report's value is in the honest assessment and identification of which components show promise, not in claiming success.
- The per-drug prediction quality figure may reveal that the model works well for some held-out drugs (structurally similar to training drugs) but poorly for others (novel scaffolds). This is expected and should be reported.
- Figure generation should use the Phylo color palette for consistency: #000000, #ECE9E2, #FAF9F3, #E9ED4C, #FF9400, #75A025, #FD9BED, #0279EE.
- The PDF report is the primary deliverable for non-technical users. It should be self-contained — a reader should understand the project without reading the step files.
