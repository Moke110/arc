# Stage 2 map #261 — coverage and drug-specific increment

November 21, 2026 · Step 03 ADR0078 paired G_5000 views. No fitting, checkpoint selection or outer-test tuning was performed. The audit reaggregates accepted outer-test Δ PCC and uses frozen validation rows to describe chemical support.

## Observations

Cold-drug-dose raw Δ PCC for latent dose embedding was 0.2431 higher in Tahoe than tenv under profile weighting. The Mean_cell gap was 0.2408. Latent dose embedding improved over the paired Mean_cell baseline by +0.0304 in Tahoe and +0.0281 in tenv. The resulting cross-domain increment difference was **+0.0022** (95% parent bootstrap interval [−0.0204, +0.0235]), or **+0.0034** [−0.0192, +0.0242] with parent weighting.

Cold-drug held all doses for each test parent. The raw profile-weighted latent embedding gap was +0.2658, and the Mean_cell-adjusted gap was +0.0339 [+0.0164, +0.0499]. This difference is descriptive across domains. For both scenarios the no-dose method's gain over Mean_cell was smaller in Tahoe than in tenv. Full domain estimates and paired contrasts are in [published_scores.csv](published_scores.csv) and [published_domain_gap.csv](published_domain_gap.csv).

Tahoe had more fit profiles but fewer fit parents in every fold. Cold-drug fit partition sizes were 39,024–39,826 profiles and 210–211 parents across 50 cell lines in Tahoe, compared with 2,511–2,529 profiles and 275 parents in one hiPSC-CM context in tenv. The median was one treated well per Tahoe parent × dose × cell line and three per tenv condition. Cold-drug-dose fit partitions contained 24,124–27,797 Tahoe profiles and 1,674–1,686 tenv profiles with the same parent and cell counts. Each frozen test profile was counted once. Cold-drug test covered 65,318 profiles / 352 parents in Tahoe and 4,206 / 459 in tenv. Cold-drug-dose covered 21,243 / 352 and 1,402 / 459, respectively. See [coverage.csv](coverage.csv), [condition_coverage.csv](condition_coverage.csv) and [replicate_counts.csv](replicate_counts.csv).

The actual dose grids were 0.05/0.5/5 µM in Tahoe and 0.2/1/10 µM in tenv. Validation folds share neither an exact measured dose nor a cell line. Chemical support overlaps in some coarse ECFP4 similarity bands. For example, the cold-drug validation fold-2 Tahoe >0.75–1 band is empty. The occupied common bands are recorded by fold, and within-band scores remain tied to each domain's cells and doses. Counts use parent identities, and nearest-fit similarity includes exact fingerprint matches where present. Training rows have their own fingerprint in the reference set. See [development_similarity.csv](development_similarity.csv) and [development_joint_support.csv](development_joint_support.csv).

A matched estimate conditional on dose and cell context is unidentified with these splits. The baseline-adjusted differences isolate an increment over this train-only baseline algebraically. They cannot attribute the raw difference to training sample size, independent replicates, assay properties, or control reliability. Those factors vary together across the observed domains.

## Replay and estimands

The audit verifies the paired input manifest and used metadata/fingerprint digests, the published pooled CSV, and all used run score files against the tracked baseline acceptance receipt. It reconciles accepted outer-test profile means within absolute tolerance 1e−12 for all three methods before estimating intervals. It checks original profile/source-row identity and exact score alignment on profile, parent, dose and cell context. Expression and checkpoint payloads are unnecessary.

Profile weighting averages per-profile Δ PCC. Parent weighting first averages the same per-profile scores within parent, then gives each parent equal weight. All contrasts use the intersection of finite score rows and preserve pairing across methods. [published_parent.csv](published_parent.csv) records the scores and contrasts per parent. Per-condition well counts report source profiles within a cell context. They do not count independent biological replicates. For validation strata only the two existing TranSiGen scores are used because the accepted Mean_cell run has no validation scores.

The percentile bootstrap draws 4,000 independent resamples of parent clusters within each source/scenario (seed 261), retaining all profiles of each selected parent. Cross-domain intervals subtract independently resampled domain summaries. These intervals describe conditional parent-sampling variability with the released fits held fixed, including no fit seed, split or cross-domain pair variance. Split seed 226 and fit seed 42 are frozen.

Run from the repository root with the `arc` environment. Substitute `--step` with the materialized Step 03 path when running from an isolated checkout.

```sh
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 /usr/bin/time -v python \
  stage2/03_tsg_dose/execution/axis_local/coverage_audit.py \
  --step stage2/03_tsg_dose --output /tmp/stage2-261-replay
python -m unittest discover -s stage2/03_tsg_dose/execution/tests -p test_coverage_audit.py
```

The real-input run took 4.80 s at 249,284 KiB peak RSS. [run.json](run.json) pins the generated CSV digests and seeds.
