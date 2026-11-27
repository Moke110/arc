# Stage 3 预测优化 — recorded option space for a Δ-based viability-direction model

Read-only survey assembled 2026-11-13 for the Stage 2 → Stage 3 decision map. Only this
file was written; no corpus, model, ADR, branch or commit was touched. Every quantitative
claim is tagged:

- **[verified]** — re-read or recomputed from a tracked artifact during this survey.
- **[lead]** — supplied by the Lead, not re-verified here (reason given where relevant).
- **[inference]** — reasoned from code or artifacts, not itself a measurement.
- **[not measured]** — no artifact exists in this repo.

This is a survey of what the repository already fixes, measures and has decided. It is not
a new ADR and does not select a model. Scope split with the companion decision map
(`docs/research/stage2-stage3-prediction-optimization.md`): that note owns the **Stage 2
Δ-shape** question (conditioning, effective rank, metric fix); this note owns the **Stage 3**
consumption question (Δ processing, model class, evaluation, leakage) and cites it rather
than re-deriving it.

---

## Fixed at the boundary

These are frozen by accepted decisions or by the completed formal benchmark. Re-opening any
of them invalidates comparability with the frozen Morgan floor.

| # | Fixed item | What is frozen | Source |
|---|---|---|---|
| 1 | Source and assay | Heliste 2021 Dataset S1, HL-1 mouse atrial cardiomyocytes, normalized CTG ATP-luminescence `cell_viability`; 24 h pretreatment + 48 h hypoxia + 24 h reoxygenation, duration-matched normoxia arm, ~96 h elapsed, 0.01% DMSO controls. Persistence and the exact normalization denominator are unresolved. | `docs/adr/0044-stage3-source-first-corpora.md`; `cpsp/01_corpus/README.md`; `docs/research/stage3-heliste-normalization-label-audit.md` |
| 2 | Corpus shape | 689 source-label rows, 3 unresolved identities, 686 resolved labels, 663 unique desalted parents, 4 doses (0.04/0.2/1/5 µM), 2,756 cells per condition, 2,744 admitted per condition, 92 extra rows from 22 multi-alias compounds (21 with 8 rows, 1 with 12). One record = one supplied numerical cell; `rep: null`. | `cpsp/01_corpus/README.md`; `cpsp/01_corpus/output/{baseline,H-R}/corpus.csv` **[verified]** |
| 3 | Label rule | `direction ∈ {−1,0,+1}`: <80 → −1, 80–120 inclusive → 0, >120 → +1. Magnitude-only, under the unconfirmed assumption that each arm's own 0.01% DMSO vehicle = 100. Class 0 is within-band, **not** demonstrated equivalence; no significance test. | ADR 0044 "Magnitude direction labels amendment"; corpus README |
| 4 | Model interface | One model per `(condition, phenotype)`; three-class softmax over [−1,0,+1]; no shared multi-phenotype output; target is magnitude direction, not therapeutic benefit. | `docs/adr/0045-stage3-morgan-direction-baseline.md` |
| 5 | Splits and seeds | Split seed `20261105`; cold-drug 3 outer folds (test parents never in training, doses may overlap); cold-dose 4 global held-dose folds (0.04 and 5 µM = extrapolation, 0.2 and 1 µM = interpolation); warm-start 3 folds; seeds 17/42/2026; all aliases and doses of a parent share one role across conditions; 60 fits for 2 conditions. | ADR 0045; `cpsp/02_base_morgan/execution/split_manifest.json`; `cpsp/02_base_morgan/output/formal/summary.json` |
| 6 | Baseline envelope and schedule | Morgan radius 2 / 2048 binary bits / chirality + `log10(dose_um)` standardized on training-partition statistics → `Linear(input,128) → ReLU → Linear(128,3)`; unweighted cross-entropy; AdamW; lr 1e-3, wd 1e-4, batch 64, ≤300 epochs, patience 30, restore best validation macro-F1. | ADR 0045 |
| 7 | Metric convention | Macro-F1 primary; per-class precision/recall, 3×3 confusion, log loss; each held dose reported before its average; cold-drug additionally reports pooled out-of-fold; `class_order [−1,0,1]`, argmax ties first in class order, `zero_division 0`, log-loss clip to `[eps, 1−eps]` then row-normalize. | `summary.json` `metric_convention`; `cpsp/02_base_morgan/execution/metrics.py` |
| 8 | Frozen floor numbers | H-R cold-dose macro-F1 0.5991 / 0.5881 / 0.5710 / 0.4910 at 0.04 / 0.2 / 1 / 5 µM; H-R cold-drug pooled 0.4123 ± 0.0230; H-R warm-start 0.5301 ± 0.0149. Baseline cold-dose 0.5687 / 0.5782 / 0.5693 / 0.5377; cold-drug pooled 0.4081 ± 0.0044; warm-start 0.5504 ± 0.0448. Class +1 F1 = 0 in every cold-dose and cold-drug aggregate. | `cpsp/02_base_morgan/output/formal/report.md`, `summary.json` **[verified]** |
| 9 | Delivery boundary | Stage 3 Δ work is a later design step; the comparison to the Morgan baseline must be paired on jointly covered samples if Δ coverage is incomplete, so that model comparison does not depend on different evaluation populations. | ADR 0045 "Delivery boundary" |

### What is genuinely open

1. **The Δ release contract (Lead's Q1–Q5, unconfirmed by the user).** Not written down in
   any tracked doc; the five items are: the 663×4 grid, the 4-run ensemble, the corpus
   control mean as baseline, the 0.04 µM dose extrapolation, and the fixed gene axis with
   zero-variance genes dropped.
2. **The processing ladder** (`centre → PCA-32 → amplitude → unit-L2 → dose`) — recommended,
   not accepted.
3. **Model class** — explicitly undecided in the user's brief.
4. **Primary generalisation setting for the decision.** The frozen report headlines cold-dose,
   but the deployment case is a compound the model has not seen; cold-drug is the honest
   setting and is not the frozen headline.
5. **The model interface.** ADR 0045 says "Future Δ predictors concatenate Δ and dose and use
   the same hidden and output architecture, changing the first-layer input dimension." A
   18,164-dim concat is 2.33 M parameters; whether the "input" may be a processed Δ
   representation is a decision to record, not to assume.
6. **Gene-set provenance** for pathway/WCS features (no accepted Stage 3 gene sets exist).
7. **Acceptance bar and metric additions** for a class-+1-starved problem.
8. **Whether baseline or H-R is the headline target** for the Stage 3 optimisation claim.
9. **Assay wording.** The user's example is "CCK8/MTT"; the only frozen corpus is CTG
   (`cell_viability`, ADR 0044). The label is assay-class-agnostic but the training data are
   not CCK8/MTT.

### Evidence status: three Lead-supplied facts that do not reproduce

Verified from the frozen CSV, and therefore used everywhere below instead of the Lead's
versions:

| Lead-supplied | Recomputed **[verified]** |
|---|---|
| change rate rises 4.2% → 11.8% (baseline), 3.8% → 11.4% (H-R) with dose | baseline 4.4 / 5.8 / 9.0 / 17.8 %, H-R 3.9 / 4.7 / 7.0 / 14.3 % of rows per dose |
| 140/663 (baseline) and 110/663 (H-R) compounds have any non-zero label; 523/553 are all-zero | 146/663 (baseline) and 112/663 (H-R) compounds have any non-zero label; 517/551 are all-zero |
| 46 compounds labelled −1 only at 0.04 µM | only 4 baseline / 3 H-R compounds are −1 only at 0.04 µM; 29 / 26 rows are −1 at 0.04 and all but 4 / 3 of those compounds are also −1 at a higher dose |

The real asymmetry is different and more decision-relevant: −1 outnumbers +1 at every dose
(baseline 194 vs 60 rows; at 5 µM 81 vs 41), singleton labels cluster at the **top** dose
(baseline 73 of 95 singleton labels are at 5 µM, versus 4 at 0.04 µM), and the asymmetry is
structural — +1 requires >120% of vehicle while −1 only requires <80%, and cytotoxicity is
easier to produce than a 20% viability gain. Do not build the analysis around low-dose
label noise; build it around class +1 being rare and structurally hard.

---

## Non-learned Δ processing: candidate list

### The measured shape of the predicted Δ (what any processing must respect)

The only quantitative comparison is the fold-0 probe on 1,401 test wells against the
measured tenv Δ (`tmp/stage3_delta_probe/probe.py`, `pca.npz`, `summary.npz`).

| Quantity | Predicted Δ̂ | Measured Δ | Source |
|---|---|---|---|
| Participation ratio (effective rank) | **4.00** | **14.85** | `pca.npz` **[verified]** |
| Exp-entropy effective rank | 7.53 | 202.1 | `pca.npz` **[verified]** |
| PC1 variance share | 45.48% | 24.91% | `pca.npz` **[verified]** |
| PC2 / PC10 / PC16 / PC32 share | 17.80 / 90.52 / 98.19 / 99.91% | 6.07 / 38.45 / 40.63 / 44.16% | `pca.npz` **[verified]** |
| Median per-gene std across wells | 0.0655 | 0.1845 (ratio 0.355) | `summary.npz` **[verified]** |
| Median squared cosine of a row with its own mean vector | 0.425 | 0.101 | `summary.npz` **[verified]** |
| ‖global mean Δ‖ | 12.07 | 4.26 | `summary.npz` **[verified]** |
| cos(predicted mean, measured mean) | 0.350 | — | `summary.npz` **[verified]** |
| Cross-dose direction correlation | 0.945 (Lead, all-well aggregation) / **0.979** (group-level, median over 218 within-drug dose pairs) | 0.313 (Lead) / **0.138** (group-level) | **[lead]** for 0.945/0.313; the group-level pair is from the companion note `docs/research/stage2-stage3-prediction-optimization.md` §9, which re-checked it from `prediction.npy` + bounded `layers/delta` slabs. Aggregation level explains the difference; both estimators agree that dose barely changes the predicted direction. `probe.py` prints the all-well version but did not persist it into `summary.npz`, and this survey did not load the truth layer |
| Within-group noise as fraction of between-group | 0.5% | — | **[lead]** |

Code-level cause **[inference]**: `stage2/02_models_tenv/execution/script/prnet_fit_predict.py:427`
computes `chemical = adapter(Morgan) * log10(dose + 1)` and feeds it to the native PGM. The
prediction is therefore approximately a rank-1-in-dose chemical signature scaled by
`log10(dose+1)`: direction is nearly dose-invariant, and amplitude is nearly a deterministic
function of (compound, dose). This is consistent with the recorded predicted cross-dose
correlation of ≈0.95–0.98 and with PC1 = 45% of variance, and it means the delivered Δ cannot
express the measured change in *which* genes respond across doses (recorded measured
correlation 0.14–0.31).

Consequence for the whole decision: **the predicted Δ has at most ~4 effective dimensions,
~3 after removing the shared mean direction.** The Stage 3 question is whether those
directions carry phenotype-direction information, not whether a model has enough capacity.
Any processing that produces more than a handful of dimensions is fitting noise amplitude,
and a 32-dimensional PCA is already 28 dimensions more than the predicted data supports.

### Candidate table

Definitions used: a **transformation** is a fixed per-row mapping of the Δ matrix (no external
information, no labels); a **feature** is a new scalar/vector channel added to the model input.
"Label-free" means no Stage-3 `direction` is used to define the mapping; choosing `k`, `λ` or
a cutoff **using validation labels inside a fold** is leakage-safe but must be declared.

| # | Candidate | What it does | Why it might help | Cost | Label-free | Type | What it destroys | Redundancy |
|---|---|---|---|---|---|---|---|---|
| 1 | Centring / shared-component removal | subtract the corpus mean Δ vector, or project out the top 1–2 PCs of the prediction matrix | predicted rows share 42.5% squared cosine with their own mean (measured 10.1%); the shared axis otherwise dominates any linear/MLP fit and is collinear with amplitude | trivial | yes | transformation (linear projection) | the shared axis — which may itself carry the generic viability signal (Szalai's global viability axis warning) | centring is the rank-1 case of PC removal; "remove top-2 PCs then PCA-k" ≡ PCA-(k+2) with the first two dropped |
| 2 | PCA / truncated SVD | fit an orthogonal basis on the predicted Δ population, keep `k` scores | label-free compression; 99.91% of predicted variance in 32 components, so it is lossless for the predicted corpus; decorrelates the model input | trivial | yes | transformation | only the discarded 0.09% of predicted variance; with `k<4` it destroys real predicted structure | subsumes #1 if centred; whitening subsumes most of #7; unwhitened PCA keeps amplitude in the scores |
| 3 | Pathway / gene-set activity | fixed gene-set × Δ projection (linear activity) or ssGSEA/rank activity | denoises gene-level instability; standard escape hatch in the literature (GRIP-Lung ssGSEA; PROGENy/DoRothEA) | needs a gene-set collection mapped to the frozen 18,164 axis; low compute | yes | feature | within-pathway gene detail; resolution collapses to set membership | **linear activity ≡ a linear projection** and adds no inductive bias beyond dimension reduction; only rank-type (ssGSEA/KS) scores add bias |
| 4 | WCS-style projection | repo's Weighted Connectivity Score: GSEA/weighted-KS enrichment of an up-set and a down-set on the Δ ranking, `WCS=(ES_up−ES_down)/2` when signs oppose, else 0 | signed, rank-based, robust to the distorted dispersion of predicted Δ; the repo already has an audited implementation | trivial once gene sets exist; needs H-R/MIRI up/down sets whose provenance is undecided | yes **if** the gene sets are not derived from Stage-3 labels | feature (2–3 scalars) | everything except enrichment relative to the chosen sets | a nonlinear special case of #3; WCS with weight 1 is a KS statistic, not a linear projection |
| 5 | Amplitude channel | `log‖Δ̂_centred‖` (one scalar per row) | change probability rises strongly with dose/effect size; the Lead records r ≈ 0.43 with measured effect size **[lead, not re-verified]** | trivial | yes | feature | nothing; it is a summary of the row norm | redundant if the main representation is amplitude-preserving and unnormalized; **essential** if #6 is applied |
| 6 | Unit-L2 / cosine normalisation | divide each row by its L2 norm | makes rows comparable; separates pattern from size | trivial | yes | transformation | all amplitude information (which is phenotype-relevant) | complementary to #5; after #6, #5 is the only amplitude left |
| 7 | Per-gene standardisation | z-score each gene across the Δ population | stops high-variance genes dominating; predicted per-gene std is fairly uniform (IQR 0.048–0.087 around a 0.0655 median), so the effect is mild | trivial | yes | transformation | nothing material for these data; changes the PCA metric from covariance to correlation | absorbed by PCA-with-whitening; otherwise a diagonal variant of #2 |
| 8 | Rank / quantile transform | replace a row's values by their ranks or normal scores | monotone-robust to outliers and to the predicted dispersion distortion | trivial | yes | transformation | magnitude and linearity; incompatible with a linear model | useful only as the input to #3/#4, where it is already implicit |
| 9 | Sign×\|Δ\| sparsification / top-k | keep the k largest \|Δ\| genes, zero the rest | reduces to a workable dimension without a learned basis; matches the field's "rank by effect size, do not p-value" rule | trivial | the fixed-k map is label-free; **choosing k must be declared** | transformation | the small-`abs(Δ)` tail — which may be the only non-shared signal | a hard-threshold version of #10; different only in that it zeroes rather than scales |
| 10 | Shrinkage | global λ·Δ, or per-gene empirical-Bayes shrinkage using the PRnet predicted variance head | the literature's recommended substitute for significance filtering on predicted data; the variance head exists in the checkpoint | low | yes (λ may be tuned on validation) | transformation | magnitude if λ is small; per-gene shrinkage can erase the shared axis | shrink-to-mean = #1; global λ is redundant after #6; only the per-gene variance-weighted version is distinct |
| 11 | Dose-relative Δ | predict Δ at a fixed reference dose (e.g. 1 µM) and give only actual dose to the head; or supply both plus the difference | decouples Stage-2 dose behaviour at unmeasured doses (0.04 µM is below the `log10(dose+1)` training range; 5 µM is inside it but not a measured training dose) from Stage-3 dose generalisation | one extra inference pass | yes | feature / alternative input | the genuine dose dependence of the transcriptome (measured cross-dose corr 0.14–0.31) | **nearly degenerate given the conditioning**: Δ̂ is approximately `f(compound)·log10(dose+1)`, so ref-dose Δ and actual-dose Δ are collinear up to scale. The informative split is pattern-only vs pattern+amplitude (#5/#6), not ref vs actual |
| 12 | Ensemble-dispersion / uncertainty channel *(extra, not in the Lead's list)* | variance across the 4 accepted runs (or across noise draws) per compound–dose, as a scalar or per-gene weight | the 4 runs differ by fold and seed, and PRnet's decoder takes a random noise vector at inference (`prnet_fit_predict.py`, `predict_deltas`), so dispersion is a label-free reliability signal | requires repeated inference | yes | feature | nothing | orthogonal to every other candidate; cheap and currently absent |
| 13 | Signed log / variance-stabilising transform | `sign(Δ)·log1p(abs(Δ))` | none specific; the target is already on a log1p CP10K Δ scale | trivial | yes | transformation | magnitude ordering | subsumed by #7; not recommended |

### Redundancy and the minimal non-redundant set

The candidate list collapses to five orthogonal choices:

1. **Remove or keep the shared mean/leading axis** (#1, and the leading components of #2).
2. **Choose the Δ basis**: PCA-k (#2), per-gene scaling (#7), top-k (#9), or shrinkage (#10) —
   all are linear or diagonal re-encodings of the same matrix; do not stack them as if they
   were independent evidence.
3. **Split amplitude from pattern** (#5, #6). These are two halves of one decision: if you
   normalise, you must add amplitude back; if you do not, amplitude is implicit.
4. **Add a nonlinear gene-set statistic, or not** (#3 linear = redundant; #4/#3-rank are the
   only genuinely new inductive biases).
5. **Choose the dose parametrisation** (#11), which is nearly degenerate under the current
   Stage-2 conditioning.

#12 adds a reliability signal that no single Δ contains; #3/#4 add external gene-set
knowledge. Everything else is a re-encoding of the same Δ matrix.

### What the shape evidence supports

- **Supported:** centring (or dropping PC1–PC2), PCA with k ∈ {4, 8, 16, 32}, an explicit
  amplitude scalar, and a per-gene/ensemble reliability weighting.
- **Weak/unsupported as a primary path:** unit-L2 alone (discards the dose-scale signal that
  drives the label), pathway-40 alone (a rank-40 linear map of a rank-4 Δ spans at most four
  directions, so it adds no information over those four), full-gene input without reduction
  (18,164 features for 1,648 cold-dose training rows, i.e. ~11 features per row; see
  arithmetic), and Δ@reference-dose treated as a *separate* hypothesis from amplitude.
- **Genuinely open:** whether removing the shared axis helps or hurts, because the shared
  axis may be the generic viability/proliferation signal that the literature says is the
  dominant phenotype axis (`tmp/lit-review/transcriptome-to-phenotype-landscape.md` §2.1,
  Szalai's global viability factor and its 0.63 → 0.74 AUC removal result).

---

## Model selection: options and arithmetic

### Sample-size ledger [verified from corpus CSV and `summary.json`]

| Population | Baseline | H-R |
|---|---|---|
| Rows / unique structures / parent–dose combos | 2,744 / 663 / 2,652 | 2,744 / 663 / 2,652 |
| Label-bearing rows (−1/+1) | 254 (194/60) | 205 (178/27) |
| Compounds with any non-zero label | 146 | 112 |
| Compounds ever −1 / ever +1 / always 0 | 93 / 54 / 517 | 87 / 26 / 551 |
| Cold-drug training fold rows, label-bearing rows (3 folds) | 1,468/1,444/1,476; 148/129/137 | 1,468/1,444/1,476; 118/102/113 |
| Cold-dose training rows, label-bearing rows (4 folds) | 1,648/1,643/1,647/1,645; 175/169/149/101 | 1,648/1,643/1,647/1,645; 140/130/126/83 |
| Cold-dose test rows per held dose | 686 | 686 |
| Cold-dose 0.04 test class counts (−1/0/+1) | 29/656/1 | 26/659/1 |
| Cold-dose 0.2 test class counts | 34/646/6 | 32/654/0 |
| Cold-dose 1 test class counts | 50/624/12 | 42/638/6 |
| Cold-dose 5 test class counts | 81/564/41 | 78/588/20 |

Floors computed analytically from these counts:

- all-zero predictor macro-F1 = **0.3172** (baseline) / **0.3204** (H-R); accuracy 90.7% / 92.5%.
- Three-class macro-F1 with F1(+1) forced to 0 has ceiling **0.6667**. The frozen H-R
  cold-dose scores at 0.04/0.2 µM (0.588–0.599) are therefore 88–90% of that ceiling, the
  5 µM score (0.491) is 74%, and H-R cold-drug (0.412) is 62%.
- A dose-only argmax model is the all-zero predictor at every dose (0 is the majority at all
  four doses), so **dose is a strong predictor of change probability and of log loss, but not
  of direction**. Any direction gain must come from the Δ interaction.
- The frozen baseline is already overparameterized: 262,787 parameters for 1,468 cold-drug
  training rows (179 parameters/row), and 1,776 parameters per label-bearing training row.

### Parameter ledger

`Linear(d,128) → ReLU → Linear(128,3)` has `128d + 515` parameters; multinomial logistic has
`3(d+1)`. Input dimensions follow the candidates above.

| Input | d | Logistic params | MLP-128 params | MLP params per label-bearing cold-drug row (148) |
|---|---|---|---|---|
| dose only | 1 | 6 | 643 | 4 |
| PCA-8 + amplitude + dose | 10 | 33 | 1,795 | 12 |
| PCA-32 + amplitude + dose | 34 | 105 | 4,867 | 33 |
| pathway-40 + dose | 41 | 126 | 5,763 | 39 |
| Morgan-2048 + dose (frozen floor) | 2,049 | 6,150 | 262,787 | 1,776 |
| Δ (18,164) + dose, literal ADR 0045 concat | 18,165 | 54,498 | 2,325,635 | 15,714 |
| Δ + Morgan + dose, early concat | 20,213 | 60,642 | 2,587,779 | 17,486 |

Reading: the ≤41-dim **linear** path has fewer parameters than the number of label-bearing
training rows; the matching MLP has ~33 parameters per label-bearing row, borderline but
regularisable. The literal ADR 0045 concat is ~15,700 parameters per label-bearing row and
~848 parameters per training row overall — not a model-selection question but a variance
question.
The 18,164-dim *linear* model (54k parameters) is defensible **only** with strong L1/L2
regularisation and a declared validation-selected penalty.

### Options

**(i) Regularised multinomial logistic regression on PCA features.** `d ≈ 34`, 105 parameters.
- Risk: underfits interactions; linear in dose so cold-dose extrapolates smoothly; the PCA
  basis is fixed and label-free; no seed variance (one fit per fold).
- Protocol: same folds; L2 penalty chosen on each fold's validation partition by macro-F1,
  declared grid; report the selected penalty. Same metric module.
- Verdict: **the statistical default.** Ordinal (proportional-odds/cumulative-link) regression
  on the same features is the natural variant: the label really is a two-cutpoint
  discretisation of a continuous ratio, so the ordering carries information, and the rare +1
  class borrows from the −1 side of the latent score. Risk: the proportional-odds assumption
  on a band with asymmetric cutpoints (80/120) may be violated; fit as a second arm, compare
  by validation log loss and macro-F1.

**(ii) Gradient boosting (XGBoost/LightGBM).**
- Risk: with only four dose values, a tree cannot extrapolate `log10(dose)` to a held-out dose
  — cold-dose generalisation is structurally weak for trees; small n plus class imbalance
  needs class weights and shallow depth; temperature of the validation macro-F1 is high.
- Protocol: depth ≤ 3, high `min_child_weight`, class weights, early stopping on the inner
  validation; report as a **diagnostic** for interactions, not as the primary model. Do not
  tune on test; report feature importance only as description.

**(iii) Shallow MLP (ADR 0045 shape).**
- Mandated for the cross-stage comparison: same hidden layer and head as the frozen baseline,
  same training schedule, same early stopping. On 34 inputs it is 4,867 parameters (33 per
  label-bearing row) and is defensible.
- Risk: seed variance on a target with 60 (+1) rows in the whole baseline corpus (27 in H-R);
  early stopping on validation macro-F1 with a 0/1/2-count +1 class can select degenerate
  all-zero-+1 checkpoints.
- Protocol: run all 30 fits per condition for a real comparison; report seed SD and paired
  per-fold deltas; do not compare a 1-seed Δ MLP to a 3-seed frozen baseline.
- The literal 18,165-dim concat arm should be run **once**, declared in advance as a
  negative control, because ADR 0045's wording names it; its parameter count must be reported
  next to its score.

**(iv) kNN / similarity-to-annotated-compound.**
- Two versions: Tanimoto kNN in Morgan space (chemistry-only prior) and cosine kNN in
  Δ-PCA space (does Δ position predict direction?). No training.
- Risk: cold-drug test compounds may have no close neighbour; kNN scores will be dominated by
  the 78–83% of compounds that are always 0; and any distance in Δ space is dominated by
  amplitude/PC1 unless normalised — which reintroduces candidate #6.
- Protocol: report as a floor/diagnostic; choose k on validation only; never let a test
  compound's label enter its own neighbour set.

**(v) Hybrid Δ + Morgan + dose.**
- Published evidence favours late (probability-level) fusion over early concatenation
  (Moshkov et al. 2023, `tmp/lit-review/transcriptome-to-phenotype-landscape.md` §2.2/table
  row 7). Early concat is 20,213 dims and 2.59 M MLP parameters — not defensible here.
- Defensible form: two small encoders (PCA-34 logistic/MLP; Morgan-dims logistic/MLP) and an
  average of logits, or a 2-parameter stacker fit on inner-validation predictions only.
- Risk: double use of validation; a stacker fit on the same validation that early-stops the
  encoders is optimistic. Use nested (inner) CV or a fixed 50/50 logit average.
- Expected gain is small: the Δ side has ~4 effective dimensions, and measured L1000 added
  only +2 well-predicted assays over structure alone in Moshkov's 270-endpoint study
  (`tmp/lit-review/transcriptome-to-phenotype-landscape.md` §3 Q2).

**(vi) dose-only and Morgan-only floors.**
- dose-only: equals the all-zero predictor for argmax (macro-F1 0.317/0.320); useful only to
  show that a Δ model is not merely reading dose.
- Morgan-only: not measured anywhere in this repo — the frozen baseline is Morgan **+** dose.
  If the Δ model is compared only to Morgan+dose, the dose contribution is not isolated.
  Running Morgan-only with the same folds/schedule is cheap and should be added.
- The official floor remains the frozen Morgan+dose benchmark; the two above are ablations,
  not replacements.

### Fair-comparison rules

1. Reuse the split manifests and metric module; do not re-split. `cpsp/02_base_morgan/execution/split_manifest.json`, `metrics.py`.
2. Compare on the **jointly covered population** (ADR 0045 delivery boundary). If the Δ
   release covers 663×4, that is the full 2,744 rows; if it covers only the published
   artifacts, it is only 72 rows (35 Stage-3 compounds at the two dose-matched levels; 36
   alias-inclusive rows per dose) — too small for any comparison. The coverage scenario must be fixed before any Δ run is scored.
3. Paired, per-row comparison on identical rows; report per-fold deltas and a paired bootstrap
   CI, not two independent means.
4. Same inner-validation procedure and early-stopping criterion for the MLP arm; logistic and
   ordinal arms declare their penalty grid and select it on validation.
5. Seeds: the MLP arm uses 17/42/2026; deterministic arms produce one fit per fold and must
   say so rather than reporting three identical "seeds".
6. Compute budget: the 34-dim arms are seconds per fit; the 18k arms are still CPU-feasible
   but must follow `docs/agents/execution-environment.md` (10 GiB, 2 threads, one heavy job).

### Conclusion of the arithmetic

n is not 2,744 independent samples: it is **663 compounds / 2,652 parent–dose combinations**,
of which **146 (baseline) or 112 (H-R) compounds carry any non-zero label**, and only
**54 / 26 compounds ever carry a +1**. At that effective n, only the low-dimensional,
strongly-regularised path is statistically defensible. Full-Δ MLP is not a candidate; full-Δ
linear is a declared stress test, not a default; PCA-k + amplitude + dose with a linear or
ordinal head is the only option whose parameter count is smaller than the label-bearing
sample. The comparison the repo actually needs is therefore **not "which model class wins"
but "do the ~4 predicted-Δ dimensions add anything over Morgan + dose"**.

---

## Evaluation protocol

### What the current formal summary does (do not change any of it)

- 60 formal fits = 2 conditions × (3 cold-drug + 4 cold-dose + 3 warm-start) folds × 3 seeds;
  complete and frozen (`summary.json` `complete: true`, `formal_benchmark_complete: true`).
- Per-fold and per-seed macro-F1, per-class precision/recall, 3×3 confusion, log loss; cold-dose
  reported per held dose then averaged; cold-drug pooled out-of-fold per seed and across seeds.
- Seed/fold summaries are descriptive; seeds are not biological replicates (ADR 0045).
- Metric convention fixed in `metrics.py` and echoed in `summary.json`.

### What must be added so a Δ model can be judged honestly

**P1 — must exist before any Δ score is reported.**

1. **A pre-declared primary setting.** Recommend **cold-drug** (the deployment case: a new
   compound), with cold-dose as confirmation. The frozen report's cold-dose headline lets a
   compound be memorised across doses (all 663 parents appear in training), so it flatters any
   Δ model that encodes compound identity. Cold-drug also has more +1 test rows (19–21 per
   fold baseline vs 1 at the 0.04 cold-dose fold).
2. **Class-+1 reporting that is not structurally zero.** Add, per fold: one-vs-rest AUROC for
   +1 (defined even when argmax never picks +1), +1 recall at a fixed precision budget, and
   the binary "changed vs within-band" F1 (`|y|>0`). Keep the 3-class macro-F1 as the
   headline for comparability, but never decide on it alone.
3. **Per-dose breakdown inside cold-drug.** Dose drives change probability (4.4% → 17.8%), so
   a pooled cold-drug number mixes four different base rates. Report cold-drug metrics per
   dose level too.
4. **Jointly-covered-population statement.** State exactly which rows a Δ model was trained
   and scored on, and re-score the frozen Morgan baseline on the same rows (or use the full
   2,744 if the release is 663×4). Otherwise the two models are not comparable (ADR 0045).
5. **Paired uncertainty.** Per-row paired bootstrap or per-fold paired deltas with a CI;
   3-seed SD alone is not an uncertainty statement for the Δ-vs-Morgan question.

**P2 — should exist.**

6. **Label-noise statement from the band.** Verified boundary density: fraction of readouts
   within ±0.5 / ±1 / ±2 / ±5 / ±10 percentage points of the 80 or 120 cut is
   0.47 / 1.13 / 2.84 / 7.54 / 22.08% (baseline) and 0.40 / 0.95 / 1.93 / 5.54 / 18.40%
   (H-R) **[verified]**. There are no replicate SDs in the corpus, so the exact flip rate is
   **not measured**; but any assay noise of a few percent flips a few percent of labels, and
   this bounds the achievable macro-F1. Report the band and this density next to the score
   instead of treating the labels as noiseless.
7. **All-zero compounds.** 517/663 (baseline) and 551/663 (H-R) compounds are always 0. A
   model that predicts 0 for them is correct and gets most of the accuracy. Report a
   non-zero-compound subset metric as a **descriptive** secondary number; never make it the
   headline (it changes the population and can be gamed).
8. **The +1 structural limit at low dose.** 0.04 µM has exactly 1 +1 row in each condition and
   0.2 µM has 6 (baseline) / 0 (H-R). The +1 F1 at those folds is structurally 0 or 1; report
   it as undefined evidence, not as a model failure, and do not let it drive model selection.
9. **Ordinal-consistency check.** Because the label is a two-cutpoint discretisation of a
   continuous ratio, report mean predicted P(+1)−P(−1) (the latent direction score) against
   the continuous readout on the non-zero rows. This is a stronger test than argmax F1 and is
   the natural way to show a simple assay-detectable direction signal.

**P3 — nice to have.**

10. Calibration (reliability / ECE) per class, since probabilities will feed a prioritisation.
11. A measured-Δ sanity overlay on the 45-compound chemical intersection (see traps): sign
    agreement of the top-|Δ| genes, reported as a descriptive gate, never as a training signal.

### Minimum acceptance bar

A Δ-based model is worth its Stage 2 dependency only if **all** of the following hold. These
are proposals for the user to confirm; they are not accepted.

1. On identical folds/rows/seeds, mean paired Δmacro-F1 vs the frozen Morgan+dose model
   **≥ +0.02** with a paired 95% CI excluding 0 in **cold-drug**, and no single cold-drug fold
   worse by > 0.02.
2. Cold-dose: mean Δmacro-F1 **≥ 0** (non-inferior) and better on **≥ 2 of 4** held doses.
3. Hard-class evidence: **+1 one-vs-rest AUROC ≥ 0.55** with CI excluding 0.5 in warm-start
   and at the 5 µM cold-dose fold, or +1 F1 > 0.10 at 5 µM cold-dose.
4. Beats the all-zero/dose-only floor by ≥ +0.10 macro-F1 in every reported setting (the
   frozen baseline already clears this, so this only catches broken Δ models).
5. Does not lose −1 F1 by more than 0.02 vs the frozen baseline.
6. Reproducible: seed SD of the paired Δ ≤ 0.02, and the sign of the paired Δ agrees in both
   conditions.
7. **Adds value beyond amplitude+dose.** A pattern-based Δ model must beat an
   amplitude+dose-only model by ≥ +0.02 macro-F1. Because Δ̂ direction is nearly dose-invariant
   and amplitude is nearly dose-determined, a Δ model that ties amplitude+dose has demonstrated
   no directional transcriptomic information.
8. Never accepted on warm-start alone, and never accepted on a plus-one-only or
   log-loss-only improvement.

If only bar 4 and 5 pass, the honest conclusion is "dose + a generic effect-size axis explains
the label; the Stage 2 Δ adds nothing measurable", which is a legitimate and useful outcome.

---

## Leakage and protocol traps

**T1 — PCA fitted on the wrong population.** Fitting PCA on train folds only is the clean
choice; fitting on all compound–dose predictions without labels is label-free but
transductive and must be declared. Never fit on Stage-3 test rows *together with* their labels
(e.g. selecting components by how well they separate classes). A 34-dim input with 146
label-bearing compounds leaves no room for a "components chosen by label" step.

**T2 — Gene/component selection using Stage-3 labels.** Selecting top-|Δ| genes, pathway sets,
PCA components or the amplitude cut using correlations with the Stage-3 `direction`, even
"just to see", is leakage and will not survive an independent split. All of #1–#13 in the
candidate table are label-free by construction except the *choice of k/λ/cutoff*, which must
live inside the fold.

**T3 — The 45-compound measured-Δ gate is a gate, not a training set.** Chemical intersection
between tenv and the frozen Stage-3 corpus is **45 compounds** (InChIKey) **[verified]**; 10 of
them appear only in tenv fold 2, which has no runs, so only **35** have any published
prediction. Of those 35: baseline 144 rows / 20 non-zero rows / 11 compounds with any non-zero
label; H-R 144 / 16 / 8. At the only dose-matched levels (0.2 and 1 µM) it is **72 rows, 8
non-zero rows, 6 compounds (baseline)**. Using this intersection to train, to select features
or to choose a model is overfitting 35 compounds. It can only be a descriptive sanity gate,
and it must use the fold-appropriate tenv truth (a compound in a tenv test fold may not use
that fold's truth outside evaluation).

**T4 — Ensemble exposure is heterogeneous by construction.** The four accepted runs are
2 cold-drug folds × 2 seeds. The folds' test drug sets are disjoint (154 and 153 drugs), so a
compound in fold 0's test set was in fold 1's **training** set. Every predicted compound has
been seen in training by at least one ensemble member and held out by another. This is not
Stage-3 label leakage (the target is a different phenotype), but the ensemble is not a clean
held-out predictor for any compound, and the mixing must be recorded per compound. Compounds
outside tenv (618 of 663) are pure chemical extrapolation.

**T5 — Using H-R labels when building the H-R model (or baseline labels for H-R).** ADR 0045
requires a separate model per (condition, phenotype) and aligns parent roles across conditions.
Cross-condition transfer (train on baseline, fine-tune on H-R) would put the same compound's
baseline label in the training pool while its H-R label is in the test pool for cold-drug and
warm-start folds — legal per target but a protocol change that must be decided explicitly, not
slipped in. Multi-task heads are ruled out by the "no shared output model" decision unless
re-opened.

**T6 — Row-level metrics are alias-inflated.** 22 compounds contribute 2–3 alias rows
(92 extra rows). Row-level macro-F1 and confusion matrices therefore weight those compounds
2–3×. Report compound-level counts alongside row-level metrics; do not switch the metric
population silently.

**T7 — The 0 class is not equivalence.** The ±20% band labels real 10–19% effects as 0. A model
trained on these labels learns "within band", not "no effect"; do not interpret a 0 prediction
as evidence of no biological effect, and do not build a "detect no-change" claim on it.

**T8 — Dose standardisation leakage.** `log10(dose_um)` must be standardised with
training-partition statistics only, exactly as the frozen workflow does; in cold-dose this
means statistics from the three remaining doses. Refitting the scaler on all doses including
the held dose leaks the held dose's position.

**T9 — Early stopping on a tiny +1 class.** Validation macro-F1 with 1–6 +1 rows can select a
degenerate model. Keep the frozen criterion for fairness, but also report fixed-epoch
checkpoints and the +1 one-vs-rest AUROC of the selected checkpoint.

**T10 — Reading the truth layer during feature work.** The measured tenv truth
(`layers/delta`) is evaluator-only per ADR 0063/0064. Any per-gene or per-component selection
that uses it is a cross-stage violation even though it is not a Stage-3 label. This survey
deliberately did not load the h5ad.

**T11 — Δ@reference-dose vs Δ@actual-dose is not an independent hypothesis.** Since the
Stage 2 conditioning is approximately `f(compound)·log10(dose+1)`, the two arms are collinear
up to scale. Presenting them as two tests of "does dose-specific Δ matter" would overstate the
evidence; the real test is amplitude vs pattern (candidates #5/#6).

**T12 — Coverage-driven model comparison.** There are three coverage scenarios and they are
not interchangeable: (a) published artifacts only = 35 Stage-3 compounds at the two tenv-matched
doses = 72 rows (36 alias-inclusive rows per dose); (b) new inference on those 35 compounds at
all four Stage-3 doses = 144 rows;
(c) the 663×4 contract = 2,652 rows, which requires inferring 628 compounds that never appear
in tenv. If the Δ model trains on (a) or (b) while the frozen Morgan model trains on 2,744
rows, a lower Δ score is uninformative. Either deliver (c), or re-train/re-score the Morgan
baseline on the same rows, or restrict both to those rows. ADR 0045 already requires this;
make it an executable check, not a sentence.

---

## Dependency on Stage 2 decisions

The Lead's labels "A2" and "Q1–Q5" are **not defined anywhere in the tracked repo**
(grep found no tracked use). This survey therefore treats **Q1–Q5 = the five
unconfirmed Δ-release-contract items** listed in the brief, and **A2 = "the Stage 2
metric fix"** as defined in the companion decision map
(`docs/research/stage2-stage3-prediction-optimization.md`, "What is blocked on the
metric fix"): per-group normalised threshold `λ = d / max|Δ̂_group|` instead of one
run-wide `M`, macro average over **changed groups only**, and `@λ` naming, keeping
ADR 0064's truth side. That note also records why the current family cannot rank
models: PRnet's peak `Reliability@d*` sits at the strictest grid point with a mean
window of 0.03 genes, and changed-group `Reliability@d* = 0.000` on all four
accepted runs (`stage2/02_models_tenv/README.md`).

| Decision | Blocked by A2? | Blocked by Q1–Q5 (Δ release)? | Can be decided now? |
|---|---|---|---|
| Corpus, label rule, phenotype target, ±20% band | no | no | **yes — frozen** |
| Splits, seeds, metric convention, frozen floor numbers | no | no | **yes — frozen** |
| Processing candidate list and redundancy map | no | no (the list is label-free by construction) | **yes — this document** |
| Which processing arm is primary | no | partly (amplitude/centring need the actual Δ matrix) | **recommend now, run later** |
| Model class short-list and parameter budget | no | no | **yes** |
| Evaluation protocol additions and acceptance bar | no | no | **yes** |
| Running any Δ-based training/scoring | no | **yes** — needs the released Δ matrix | no |
| PCA-k / amplitude / centring empirical choices | no | **yes** | no |
| Measured-Δ sanity gate | no | **yes** (needs the 663×4 Δ to compare on the 45 compounds) | no |
| Δ@reference-dose arm | no | **yes** — needs extra inference | no |
| Ensemble composition and per-compound exposure record | no | **yes** (Q2) | no |
| 0.04 µM extrapolation handling | no | **yes** (Q4) | no |
| Gene axis / zero-variance filter | no | **yes** (Q5) | no |
| Whether the model input may be a processed representation (ADR 0045 wording) | no | no | **yes — needs a recorded decision** |
| Whether to amend ADR 0045 or run the literal 18k concat as a declared negative control | no | no | **yes** |
| Which Stage 2 runs/checkpoints constitute "the Δ" | **yes** if A2 changes the accepted product or the reliability gate | **yes** (Q2) | no |

The practical split: **everything in this document except the empirical Δ runs can be decided
now**; the empirical program waits on Q1–Q5, and the *identity* of the accepted Stage 2
product waits on A2.

---

## Recommended default if no further input

A concrete, reversible default that does not require the user to answer anything before work
starts. Every step is label-free until the Stage 3 head.

1. **Contract (provisional, to be confirmed):** consume Δ on 663 compounds × 4 doses from the
   ensemble of the 4 accepted runs; baseline expression = tenv corpus control mean; drop only
   zero-variance genes; record per compound (a) its tenv fold membership and (b) which
   ensemble members saw it in training; record that under the Stage 2 `log10(dose+1)`
   conditioning 0.04 µM is below the training minimum (0.017 vs 0.079) and 5 µM is inside the
   range but not a measured training dose (0.778 vs min 0.079 / max 1.041) **[verified algebra]**.
2. **Representation:** centre the predicted Δ by subtracting the label-free corpus mean, then
   PCA fit on the label-free predicted corpus, and carry `k ∈ {4,8,16,32}` as a declared
   sensitivity grid (not a tuned hyperparameter). Keep the uncentred PCA-32 arm as the
   control for candidate #1. Add amplitude `log‖Δ̂_centred‖`. Do **not** unit-normalise in the
   primary arm; keep a pattern-only arm as candidate #6. Dose enters as standardised
   `log10(dose_um)`, matching ADR 0045 — never four copies of Δ.
3. **Model:** primary = L2 multinomial logistic on PCA-k + amplitude + dose (or the
   proportional-odds ordinal variant); secondary = ADR 0045's 128-hidden MLP on the **same**
   low-dimensional input, same schedule, all 30 fits per condition. Run the literal 18,165-dim
   concat MLP once as a declared negative control and report its 2.33 M parameter count. Add
   Morgan-only and dose-only floors.
4. **Protocol:** primary setting cold-drug; cold-dose as confirmation; report per-class,
   per-dose, one-vs-rest +1 AUROC and changed-vs-not F1; paired comparisons on joint rows;
   3 seeds for the MLP, one fit per fold for deterministic heads.
5. **Acceptance:** the eight-bar checklist above. If bars 1–7 fail but 4–5 pass, record the
   negative result ("Δ adds nothing over dose + amplitude") and stop; that is the expected
   outcome given the literature and the predicted Δ's rank.
6. **Ablation order (cheapest first):** dose-only → Morgan-only → Morgan+dose (frozen) →
   amplitude+dose → Δ_PCA-k+amplitude+dose → +centring off → +pattern-only → +Morgan (late
   fusion) → pathway/WCS feature → full-Δ L1 logistic. Stop at the first arm that meets the
   bar; do not run the full ladder speculatively.
7. **Sanity gate:** on the 35 predicted Stage-3 compounds that intersect tenv, compare predicted
   vs measured group-mean Δ descriptively (top-|Δ| sign agreement, low-dim correlation). It is
   a gate, never a training input. State the 45-compound / 72-dose-matched-row limit explicitly.

---

## Citations

**Frozen Stage 3 decisions and artifacts**

- `docs/adr/0044-stage3-source-first-corpora.md` — corpus selection, CTG scope, ±20% direction rule and the vehicle=100 assumption.
- `docs/adr/0045-stage3-morgan-direction-baseline.md` — three-class interface, splits/seed, training schedule, delivery boundary and the "exactly one model per (condition, phenotype)" rule.
- `docs/adr/0043-stage3-source-measurement-corpus.md` — measurement-preserving corpus and `rep` semantics.
- `docs/adr/0039-condition-corpus-and-model-view-separation.md` — condition/phenotype partition and alias handling context.
- [`cpsp/README.md`](../../cpsp/README.md), [`cpsp/01_corpus/README.md`](../../cpsp/01_corpus/README.md), [`cpsp/02_base_morgan/README.md`](../../cpsp/02_base_morgan/README.md) — current public contracts.
- [`cpsp/01_corpus/output/baseline/corpus.csv`](../../cpsp/01_corpus/output/baseline/corpus.csv), [`cpsp/01_corpus/output/H-R/corpus.csv`](../../cpsp/01_corpus/output/H-R/corpus.csv) — the 2,744-row corpora; all class, dose and boundary statistics recomputed here.
- `cpsp/02_base_morgan/output/formal/summary.json`, `report.md`, `benchmark.json` — frozen 60-fit benchmark, per-class metrics, split role counts. [Release record](../../cpsp/02_base_morgan/execution/release_record.json)
- [`cpsp/02_base_morgan/execution/metrics.py`](../../cpsp/02_base_morgan/execution/metrics.py) — fixed metric convention.
- [`cpsp/02_base_morgan/execution/split_manifest.json`](../../cpsp/02_base_morgan/execution/split_manifest.json) — persisted partitions.
- `docs/research/stage3-heliste-normalization-label-audit.md` — why the vehicle=100 denominator is an assumption, and the missing replicate/variance evidence.

**Stage 2 artifacts that define the Δ**

- `stage2/02_models_tenv/execution/script/prnet_fit_predict.py` — the `adapter(Morgan) × log10(dose+1)` conditioning (line 427) and the `mean − baseline` prediction semantics (`predict_deltas`).
- `stage2/02_models_tenv/README.md` — the accepted 2-fold × 2-seed matrix, the 1,377-group / 918-test-group population, and the current `@d` reliability numbers.
- `stage2/02_models_tenv/output/inputs/fold_0/seed_17/test_cases.csv`, `.../fold_1/seed_17/test_cases.csv` — the 307 predicted drugs and their doses.
- `stage2/01_bases/output/tenv-truth/groups.csv` — 459 tenv drugs × 3 doses (0.2/1/10 µM); used for the 45-compound InChIKey intersection.
- `stage2/01_bases/README.md` — tenv = single-cell-line hiPSC-CM TempO-Seq at 24 h; the domain gap to HL-1/CTG.
- `docs/adr/0063-tenv-comparison-direct-and-indirect-metrics.md` — Δ = well expression − matched same-plate DMSO control; evaluator-only truth.
- `docs/adr/0064-tenv-reliability-threshold-scan.md` — the Stage 2 metric family that A2 replaces.

**Predicted-Δ shape evidence**

- `docs/research/stage2-stage3-prediction-optimization.md` — companion decision map: the Stage 2 Δ-shape options (A1/A3/A4/A5), the A2 metric-fix definition, and the group-level re-check of effective rank and cross-dose correlation.
- `tmp/stage3_delta_probe/probe.py`, `tmp/stage3_delta_probe/pca.npz`, `tmp/stage3_delta_probe/summary.npz` — 1,401 test-well probe; participation ratios, PCA shares, per-gene std, mean-vector statistics recomputed here.

**Literature (do not redo)**

- `tmp/lit-review/predicted-transcriptome-as-input.md` — no paper significance-filters predicted Δ before a downstream model; downstream rules use full vector, latent, |Δ| ranking or shrinkage; pseudo-replication caveats.
- `tmp/lit-review/transcriptome-to-phenotype-landscape.md` — 26-work landscape; Q1–Q3 answers, Szalai's global viability axis and the 0.63→0.74 AUC gene-removal result, McFarland's β₀/β₁ decomposition, Moshkov's +2-assay ceiling for measured L1000 over structure, late-over-early fusion.
- `tmp/lit-review/xixi-zhang-2025.md` — the regX verification; it is not a precedent for the Stage 3 chain.

**WCS implementation**

- `stage1/xpert-fda-miri-v3/_lib/scoring.py` — the repo's Weighted Connectivity Score (GSEA/weighted-KS enrichment of up/down gene sets on the Δ ranking).
- `GLOSSARY.md` — "MIRI efficacy score" definition and the WCS+RGES/Zhang/pathway scoring family.

**Operating constraints**

- `docs/agents/execution-environment.md` — 10 GiB `RLIMIT_AS`, two threads, one heavy job at a time.
