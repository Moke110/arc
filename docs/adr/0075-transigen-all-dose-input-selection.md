# TranSiGen all-dose input selection

Status: accepted protocol (2026-11-17). The user confirmed shared understanding
and all 147 fits, then authorized specification/ticket publication only;
implementation and training must wait for a later execution session. Source:
user grill-with-docs session and
[decision memory](https://app.notion.com/p/41bab6347e9681e0b03bf98b05a4a2fe).

Published specification: [#232](https://github.com/Moke110/arc/issues/232).
Implementation/experiment tickets: #233 (no-dose vertical slice), #234 (four
100-dimensional concatenation arms), #235 (scaling/FiLM), #236 (primary63 fits),
#237 (auxiliary63 fits), #238 (final-selection21 fits), and #239 (release
acceptance). All remain future work; issue publication is the endpoint of this
session.

Create Stage 2 `03_tsg_dose` to compare dose inputs for the selected TranSiGen
architecture on the existing tenv G_5000 view. Use all 4,206 treated wells,
459 drug parents and three doses (0.2, 1, 10 µM); each dose has 1,402 wells
and all 459 parents. Retain Bases normalization, genes, controls and identities.
Keep the accepted single-dose comparison independent and Stage 3 frozen.

## Seven fixed candidates

For the four concatenation arms, let s be log10(dose/1 µM), standardized
using training-only mean and standard deviation. Every dose vector has 100
dimensions, superseding the earlier 64/32-dimensional proposal.

| Candidate | Dose interface |
|---|---|
| No dose | Original control and Morgan2048 inputs, trained on all doses |
| Fingerprint scaling | Morgan × log10(1 + dose/1 µM), without standardization |
| Input numeric | Append 100 copies of s to Morgan; encoder input2148 |
| Input embedding | Linear1→100, ReLU, Linear100→100; append to Morgan, input2148 |
| Latent numeric | Append 100 copies of s to [control latent100, chemical400]; heads input600 |
| Latent embedding | Linear1→100, ReLU, Linear100→100; append at the same heads, input600 |
| FiLM | s→Linear32→ReLU→Linear800; chemical embedding becomes (1+gamma)×embedding+beta |

Neither embedding has a final activation. Numeric repetition has no square-root
rescaling. Both mean and logvariance heads receive the concatenated vector.
FiLM's final layer is zero initialized, starting from identity modulation;
its 400-dimensional scale and shift are not concatenation dose vectors.
Preserve native initialization elsewhere. New dose Linear layers use the seeded
PyTorch default initialization; FiLM overrides its final weight and bias to zero.
Record constructor order, software version and initialization in execution identity.

Keep native loss, sampling, control encoder, latent100, decoder800,
beta0.1/dropout0.1, batch64, Adam1e-4/weight_decay1e-5, no scheduler,
float32 and reviewed CUDA policy. Modify project adapters, never immutable source.
Report parameter counts and runtime; the result compares fixed configurations,
not optimally tuned architecture families.

## Primary selection and evaluation

Cold-drug generalization is primary: isolate every parent with all its aliases,
doses and repeats. Use three complete outer folds, seeds17/42/2026, shared
partitions across seven candidates. Within outer training, use a parent-held
validation partition (ceil10% of parents); no outer-test data may select epochs,
methods or preprocessing. Keep split seed226 and freeze all identities before
training. This nested design uses a single inner holdout, not inner k-fold CV.

Each fit selects its checkpoint by complete fixed-randomness validation
treated-expression MSE, strict improvement, patience50/max1000, best restoration.
Compare restored checkpoints by validation parent-equal Delta PCC, averaging
the three seeds' scores before choosing one method for that outer fold.
Use those restored checkpoints directly for that
fold's outer testing, without a second refit on inner validation data.

For each seed, average predicted and measured Delta within parent–dose repeats,
compute each group's gene-wise metric, average doses within parent, then average
parents. Average seed metrics, not predictions, for reported comparison values.
PCC is primary; Spearman, RMSE, R², per-dose and per-well results are auxiliary.
Preserve explicit constant-vector, undefined-metric and failure accounting.
Use parent-paired bootstrap (10,000 draws, seed2262026), retaining all doses and
seeds of each sampled parent. Seeds are not independent biological observations.

The primary OOF estimate evaluates the fold-specific selection procedure.
Candidate-wise outer results are descriptive and cannot choose the final method.
Separate a highest-validation-score deployment choice from evidence of scientific
superiority; no-dose may win and insufficient evidence must remain inconclusive.
Exact validation ties prefer fewer parameters, then a preregistered simple-to-
complex ordering. Do not claim universal superiority from exploratory pairwise CIs.

## Auxiliary dose generalization

Compare all seven candidates with each dose held out completely: 1 µM tests
interpolation, 0.2 and 10 µM test lower/upper extrapolation. All evaluated parents
must occur at another training dose. On the two available doses, split validation
by complete parent–dose groups, retaining another dose in training for every
validation parent. Do not fit transforms, early stop or select methods using
the held dose. Report all methods and highlight methods selected by primary
cold-drug validation; auxiliary outcomes never change the primary winner.

## Final model and execution boundary

Final deployment selection: train all seven methods on one frozen
full-population development split with parent-held validation and three seeds;
choose by the same validation rule and publish its three best checkpoints.
This deployment choice has no additional independent test score. No training
uses test results to alter candidates or numerical settings.

The approved complete workload is 63 primary fits +63 auxiliary fits +21
deployment-selection fits =147 fits, excluding bounded smoke tests and failed
attempts. This session publishes issues only and does not start these fits.

Nonfinite loss/prediction terminates a fit and preserves evidence. No ad hoc
seed replacement, LR reduction or loss change; incomplete comparisons are not
successful releases. New candidates have independent splits, method/execution
identities, exact fresh-process sampled replay and complete prediction coverage.
Apply ADR0062 candidate-first publication and CSV-only figure reproduction.
Run long training singly with repository memory caps and resource monitoring.

This separates dose information from additional training samples, avoids averaging
away dose response, and prevents method selection from consuming outer test data.
Three measured dose levels constrain the scope of interpolation/extrapolation
claims. Repository plans and implementation details remain in the repository.

## Execution authorization amendment — 2026-11-17

In a later session the user invoked `implement-spec` with “部署stage2 03”,
authorizing implementation and the complete147-fit workload through acceptance
and deployment. This supersedes the publication-only execution boundary above;
the scientific protocol remains unchanged. Formal acceptance and resulting model
selection are recorded in the step README and acceptance log when complete.
