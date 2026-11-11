# HVG and shared indexed evaluation

These interfaces belong to 01_bases and require no deleted Stage 2 consumer.
Add this directory (execution) to the Python import path.

## HVG release

build_hvg.py defaults to output/training-v2/tenv.h5ad and seals that exact
file hash in its release. It reads only the selected training h5ad's small metadata and X
in 64-row windows. It never discovers or opens significance files. Its recipe
recovers linear CP10K with expm1 from the published float32 natural-log1p
expression, computes sample variance/mean, and standardizes log dispersion in
20 equal-width bins of log1p mean. Singleton bins score one; undefined
dispersion sorts last. Ties use the fixed gene index. Rankings have no
perturbation direction.

    python stage2/01_bases/execution/build_hvg.py
    python stage2/01_bases/execution/build_hvg.py --verify

The default new bases-cold-parent-v1 cohort does not reuse historical
experiments. SHA256(seed|desalted parent) orders parents, then three equal
contiguous buckets define test folds. The first ceil(10%) of remaining
parents is validation; the rest is training. All doses and repeats of a
parent remain in one role. The exact sample_id, drug_parent, dose_um,
fold, role rows, seed, algorithm and population hashes are published.
The --cohorts option instead accepts an explicit complete partition; no
parent-dose contrast may cross roles.

The descriptive full-tenv population is for inspection only. Model selection
must declare the exact training population:

    from bases.hvg import HVGRelease
    release = HVGRelease("stage2/01_bases/output/hvg-tenv-v1")
    selected_gene_indices = release.select(
        "fold_0", 2000, training_sample_ids=train_sample_ids
    )

A descriptive population or a different training population raises an error.
verify_hvg seals the training file and independently reruns all published
ranks and population membership. Fixture tests mutate held-out expression
while preserving the training rank and reject split contrasts.

## Numeric metrics and measured truth

bases.metrics.retrieval(prediction, signed_truth, covered_gene_indices)
retains the original exact AP definition. Core arrays use integer gene
indices. Predictions rank by absolute group delta for combined retrieval,
and by positive/negative magnitude for directional retrieval. Wrong signs
do not count as hits, equal magnitudes enter together, zeros abstain, and
recall denominators include significant truth outside prediction coverage.
A tested group with no positives is undefined; an untestable contrast is a
separate state.

score_groups averages predicted well deltas within each full desalted
parent-dose group before scoring. It accepts the evaluator-only
PartitionTruth interface, which verifies complete group membership for the
requested role. There are no prediction significance tests or filters.
Returned summaries are equal-group macro AP; summarize_runs then uses
equal run weights. continuous(truth, prediction) exposes PCC, Spearman,
RMSE and R2 on the declared numeric axis, preserving constant-axis undefined
correlations and R2.

The standalone evaluator reads prediction .npy arrays through a memory map:

    python stage2/01_bases/execution/evaluate_predictions.py \
      --predictions predictions.npy --samples prediction_samples.csv \
      --genes prediction_genes.csv \
      --truth stage2/01_bases/output/truth-tenv-v1 \
      --assignments stage2/01_bases/output/hvg-tenv-v1/cohorts.csv \
      --fold 0 --role test --label Model --output evaluation-output

The samples CSV gives sample_id, drug_parent, dose_um in prediction row
order; the genes CSV gives integer gene_index in prediction column order.
The assignment file defines complete groups. Predictions must cover every sample
in the selected fold and role exactly once; omitted groups, repeated samples
and samples from another role are rejected before scoring. The evaluator publishes
per-group AP, macro summaries, input seals and the plotted numeric source.
It does not train a model or restore historical predictions.

The long-running truth release seals its original frozen training-v1 input.
Using its labels with training-v2 cohorts requires an accepted equivalence
envelope proving the same raw counts, exact treated/control identities,
parent-dose membership and normalization. The original truth manifest is not
rewritten to claim a different input hash. This implementation does not
claim that equivalence until the v2 source acceptance is complete.

## Plotting and replay

bases.plots.render_retrieval writes a tenv-only three-panel PNG, exact
plotted CSV and a sealed JSON summary. The drawn curve pools a fixed
eligible group population at shared magnitude thresholds; its displayed AP
is the exact group macro, not area under the sampled pooled curve. No
Tahoe significance or PR publication is supported.

    python stage2/01_bases/execution/replay_pr.py evaluation-output/pr.json replay.png

The contract example rendered during implementation is synthetic hand-worked
evidence, not a model result. Its CSV replay reproduced identical PNG bytes.
Actual model reports require real predictions and an accepted truth release.

## Actual-data diagnostic acceptance

After training-v2, the complete truth release and its accepted successor
envelope exist, the following commands replay actual HVG populations and
export a diagnostic panel:

    python stage2/01_bases/execution/build_hvg.py
    python stage2/01_bases/execution/build_hvg.py --verify
    python stage2/01_bases/execution/verify_hvg_metrics.py       --training stage2/01_bases/output/training-v2       --hvg stage2/01_bases/output/hvg-tenv-v1       --successor stage2/01_bases/output/truth-tenv-v2       --output stage2/01_bases/output/hvg-metrics-diagnostic-v1

Run these only in an approved dataset slot, with the documented memory and
thread caps. The isolation verifier uses temporary HDF5 virtual datasets:
original training rows are referenced unchanged, and every held-out row is
replaced by a legal single-gene CP10K spike. It recomputes each fold and
requires bit-exact training statistics and ranks, changed descriptive
statistics, and rejection of descriptive or wrong-population selection.
It never copies or materializes the whole expression matrix.

There are no surviving model predictions in this scope. The diagnostic
therefore ranks the measured delta layer against measured significance
truth through SuccessorTruth. All three test-role populations must cover
all 4206 wells and 1377 parent-dose groups exactly once. Each fold rejects
truth access for a training-role group. Exported per-group CSV, fixed-group
PR source, PNG and replay are marked diagnostic; the title explicitly says
not model performance. These numbers are mechanical acceptance evidence,
not a trained model result. Numeric success remains
NUMERIC_PASS_VISUAL_PENDING until the rendered PNG is inspected and that
separate visual QA is recorded.

The generic evaluator also supports an accepted successor:

    python stage2/01_bases/execution/evaluate_predictions.py       --predictions diagnostic.npy --samples samples.csv --genes genes.csv       --successor stage2/01_bases/output/truth-tenv-v2       --training stage2/01_bases/output/training-v2       --assignments stage2/01_bases/output/hvg-tenv-v1/cohorts.csv       --fold 0 --role test --label "Measured delta" --diagnostic       --output diagnostic-output

The --diagnostic flag sets both report metadata and the visible title.
Ordinary model prediction evaluation omits that flag; neither mode trains
a model or alters any significance release.
