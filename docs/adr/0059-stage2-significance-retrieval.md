# Use three full-CM significance retrieval metrics throughout Stage 2

**Status:** accepted (2026-11-11). Amends the retrieval/reporting portions of
ADR 0030 §3 (the scoped output placement below), and ADRs 0042, 0046–0049, 0057 and 0058 following the user's metric consolidation.

Evaluate retrieval only on **tenv**, as clarified by the user. Tahoe-trained
models predicting tenv remain in scope; the Tahoe test set has no PR/AUPRC
report and requires no new significance build.

Use only `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` for Stage 2 retrieval
statistics. All three use measured significant genes on the same 18,164-gene
CM axis, including genes outside a model's coverage in the recall denominator.
Match the measured DESeq2 log2 fold-change direction. Retain the Step 09
DESeq2 two-sided Wald and within-contrast BH convention over CM plus HDGFL3,
then exclude the supplemental gene from retrieval. Use matched same-plate
controls, additionally matched by cell line for Tahoe. Missing measured CM
genes remain untestable slots rather than invented significant genes.

Combined retrieval ranks every nonzero prediction by absolute Δ. Positive
retrieval ranks only positive predictions and uses every significant positive
CM gene in its recall denominator. Negative retrieval does the corresponding
operation for negative predictions and significant negative genes. Equal
magnitudes enter together. Zero predictions abstain. A group with zero
relevant significant genes is undefined; a group with positives but no
predictions scores zero. Combined AUPRC is not the mean of directional AUPRC.

Compute exact non-interpolated AP within each drug–dose group, then macro
average defined groups within run. Summaries give runs equal weight and report
defined counts. Display each PR curve by pooling hits and predictions at
common magnitude thresholds over a fixed group set. Curves stop at their
attainable recall; no extrapolation or recall-conditioned group removal.
Bounded threshold sampling affects display only, never scalar AP. The pooled
PR area is a different statistic from macro AUPRC.

`stage2/evaluation/` owns the single scoring and rendering implementation.
Steps 03–06 and 09 consume it; setup, corpus and reconstruction steps do not
acquire artificial perturbation retrieval scores. Preserve existing continuous
metrics and frozen prediction/checkpoint identities. This is retrospective
rescoring, not a new checkpoint-selection decision.

Remove superseded scores and PR/metric figures, including historical copies,
from the working tree as explicitly requested. Git history is their retention
mechanism. Keep input provenance, significance labels, trained checkpoints,
loss histories and non-retrieval figures. Missing significance labels must be
reported as pending, never converted to zero or replaced by rank-defined truth.

## See

- [Earlier direction-matching correction](0058-step09-direction-matched-pr.md)
- [Shared evaluator](../../stage2/evaluation/README.md)

For these canonical retrieval releases only, amend ADR 0030 §3: Steps 03,
04, 05 and 09 keep the current generated tables, figures and manifest together
in `output/retrieval/`. Step 06 keeps PNG/SVG figures directly in `output/`
and numeric tables plus manifest in `output/retrieval/`. The shared evaluator's
label matrices remain ignored large assets under `stage2/evaluation/output/labels/`.
This exception does not change input provenance or other step layouts.

Removed admission and prior-exposure metadata may be read from the pinned Git
revision through `evaluation/history.py`, without restoring old score files.
Step 03's original completed run verification has a narrowly authenticated
compatibility path: committed inventory/run hashes, exact non-source identity,
and the reviewed current five-file runtime must all match the committed
compatibility record. Any further predictor code change invalidates this path.
