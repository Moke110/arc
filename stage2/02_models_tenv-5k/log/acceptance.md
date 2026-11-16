# #227 shared-input and Mean slice acceptance

2026-11-16, arc conda Python. Bases accepted view read from the main checkout;
no source bytes modified. Candidate inputs `/tmp/arc-tenv5k-inputs` and Means
`/tmp/arc-tenv5k-means` are retained for downstream execution. Frozen compact
identity/config/digest records are in `frozen-inputs/`; all payloads reconstruct
with `execution/prepare.py` against that accepted view. Formal stable publication
remains #231's responsibility.

Freeze verified 1,402 1 µM treated wells, 459 normalized parents, 5,000 genes and
192 distinct original controls. All ten fold partitions completed; cold parents
are disjoint, warm training covers validation and test parents. Mean slice has
25 complete outputs and five explicitly unavailable cold Mean_drug entries.
Each complete output retains ordered identities, expression predictions and both
well-level and parent-aggregated Δ metrics. Mean_cell is identical to Mean.

Five public-boundary tests pass, including real Bases fixture reader through
Mean/scoring, original-control deduplication, test-label perturbation isolation,
scenario separation, hand-worked parent aggregation, constant/nonfinite handling
and reordered prediction rejection. Identity rejection was red (unexpected keyword)
then green; `identity-red.txt` retains that evidence. Earlier implementer recorded
Mean red before its implementation. The accepted view's true delta equals X minus
its matched control; no normalization is refitted.

Both real commands ran singly under 10 GiB virtual-memory cap and time-v, with
coordinator monitoring every 10 seconds; preflight RAM available was 41 GiB.
Peak RSS: preparation 391,612 KiB (3.36 seconds), Means 325,240 KiB
(14.33 seconds). Verified all 459 parent-to-Morgan mappings, exact frozen
Δ = X − control, and Mean_cell = Mean for all ten folds. Full resource receipts: `prepare.time`, `means.time`. No GPU used. No final model
or selection claim is made at this intermediate acceptance boundary.
