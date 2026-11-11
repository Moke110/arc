# Stage 2 significance retrieval

Retrieval is evaluated **only on tenv**. Models trained on Tahoe are still
scored when their target observations are tenv. No Tahoe test PR/AUPRC is
computed. The active reporting metrics are `AUPRC`, `AUPRC_positive`, and
`AUPRC_negative`. All use measured significance on the full 18,164-gene CM
axis. Recall includes significant genes outside prediction coverage. Correct
hits require the measured DESeq2 log2 fold-change direction. Zero predictions
abstain; equal magnitudes enter together. AUPRC is exact non-interpolated
average precision within each drug–dose group, averaged equally over groups
with a defined score. A zero-positive group is undefined; a group with
positives but no predictions has score zero.

Combined retrieval ranks nonzero predictions by absolute Δ and counts
significance with the same direction. Positive retrieval ranks only positive
predictions against every significant positive CM gene; negative retrieval
ranks only negative predictions against every significant negative CM gene.
Combined AUPRC is not the mean of the two directional AUPRC values.

Each PR curve pools hits and predictions at common magnitude thresholds over
a fixed eligible group set. Stop at attained recall; do not extrapolate or
remove groups at high recall. Display curves sample thresholds, while scalar
AP is computed at every distinct score. The PR curve area is not macro AUPRC.

`build_labels.py` builds checkpointed evaluator-only DESeq2 releases from raw
counts and matched controls. It keeps the CM-plus-HDGFL3 supplemental test
family used by Step 09 (the supplemental gene is excluded from retrieval).
Predictions and trained checkpoints remain unchanged. Historical scoring and
figures are removed from the working tree; Git retains the previous releases.


`inventory.json` binds the existing selected predictions, ordered case tables,
gene axes, and deterministic reference inputs. Verify these bindings with
`python stage2/evaluation/inventory.py`. Scoring rejects a mismatched case
order, gene axis or significance-contrast membership. Run scoring under the
repository's capped resource policy; at most two heavy jobs may run repo-wide.

```bash
python stage2/evaluation/evaluate.py score --step 9
python stage2/evaluation/evaluate.py replay --step 9
```

The other supported step numbers are 3, 4, 5 and 6. `--dataset tenv` is the
only supported retrieval dataset. Frozen Tahoe continuous metrics remain
separate from retrieval reporting. Complete
releases use `output/labels/<dataset>/`; a checked partial label release can
serve a cohort only when every group's exact membership is present. Missing
groups fail closed. Per-group tables retain undefined reasons and counts.

Outputs are generated in each step's `output/retrieval/`: three-metric
statistics, dose summaries, threshold-sampled pooled PR points, PNG/SVG figures
and a hash manifest. Replay recomputes every statistic and both figures from
frozen predictions and labels, then compares their hashes. Learning curves,
model weights, reconstruction scores and non-retrieval figures keep their
existing workflows. Previous rank-based metric commands are retired.
