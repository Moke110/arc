# Use a landmark reference for the Step 09 XPert ablation

**Status:** accepted (2026-11-10).

Allow the 978-gene `XPert_base` and shared-panel scoring within the internal
Step 09 cold-cm-drug ablation. This scoped exception to
[ADR 0024](0024-stage2-bringup-978-quarantine.md) anchors the gene-space and
compression comparisons while retaining the 15,000-gene minimum for the
Stage 2 delivery model.

Project Tahoe and tenv native counts onto each model's gene axis before
calculating CP10K and log1p. The complete author landmark panel and every
reduced HVG panel use their own per-well library sizes. The MLP arms and
`XPert_G18164` share the full CM axis and its normalization. Freeze one randomly generated identity
vector per gene on the existing 18,164-gene CM axis. Give HDGFL3, the landmark
outside that axis, a supplemental vector generated under the same fixed
seed and stored separately.

Reuse the three tenv folds and training seeds 17, 42 and 2026 from
[ADR 0042](0042-tenv-model-benchmark-redesign.md). Use one Tahoe plate-based
training and validation split across folds. Initialize each model from
scratch, select a pretraining checkpoint on Tahoe validation, then select
the final checkpoint during tenv fine-tuning on tenv validation. Retain the
treated-expression, baseline-reconstruction and separate Δ heads. Train the
Δ head with MSE and Pearson losses. The first round uses fold 0 and seed 17.

Use the intersection of all planned arms' prediction axes for the primary
figures. Compare the same genes, tenv drug–dose groups, expression scale
and direction labels across arms. Report native-axis scores and gene
coverage separately.

Build nested HVG panels per fold using only Tahoe training wells and that
fold's tenv training wells. Rank normalized dispersions separately by
source and combine the ranks with equal source weight. Compare each MLP arm
with the uncompressed `XPert_G18164` on the same input and output genes.

For shared-panel scoring, add the predicted Δ to the measured matched
control on the model's own axis to recover the predicted treated profile.
Project and normalize this prediction on the shared panel, then subtract
the measured shared-panel control. Derive observed treated and control
profiles on the shared panel directly from native counts. Freeze promotion
decisions using tenv validation before evaluating the exploratory round-one
tenv test results.

The round-one figures live together under Step 09's `output/round1/` as a
scoped placement exception to [ADR 0030](0030-external-data-provenance-and-step-layout.md).
