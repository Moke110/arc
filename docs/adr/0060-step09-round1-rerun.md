# Rerun Step 09 Round 1 with six models and validation early stopping

**Status:** accepted (2026-11-11). Amends ADR 0057's Round 1 arms, run matrix,
continuous scoring and artifact placement. ADR 0059 retrieval remains unchanged.

Replace the exploratory nine-arm, one-run screen with G512, G1024, G2048,
MLP512, MLP1024 and MLP2048 crossed with folds 0/1 and seeds 17/42 from
Step 03. G arms use nested fold-local training-only HVGs; MLP arms retain the
full 18,164-gene axis. Preserve existing normalization and loss definitions.
This compares smaller gene axes and compressed full-axis models within the
available workstation budget.

Use cold-cm-drug folds: Tahoe has a shared train/validation partition and may
contain tenv test drugs; tenv train/validation exclude that fold's test drugs.
Train each run from scratch, Tahoe then tenv, with effective batch 16 and fresh
AdamW per phase (learning rate 1e-4, weight decay 1e-5). Prefer microbatch 4
with four accumulation steps, falling back to 2 with eight. Tahoe budgets are
1–40 training-set passes; tenv budgets are 1–5. Validate the complete validation
cohort every 50 and 10 parameter updates respectively, plus phase termination.
A pass counts training-well presentations, not validation or backward calls.

Track validation from the beginning but permit stopping only after one pass.
Stop on either three consecutive strictly increasing adjacent losses (four
points), or ten consecutive points strictly above the historical minimum.
Equality resets the applicable counter; counters are phase-local. Termination
is at the first eligible validation point or the hard maximum update budget.
Select best weights from all validation points, including before the first
pass. Thus minimum training exposure does not constrain checkpoint eligibility.
Save one best checkpoint per model/fold/seed run; restore Tahoe best before
fine-tuning, then replace its file with tenv best. Keep both phase histories
and selection provenance, but only the final tenv-best weights at completion.

Use native-axis test loss, PCC, Spearman, RMSE and R². These intentionally do
not enforce a common gene panel. All three AUPRC metrics use ADR 0059 full-CM
truth, with uncovered genes absent from predictions rather than imputed hits.
Compute per-run pooled directional PR before aggregating four runs, equally
weighted; retain attainable recall and do not extend a curve beyond its endpoint.

Round 1 output contains only checkpoints/, logs/, figure_data/,
r1_learning_curves.png, r1_metrics.png and r1_PR_curves.png. This extends the
scoped ADR 0030 placement exception to checkpoints and logs. Remove obsolete
training products and arms only after preserving historical summaries and
curves outside Round 1; retain shared data and full-axis inputs needed by MLPs.
