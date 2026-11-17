# Single-dose tenv 5k native comparison

Accepted 2026-11-16 under [#226](https://github.com/Moke110/arc/issues/226).
The successor `stage2/02_models_tenv-5k` compares randomly initialized native
PRnet and TranSiGen with Mean, Mean_cell, Mean_drug and the project MLP_Morgan
on the released G_5000, 1 µM tenv population (1,402 wells, 459 parents).
This isolates architecture from pretraining and dose adaptation. All learners
share Morgan2048, cold-drug/warm-start five-fold identities and seeds 17/42/2026;
complete fixed-randomness validation treated-expression MSE selects the strictly
lowest checkpoint with patience 50 and at most 1,000 epochs.

This supersedes ADR 0063/0064's current model-comparison and evaluation duties:
primary scoring now aggregates repeated-well Δ within parent before computing
PCC, with parent-paired bootstrap and explicit constant/failure accounting.
Historical evidence retains its original scientific context. The previous step
remains accepted until the successor's complete 90-run matrix, reload, metric,
CSV-only figure and downstream acceptance passes; only then may #231 retire it.
ADR 0062's candidate publication rules still apply. Stage 3 remains frozen.

## Approved numerical-method amendment — 2026-11-17

Following reproducible non-finite loss in the frozen native float32 TranSiGen
identity, the user approved `transigen-native-f64-v1` and continuation of #226
([decision record](https://app.notion.com/p/41bab6347e96818fafa7e11ffe71dec7)).
This is an explicitly labeled numerical variant, not a successful reproduction
of the original float32 trajectory. Preserve the original failed identity and
its evidence. Initialize with the native float32 selective initialization, then
convert all floating parameters, buffers, inputs, forward/loss/backward,
Adam moments, validation and inference to float64. Integer bookkeeping stays
integer. Draw latent noise with the original float32 CPU
`FloatTensor(size).normal_()` call and draw order, then cast before scale/add.
Retain the reviewed author CUDA policy (cuDNN disabled, benchmark false,
deterministic true), native architecture, KL algebra, LR, seeds and patience50 /
max1000. No stable-KL substitution, clamping or seed/LR rescue is authorized.

Store predictions as float64 and check finiteness after serialization conversion;
finite float64 does not imply float32 representability. Shared scoring remains
unchanged. Replay rejects missing/different numerical identity and state precision.
Rerun all 30 TranSiGen identities from zero under the variant, preserving the
three historical successful float32 runs without using them in the new comparison.
The 15 completed cold-drug PRnet runs may be reused only with explicit immutable
import receipts proving original runtime551a3c9 plan/config/environment/source /
artifact identities, fresh-process exact sampled replay and scientific equivalence.
The remaining 45 learner identities require training; Mean reuse requires manifest
and scoring compatibility. None of this amendment establishes completed full-matrix
acceptance: ADR 0062 candidate-first publication and #231's 90-run checks still apply.
