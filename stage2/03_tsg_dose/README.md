# All-dose TranSiGen dose-input comparison

Independent implementation of [#232](https://github.com/Moke110/arc/issues/232)
and [ADR0075](../../docs/adr/0075-transigen-all-dose-input-selection.md).
All4,206 treated tenv wells,459 parents and0.2/1/10µM doses use released Bases
G_5000, original matched controls and Morgan2048 without refitting upstream.
This step has no accepted formal release until63 primary+63 auxiliary+21 final
selection fits and exact fresh-process replay pass complete release acceptance.

`execution/prepare.py --view ../01_bases/output/view-5k --candidate <empty>`
freezes bounded float32 matrices and all partitions before training. The public
`tsg_dose.inputs.Inputs` reader verifies every digest and recomputes the frozen
identity-only partition rule. Cold-drug folds isolate parents and all doses;
held-dose folds retain each validation/test parent's other dose in training;
final uses a separate parent-held internal-validation split.

`execution/transigen.py --inputs <freeze> --candidate <empty> --source
../../pert_model/TranSiGen/source --device cuda` fits the native no-dose adapter,
restores the strict best validation-MSE checkpoint and replays validation in a
fresh process. Use `--smoke-rows 8 --max-epochs 2` only for a bounded smoke.
`--evaluate-only --selection <locked receipt>` accesses outer-test labels only
after21 formal validation fits and choice lock; smoke evaluation is explicit.
`--replay-only` never trains. Final checkpoints have internal validation only,
with no independent test result. Failed attempts remain evidence, never results.

Public scoring averages repeats within parent–dose before computing gene-wise
metrics, then weights doses and parents equally. Seed metrics are averaged,
never seed predictions. Constants and undefined correlations follow the shared
accepted scorer. Tests: `python -m unittest discover -s execution/tests -v`.

Foundation slice (#233): all-dose input preparation and bounded native no-dose
validation/test prediction, dose-aware scoring and exact fresh-process replay
passed on2026-11-17. See [233-acceptance.json](log/233-acceptance.json). These
smokes do not enter the147-fit comparison. All payloads remain candidates.

Four concatenation arms (#234) now preserve training-only dose transforms and
100-dimensional numeric/embedded conditioning at the input or both latent heads.
Behavioral fixtures verify dose use, ordered target-independent predictions,
finite gradients and disabled-dose native equivalence (float32 tolerance). Four
8-row/two-epoch CUDA smokes passed exact fresh-process validation and test replay;
see [234-smoke acceptance](log/234-smoke/acceptance.json). They are not formal fits.

Scaling and FiLM (#235) also passed bounded CUDA learning and exact fresh-process
validation/test replay. Scaling uses unstandardized log10(1+dose); FiLM starts
with identity chemical modulation and learns dose dependence with finite gradients.
Persisted-transform mismatch and nonfinite prediction checks fail closed. See
[235-smoke acceptance](log/235-smoke/acceptance.json) and
[rejection checks](log/235-smoke/rejections.json). These remain smoke evidence.
