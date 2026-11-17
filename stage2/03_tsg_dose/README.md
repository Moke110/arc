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

`execution/run_matrix.py plan` lists all147 preregistered identities. Use
`run --inputs <freeze> --candidate <stable candidate> --source <native source>
--scenario cold-drug --device cuda` for the63 primary fits; `held-dose` and
`final` reuse the same frozen driver for63 and21 fits. Omit `--scenario` to
execute all147. The driver runs one capped subprocess at a time. The coordinator
must keep the independent10-second resource monitor active. `--limit-runs N`
pauses after N validation fits; `--resume` verifies all existing bindings and
replays complete origins before skipping. Scientific failures block release;
interrupted attempts are preserved before restarting the same scientific slot.

Per-fold21 restored validation checkpoints must replay before a choice receipt
is sealed. The receipt binds every validation manifest. Test access is logged
only after choice lock and each evaluated checkpoint replays in a new process.
Held-dose receipts are fixed comparison gates and do not select an auxiliary
winner. Final selection uses internal validation only. `accept` verifies the
same immutable input, source, configuration, code, attempts and phase evidence;
it cannot train or grant a partial147-fit public release. Candidate/report paths
remain stable until final publication, including during worktree cleanup.

The public evaluation gate independently recomputes the validation winner and
requires the evaluated checkpoint among the21 sealed origins. Validation scores
and parameter-count tie breaks are sealed before outer access. Nonfinite
evaluation or an exact-replay mismatch records a terminal scientific failure,
preserving the prior terminal and manifest in `failure.json`; resume cannot
retry that scientific slot. A process interruption has no scientific failure
classification. Opt-in native lifecycle checks use
`TSG_DOSE_LIFECYCLE_INPUTS=<freeze> TSG_DOSE_LIFECYCLE_SOURCE=<source>` with
`capped.sh python -m unittest discover -s execution/tests -p test_failure_lifecycle.py`.

Candidate `reports/<scenario>/` contains complete identity-checked OOF arrays,
per-well, parent–dose and parent metric tables, equal-seed summaries, per-dose
results, all curves/resources and validation selections. Primary inference is
selected procedure minus no-dose, with10,000 paired parent draws (seed2262026)
and all doses/seeds carried together. Choices stay fixed: intervals omit the
uncertainty from repeating method selection. Other candidate intervals are
descriptive. The final deployment choice has no new independent test score.

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

The standalone downstream reader is `tsg_dose.release.Release(directory, source)`.
Its `predict(seed, controls=..., morgan=..., doses=..., gene_ids=...,
profile_ids=..., sampling_seed=None)` returns a float32 DataFrame with the supplied
profile order and frozen gene columns. Supply Bases-normalized matched-control
profiles on exactly the released G_5000 order, binary radius2/chiral Morgan2048,
and positive µM doses. It needs no training matrices or original fit directory.
The default sampling seed is the checkpoint's recorded inference seed; exact
sampled replay requires identical input order and native batch64 boundaries.

A minimal accepted `tsg-dose-release-v1` directory contains `genes.csv`,
`input_contract.json`, `selection_receipt.json` and `checkpoints/{17,42,2026}.pt`.
`manifest.json` binds every payload digest, complete147-fit scenario counts,
the final method, three seeds, input-manifest digest, per-seed parameter counts
and environments. The reader validates all three checkpoint configurations,
saved dose transforms, runtime/source/backend identities and parameter counts.
It rejects incomplete coverage, smoke checkpoints, reordered genes, invalid
inputs and changed payloads. This API does not establish a release by itself:
formal matrix, downstream replay, figures and candidate-first acceptance remain
required before publication. The checkpoint role is internal-validation-selected,
with no new independent test result.

After all147 origins finish, run `execution/publish_release.py publish` under
`execution/capped.sh /usr/bin/time -v`, supplying `--inputs`, `--matrix`,
`--candidate` (a new empty release path), `--source`, `--destination` (the stable
`output/model`), and `--log`. Pass `--root` as the canonical checkout containing
the protected upstream products; `--figures` controls the report destination.
This command cannot train. It reruns full scientific acceptance, assembles only
the final winner's three restored checkpoints, and compares every final validation
row against the standalone reader in a fresh process before installing anything.
The installation receipt retains origin and choice audits and preservation hashes.

`figure/script/render.py --data figure/data --output figure` reads CSV only.
Publication first renders in an isolated directory containing just this script
and the CSVs. The six panels cover primary paired inference, primary per-dose
results, held-dose characterization, final validation, all147 learning curves
and resources. Held-dose tables retain primary-fold/final choice references;
held-dose outcomes do not select a method. A failed release gate leaves the prior
accepted model available. No formal release is asserted by these instructions.
