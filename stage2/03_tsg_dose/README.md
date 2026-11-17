# All-dose TranSiGen dose-input comparison

Independent implementation of [#232](https://github.com/Moke110/arc/issues/232)
and [ADR0075](../../docs/adr/0075-transigen-all-dose-input-selection.md).
All4,206 treated tenv wells,459 parents and0.2/1/10µM doses use released Bases
G_5000, original matched controls and Morgan2048 without refitting upstream.
Accepted on2026-11-17: all147 formal origins (63 cold-drug,63 held-dose,
21 deployment-selection) passed complete validation replay; all126 outer-test
origins also passed exact fresh-process test replay. There were147 unique
scientific attempts and no failures or interrupted fits. The released method is
**no-dose**, selected independently by final internal validation. Its restored
seed17/42/2026 checkpoints live at [output/model](output/model); the public
reader is [tsg_dose.release.Release](output/api/tsg_dose/release.py).

The primary cold-drug selection procedure chose scaling for folds0/1 and
input-numeric for fold2. Its OOF parent-equal ΔPCC was **0.384380**, versus
**0.362096** for no-dose: paired difference **0.022284**, conditional95% interval
**[0.012864,0.031775]** from10,000 parent draws. This interval holds the observed
validation choices fixed and omits uncertainty from repeating method selection.
It supports the observed selection procedure on this dataset; candidate-wise
outer results are descriptive and did not choose the deployment method.

| Cold-drug configuration | ΔPCC | ΔSpearman | ΔRMSE | ΔR² |
|---|---:|---:|---:|---:|
| no-dose | 0.362096 | 0.396073 | 0.162144 | -0.318000 |
| scaling | 0.391783 | 0.392280 | 0.154393 | -0.185941 |
| input-numeric | 0.363267 | 0.378867 | 0.158520 | -0.272532 |
| latent-numeric | 0.347783 | 0.388661 | 0.164662 | -0.386563 |
| input-embedding | 0.355776 | 0.379495 | 0.160509 | -0.310050 |
| latent-embedding | 0.361287 | 0.392171 | 0.161169 | -0.296807 |
| film | 0.355036 | 0.385770 | 0.164010 | -0.364452 |
| Selected procedure | 0.384380 | 0.392271 | 0.155887 | -0.220368 |

All values average three seed metrics, with repeats averaged within parent–dose
and doses/parents equally weighted. There are459 defined parents and no undefined
PCC or R² results in these summaries. Negative R² remains explicit; improved PCC
does not imply accurate expression reconstruction for every drug.

Held-dose characterization retains the primary/final method references without
reselection. Each entry below is parent-equal ΔPCC over459 parents and three
seeds. The1µM experiment interpolates;0.2/10µM extrapolate below/above the two
training doses.

| Fixed configuration | 0.2µM lower | 1µM interpolation | 10µM upper |
|---|---:|---:|---:|
| no-dose | 0.508014 | 0.512738 | 0.453084 |
| scaling | 0.473729 | 0.432670 | 0.300679 |
| input-numeric | 0.445045 | 0.498172 | 0.288563 |
| latent-numeric | 0.418558 | 0.498632 | 0.339332 |
| input-embedding | 0.483412 | 0.500017 | 0.397225 |
| latent-embedding | 0.447456 | 0.507954 | 0.314417 |
| film | 0.483933 | 0.495939 | 0.404203 |

No-dose had the largest point estimate at each held dose; its descriptive95%
parent-block intervals were [0.481114,0.534917], [0.488030,0.536953] and
[0.424265,0.481862], respectively. These intervals quantify variation among the
observed parents, not uncertainty over new dose levels. Only three measured dose
levels were studied; the auxiliary results establish neither broad dose-response
generalization nor a universally best architecture. Full four-metric estimates,
intervals and undefined/constant counts are in [figure/data](figure/data).

Final no-dose internal-validation ΔPCC was **0.425975** (three-seed mean).
This is a deployment selection score, with **no new independent test score**.
The separate downstream CUDA process reproduced all414 final-validation rows
and5,000 genes exactly for each released checkpoint, preserving profile/gene
order, transforms and sampled batch64 inference. The147-fit audit, seven choice
gates, hashes and670-file preservation check are retained in
[release acceptance](log/239-acceptance.json),
[downstream replay](log/239-acceptance-downstream.json),
[origin audit](log/239-acceptance-origins.csv) and
[verification](log/239-formal/verification.json). Bases, accepted02, frozenStage3
and immutable native source bytes were preserved.

The sequential primary/held-dose/final drivers took1h31m48s,1h39m59s and
33m44.88s; their maximum RSS values were1,934,036/1,905,088/1,946,076KiB.
Publication and full downstream acceptance took6m23.54s with1,461,944KiB maximum
RSS. All ran singly under10GiB RAM/8GiB swap caps with two threads and independent
10-second RAM/GPU monitoring. See the scenario logs and
[publication resources](log/239-formal/resource-summary.json).

Report outputs match step02's four-metric and native learning-curve conventions:
[metrics](figure/metrics.png), [paired differences](figure/paired_differences.png),
[all learning curves](figure/learning_curve.png), larger
[cold-drug](figure/learning_curve_cold-drug.png),
[held-dose](figure/learning_curve_held-dose.png) and
[final](figure/learning_curve_final.png) learning sheets, plus
[primary](figure/primary.png), [per-dose primary](figure/primary_dose.png),
[held-dose intervals](figure/held_dose.png),
[final validation](figure/final_validation.png),
[validation overview](figure/learning_curves.png) and
[resources](figure/resources.png). All12 actual PNGs passed visual inspection
against02; an independent isolated renderer reproduced their exact hashes from
the11 CSV tables alone. Native train loss uses the blue dashed left axis and
validation expression MSE the green solid right axis; their scales differ.


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
smokes do not enter the147-fit comparison; their evidence is separate from the
accepted formal release.

`execution/run_matrix.py plan` lists all147 preregistered identities. Use
`run --inputs <freeze> --candidate <stable candidate> --source <native source>
--scenario cold-drug --device cuda` for the63 primary fits; `held-dose` and
`final` reuse the same frozen driver for63 and21 fits. Omit `--scenario` to
execute all147. Always pass an explicit absolute `--source` path: the accepted
frozen driver's unused default retains a documented [ADR0075 exception](../../docs/adr/0075-transigen-all-dose-input-selection.md#frozen-driver-path-default-exception--2026-11-17).
The driver runs one capped subprocess at a time. The coordinator
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
and the CSVs. The report follows step02 for the four Δ indicator panels and native learning
curves. `metrics.png` shows PCC, Spearman, RMSE and R² for cold-drug and
held-dose with defined/undefined and constant-vector counts;
`paired_differences.png` shows all four primary/descriptive paired intervals.
`learning_curve.png` uses scenario-by-method panels with blue dashed native
training loss on the left axis and green solid validation MSE on the right,
retaining every epoch and fold/seed run. The three
`learning_curve_{cold-drug,held-dose,final}.png` sheets provide larger panels.
The shared Stage2 palette is recorded in `figure/data/plot_style.csv`, so
isolated rendering needs no upstream code or payloads. `metrics_values.csv`
retains saved summary values and count-source hashes;
`learning_curves_all.csv` retains both complete native loss and MSE coordinates.
The original six panels remain: primary paired inference, primary per-dose
results, held-dose characterization, final validation, validation-curve overview
and resources. Held-dose tables retain primary-fold/final choice references;
held-dose outcomes do not select a method. A failed release gate leaves the prior
accepted model and report figures available. CSVs and PNGs are staged before all
gates and installed with model/report rollback on replacement failure.
Held-dose intervals are descriptive95% percentile parent-block intervals from
10,000 draws (seed2262026), retaining every parent's doses and three seed metrics
together. They quantify variation across the observed parents, not uncertainty
across new dose levels or repeated model selection. Complete values, defined-parent
counts and undefined-draw counts are in `figure/data/held_dose_uncertainty.csv`.
The actual accepted release is pinned by [239-acceptance.json](log/239-acceptance.json).
