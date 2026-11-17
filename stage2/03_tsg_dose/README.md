# TranSiGen dose-input optimization — round1

[ADR0076](../../docs/adr/0076-dose-input-round1-exploration-boundary.md) and
[#240](https://github.com/Moke110/arc/issues/240) define this step as an
optimization exploration. Round1 retains **126 formal fits: 63 cold-drug and
63 held-dose** across seven fixed methods, three folds and three seeds.
It supplies scientific comparison evidence and reports, with no downstream
model-weight dependency or deployment winner. All4,206 treated tenv wells,
459 parents and0.2/1/10µM doses use accepted Bases G_5000, matched controls and
Morgan2048 without refitting upstream. Accepted02 and Stage3 remain independent.

The scientific conclusion is limited to the four tested **100-dimensional
concatenation configurations**, which did not show stable prediction gains.
Scaling and FiLM are separate mechanisms with their own observations. This
fixed-configuration exploration does not establish universal architectural
superiority, and three measured doses limit dose-generalization claims.

## Retained scientific evidence

The values below are historical **parent-equal, equal-seed OOF summaries**.
They are distinct from the round1 report's equal-weight fold×seed run means
and descriptive sample SD. Seeds/folds are not independent biological replicates;
held-dose folds also correspond to different omitted doses.

The primary cold-drug selection procedure chose scaling for folds0/1 and
input-numeric for fold2. Its OOF parent-equal ΔPCC was **0.384380**, versus
**0.362096** for no-dose: paired difference **0.022284**, conditional95% interval
**[0.012864,0.031775]** from10,000 parent draws. This interval holds the observed
validation choices fixed and omits uncertainty from repeating method selection.
It supports the observed selection procedure on this dataset; candidate-wise outer results are descriptive and do not choose a deployment method.

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

Held-dose characterization is auxiliary and never selects a winner. Historical
CSV columns referring to final selection are audit context only. Each entry below is parent-equal ΔPCC over459 parents and three
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

## Active execution and acceptance

`execution/run_matrix.py plan` lists the126 round1 origins. `--scenario`
accepts only `cold-drug` or `held-dose`; omitting it covers both. A fresh run
uses `run --inputs <freeze> --candidate <new round1 candidate> --source
<absolute native source> --device cuda`. Inputs freeze and scientific learner
APIs retain their accepted identities. Validation-first cold-drug choice locks
all21 fold validation origins before test access; selected-procedure testing
uses those same restored checkpoints without extra fits. Held-dose choice gates
are fixed comparisons, never auxiliary winner selection.

`execution/transigen.py` also excludes final training. The strict best
validation-MSE checkpoint restoration, sealed selection binding, nonfinite
failure handling and exact fresh-process validation/test replay remain required.
Apply [resource limits](../../docs/agents/execution-environment.md), run formal
fits singly under `execution/capped.sh`, and monitor RAM/GPU/WSL every10 seconds.
Round1 report publication reuses saved predictions and metrics; it never trains
or re-predicts.

`execution/run_matrix.py accept --inputs <freeze> --candidate <matrix>
--source <absolute native source> --device cuda --receipt <new receipt>` verifies all126 retained
origins and their existing exact replay proofs. It accepts the original147-fit
plan through the preserved-source mapping without requiring final origins or
weight-reader replay, leaves historical scientific tables untouched, and writes
a separate `round1-acceptance-cold-drug-held-dose.json`. New126-fit candidates
have their own plan/code identity and produce their scientific tables only after
retained origins pass. Partial scenario receipts are not full round1 acceptance.
Tests use `python -m unittest discover -s execution/tests -p test_round1.py -v`;
existing input-freeze, selection, learner/replay and scientific comparison checks
remain applicable.

## Report contract

The active report contains exactly `round1_metrics.png` (2×4),
`round1_learning_curves.png` (2×8), and
`round1_metrics_cold_drug_dose.png` (1×4), with plotting CSVs in `figure/data`
and a CSV-only renderer in `figure/script`. Candidate-first figure acceptance
replaces the previous report only after all checks pass. The metrics sheet uses
nine fold×seed dots per fixed method/scenario, run-equal means and sample SD
(ddof=1); cold-drug selected-procedure dots follow locked fold choices. SD is
descriptive run spread, not a confidence interval. Native loss/MSE curves retain
every epoch/run; per-dose panels retain parent-equal/equal-seed OOF aggregation.

## Historical147-fit release

The completed21 final fits, old no-dose seed17/42/2026 checkpoints at
`output/model`, original147-fit acceptance and historical scientific summaries
remain unchanged for traceability. They are excluded from round1's scope and
future dependencies. `execution/publish_release.py` rejects weight publication
and downstream weight-reader replay. `tsg_dose.release.Release` is retained
solely to reconstruct the old release, not as an active downstream contract.

Historical audit records: [release acceptance](log/239-acceptance.json),
[downstream replay](log/239-acceptance-downstream.json),
[origin audit](log/239-acceptance-origins.csv), and
[verification](log/239-formal/verification.json).
The [role mapping](log/round1-role-mapping.json) resolves original frozen driver
hashes to byte-identical source snapshots. The original147-fit registration,
21 final fits, released weights, sealed records and required reconstruction
payloads stay at their original paths. Old PNG roles are superseded only after
round1 candidate report acceptance; sealed historical acceptance digests retain
their historical meaning.
