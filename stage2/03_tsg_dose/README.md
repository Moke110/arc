# TranSiGen dose-input optimization — round1

The active Round 1 report compares seven TranSiGen dose methods with
Mean_cell and ECFP4_MLP on the same ADR0077 axis-local G_5000 test rows, within
each Tahoe and tenv scenario and fold. It retains 84 TranSiGen fits and adds
12 independent ECFP4_MLP fits and 12 deterministic Mean_cell evaluations.
Pooled metrics weight outer-test samples equally. The overall figure shows
equal-fold estimates and sample standard deviations. The two baselines do not
participate in TranSiGen validation selection. Round 1 supplies comparison
evidence without a downstream model-weight dependency.

- [Overall metrics](figure/metrics.png) and [sample-pooled metrics](figure/data/pooled.csv)
- [Learning curves](figure/learning_curves.png) for the eight fitted methods
- [Cold-drug metrics by measured dose](figure/metrics_cold_drug_dose.png)

From this step directory, reproduce the figures with
`python figure/script/render.py --output <new-directory>`.
The comparison figures show Mean_cell first in dark gray and ECFP4_MLP
second in light gray. ECFP4_MLP is the first trained method in the learning
curves. The summary bar axes start at zero for Δ PCC, Δ Spearman and Δ RMSE.
The [presentation receipt](log/baselines-gray-acceptance.json) links to the
[baseline acceptance receipt](log/round1-baselines-acceptance.json), which
binds the frozen inputs, replayed fit origins and plot tables.

Mean_cell pools fit-partition treated expression and deduplicated original
matched vehicles from the same cell line across training doses. ECFP4_MLP
uses only the 2048-bit ECFP4 vector to predict treated expression with 128 ReLU
hidden units and a linear 5000-gene output. Both baselines subtract the same
per-test-sample matched original control as TranSiGen before scoring Δ.
ECFP4_MLP uses float32, Adam with learning rate 0.0001 and weight decay 0.00001,
batch size 64 and treated-expression MSE. Each fold restores its strict
validation-MSE minimum, with patience 50 and a maximum of 1000 epochs.

The [84-fit release](log/round1-axis-local-installation.json) and its
baseline successor are preserved under
[execution/historical/report_axis_local](execution/historical/report_axis_local/).
The installer's `recover` command verifies those inventories and the active
gray presentation in place.

[ADR0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md) defines the
paired 5k-local successor. Its separate input and bounded smoke entry points
are `execution/axis_local/prepare.py --view ../01_bases/output/view-5k
--candidate <fresh-directory>` and `execution/axis_local/smoke.py
--inputs <paired-freeze> --source <native-transigen-source>
--destination <fresh-smoke-directory> --device cuda`. The freeze contains
source-specific `tenv` and `tahoe_bulk` inputs and a paired digest manifest.
Its split ledger records both excluded strips of each cold-drug-dose fold.

New Stage 2 comparisons refer to the shared chemical representation as ECFP4.
The [ECFP4 contract migration](log/ecfp4-contract.md) binds the canonical fields
to the exact fingerprint and parent order of each accepted source freeze.

The smoke command runs one no-dose fit for each source and scenario and checks
fresh-process validation and test replay followed by a paired acceptance receipt.
Run both commands under `execution/capped.sh /usr/bin/time -v` with the resource
monitoring described below. Use `execution/axis_local/formal.py` for formal
planning, execution and scientific acceptance.

## Historical run protocol

[ADR0076](../../docs/adr/0076-dose-input-round1-exploration-boundary.md) and
[#240](https://github.com/Moke110/arc/issues/240) define this step as an
historical optimization exploration. The archived report retains **126 formal fits: 63 cold-drug and
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
They are distinct from the archived report's equal-weight fold×seed run means
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
intervals and undefined/constant counts are in [archived plotting data](execution/historical/report_126/data).

## Historical execution and acceptance

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

## Historical report contract

The [sealed 126-fit report](execution/historical/report_126/) contains
three PNGs with 2×4, 2×8 and 1×4 panels, plotting CSVs and its CSV-only renderer.
Its inputs use the historical full-18,164-gene denominator (`stage2-view-v2`).
The metrics sheet uses
nine fold×seed dots per fixed method/scenario, run-equal means and sample SD
(ddof=1); cold-drug selected-procedure dots follow locked fold choices. SD is
descriptive run spread, not a confidence interval. Native loss/MSE curves retain
every epoch/run; per-dose panels retain parent-equal/equal-seed OOF aggregation.

Accepted on2026-11-17 from saved evidence, without fitting or predicting:
540 run-metric rows (including selected-procedure references),60 run summaries,
15,977 complete curve coordinates and96 per-dose metric rows. An independent
isolated CSV-only render reproduced all three PNG SHA256s exactly; all33
original scientific CSVs and11 original plotting CSVs remain unchanged.
[Report acceptance](log/242-round1-report-acceptance.json) records every source,
accepted report hash and superseded historical PNG role/digest.

Reproduce from the sealed 126-fit tree using its `script/render.py` with
`--data <126-fit-tree>/data --output <new-dir>`, with
`execution/historical/report_126/` as the report tree. Its plotting CSVs preserve
run identities, deterministic jitter coordinates, source digests and counts.
Nonfinite dose metrics are excluded explicitly within parent means, undefined
parents within runs, and undefined runs within report means/SD; counts remain
visible and values are never replaced with zero. Sample SD is undefined with
fewer than two finite runs.

To reconstruct a historical 126-fit candidate, use `execution/publish_report.py prepare`
with `--matrix <saved-matrix> --candidate <new-report-directory>
--previous execution/historical/report_126 --scientific-receipt <126-origin-receipt>
--receipt <new-report-receipt>`. It verifies complete round1 acceptance and
saved metric digests, preserves historical plotting CSVs, and renders in an
isolated CSV-only directory. Inspect all three candidate PNGs, then run
`accept --candidate <candidate> --receipt <receipt> --visual-inspection`.
Tests: `python -m unittest discover -s execution/tests -p test_round1_report.py -v`.

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


## ADR0078 dual-source report release

The dual-source active report and the sealed 126-fit archive each have their
own scientific receipts and estimands.

After both 42-fit scenario receipts and the combined acceptance exist in the
formal candidate, use `execution/axis_local/publish_axis_local_report.py prepare`
with `--matrix`, `--candidate` and `--receipt`. The publisher verifies every
sealed scientific table and origin, the completed driver timing and cap
evidence, and the recorded WSL RAM and GPU samples. It reproduces all three
PNGs in an isolated CSV-only directory. Inspect each candidate PNG, then use
the same script's `accept --candidate ... --receipt ... --visual-inspection`.

Run `execution/axis_local/install_axis_local_report.py health --repo ...
--destination ...` in the materialized checkout. Use its `health.json` with
`install --candidate ... --receipt ... --health ...` in the deployment checkout.
Use a fresh candidate on the same filesystem as the target. Under the install
lock the complete previous `figure/` moves to
`execution/historical/report_126/` and the new candidate becomes `figure/`.
`log/round1-axis-local-installation.json` records both file inventories.
`recover --candidate ... --receipt ...` restores the previous active directory
if a swap stopped before that receipt, or verifies the installed state when the
receipt exists. The previous renderer and CSVs remain independently usable
from the archive.

`execution/axis_local/prepare.py` retains the exact bytes identified by the
accepted ADR0078 input manifest. Changes to that preparer need a successor
input generation and its own acceptance. The report installation consumes
existing sealed results and never regenerates inputs or fits.


## Baseline successor reproduction

With the verified paired ADR0078 input freeze, run
`execution/capped.sh /usr/bin/time -v -o <driver.time> python
execution/axis_local/baselines.py --inputs <paired-freeze>
--candidate <fresh-baseline-workdir> --device cuda` under the 10-second
RAM and GPU monitor. The driver serializes all 24 slots, verifies fold and
source identity, and runs exact prediction and per-sample score replay in a
fresh process. A rerun replays each finished slot and trains missing slots.

If fingerprint publication stopped before writing its manifest, first confirm
that the driver has exited. Move only `<baseline-workdir>/ecfp4/<source>` to a
new sibling path for audit, and rerun the driver, only if the source folder
has no `manifest.json` and `<baseline-workdir>/runs/<source>/` has no sealed
run manifest. An existing fingerprint manifest or sealed run requires
verification before any change.

Preserve the original timing-file modification times when collecting evidence.
Use `cp -p` for the driver timing file and per-unit `.time` files. The
publisher uses those timestamps to align run windows with monitor samples.

Use `execution/axis_local/publish_baselines.py prepare --matrix
<accepted-84-fit-workdir> --baselines <baseline-workdir>
--candidate <fresh-figure-candidate> --receipt <receipt.json>
--timing <driver.time> --monitor <resource-log-snapshot>`. It checks exact
outer-test identities, prior report hashes, fold coverage, metrics, and all
Mean_cell predictions against the independently implemented Stage 02 mean. It
reproduces the three PNGs in isolation using only the plot CSVs and renderers.
Inspect each candidate PNG, then run the same script's
`accept --candidate <figure-candidate> --receipt <receipt.json>
--visual-inspection`. Run the offline health checks in the materialized
checkout before deploying the accepted successor. Keep the candidate fit and
replay payloads in ignored checkpoints for audit.
