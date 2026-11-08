# Step 05 — PerturbDiff transfer from Tahoe to tenv

> **Cross-dataset dose encoding limitation.** Tahoe uses 0.05, 0.5 and 5 µM.
> Tenv uses 0.2, 1 and 10 µM. The current bridge interpolates the learned
> categorical embeddings on the log-dose axis and clamps outside Tahoe's
> range. Zero-shot Tahoe receives its 5 µM embedding for tenv 10 µM.
> Transfer starts from the same mapping and can update it with tenv data.
> Improving cross-dataset dose encoding is a priority for follow-up.

[ADR 0047](../../docs/adr/0047-perturbdiff-tahoe-transfer-to-tenv.md) defines the
warm-start, cold-cm-drug and cold-drug comparisons. Each tests tenv responses
and compares CM, Tahoe, Transfer and Mix with two deterministic references.

## Formal results

All 36 learned-arm slots completed. Transfer scored below CM in all nine
paired seed comparisons on these fixed splits. Each pair subtracts CM's
mean of positive and negative AP@100 from Transfer's corresponding mean,
with equal weight for each direction. Seed variation describes repeated
fits on each fixed split.

| Regime | Learned slots | Seed 17 | Seed 42 | Seed 2026 | Complete pairs | Mean difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| warm-start | 12/12 | -0.01737 | -0.03511 | -0.03219 | 3/3 | -0.02822 |
| cold-cm-drug | 12/12 | -0.02825 | -0.02096 | -0.01864 | 3/3 | -0.02262 |
| cold-drug | 12/12 | -0.02966 | -0.02252 | -0.02481 | 3/3 | -0.02566 |

Counts below are shared parent-connectivity blocks in the realized Tahoe
and tenv memberships. Mix uses Tahoe's training membership. Tahoe and
Transfer use the source training and validation memberships.

| Regime | Tahoe train ∩ tenv validation | Tahoe train ∩ tenv test | Tahoe validation ∩ tenv validation | Tahoe validation ∩ tenv test |
| --- | ---: | ---: | ---: | ---: |
| warm-start | 5 | 9 | 0 | 3 |
| cold-cm-drug | 1 | 2 | 0 | 0 |
| cold-drug | 0 | 0 | 0 | 0 |

- [Eight-metric summary](summary.csv), [per-run scores and valid-group counts](per_run.csv), and [metric plots](metrics.png). The cold CM and deterministic references share the `cold-shared` scope.
- [Full-precision paired contrasts](execution/formal_record/paired_ap100.json) and [undefined metric counts](execution/formal_record/undefined_reasons.json).
- [Directional precision–recall curves](pr_curves.png) and [plotted points](execution/formal_record/pr_points.json).
- [Learning curves](learning_curves.png), [plotted points](execution/formal_record/learning_points.json), and [selected checkpoint index](execution/formal_record/checkpoints.json).
- [Publication checksums](execution/formal_record/publication.json) and [execution and replay evidence](execution/formal_record/run.json).
- [Output handoff](output/handoff.json) maps all 36 logical slots to the sealed prediction values, case order and selected source and target checkpoints. Paths start at the repository root. Both cold regimes point to the same three historical CM predictions. The handoff references the original local files in place and retains their SHA-256 digests.

## Shared cold tenv comparison

Both cold regimes use the frozen Step 04 tenv fold 0. The comparison CLI
verifies all three membership lists, their parent-connectivity blocks and
case order, the matched controls and the ordered 18,164-gene axis. It verifies
the three indexed CM checkpoints and their existing held-out predictions
against the Step 04 recipe, validation replay and prediction seals. Scoring
reads held-out expression only after these checks. The two seedless
references and the three CM predictions are scored once. The replay checks
all eight metrics and both directional PR curves against the Step 04 record.

Run from the repository root using the `arc` interpreter. Submit the full
replay to the dedicated scheduling subagent as a CPU job under the repository
resource limits.

```bash
step04=stage2/04_pdiff_hyperparameter
step05=stage2/05_pdiff_transfer
python "$step05/execution/script/cold_tenv.py" \
  "$step04/execution/output/tenv-formal-ticket72" \
  "$step04/output/checkpoints/tenv" \
  "$step05/execution/output/cold-tenv-replay"
```

`comparison.json` binds the two regime names to the shared artifacts, with
CM evidence indexed by seed and a SHA-256 inventory of the output.
`split.json` gives the frozen input identity, predictor paths and hashes,
case and group counts, and each role's parent-connectivity blocks.
`predictor_request.json` is the verified Step 04 request and lists only
predictor inputs. `per_run.csv`, `per_group.csv` and `tenv/pr_curves.json`
use the Step 04 evaluator format, where the model name `perturbdiff` means
the CM arm for this comparison.

The [scheduled replay record](execution/cold_record/replay.json) confirms
exact agreement with the Step 04 fold 0 scores and directional PR curves.
The frozen train, validation and test roles contain 2,388, 417 and 1,401
cases, respectively. The CPU replay finished in 14.50 seconds and peaked
at 459,572 KiB RSS. The replay record binds the local comparison and raw
scheduler evidence by SHA-256.

## Globally cold Tahoe source release

The cold-drug source release removes the tenv fold 0 validation and test
parent-connectivity blocks from Tahoe fold 0 training and validation.
The release retains the original order of eligible Tahoe wells, their
matched-control indices, and their float32 expression deltas on the ordered
18,164-gene axis. Tahoe test case metadata retains the frozen membership and
supports the native conditioning vocabulary. Only the filtered training and
validation partitions carry expression targets.

Submit the build, smoke fit and reload as separate jobs to the dedicated
scheduler. These commands use the Step 04 Tahoe formal fold 0 source request
and the verified Step 05 cold tenv split.

```bash
step05=stage2/05_pdiff_transfer
pdiff=stage2/_envs/perturbdiff/bin/python
cli="$step05/execution/script/cold_tahoe.py"
output="$step05/execution/output"
release="$output/cold-tahoe-ticket79"
smoke="$output/cold-tahoe-smoke-ticket79"
"$pdiff" "$cli" build "$release"
"$pdiff" "$cli" smoke-fit "$release" "$smoke"
"$pdiff" "$cli" reload "$release" "$smoke" "$output/cold-tahoe-smoke-replay-ticket79.json"
```

`release.json` records the frozen identities, forbidden blocks, retained and
excluded counts and the parent overlap with each tenv role before and after
filtering. `verify` reconstructs the expected Tahoe membership and checks
the request, predictor files and all three unchanged Step 04 formal recipe
files. Cold-drug Mix can use the verified request's training files, preserving
Tahoe groups and control pairing. `formal-fit` takes the same release and a
fresh run directory plus `--seed 17`, `--seed 42` or `--seed 2026`.

The [qualification record](execution/cold_tahoe_record/qualification.json)
covers 198 forbidden tenv parent blocks. Filtering removed 500 training wells
and 450 training groups across three Tahoe parent blocks. Validation already
had zero overlap. The release retains 35,936 training wells and 6,389
validation wells. The 12-step smoke selected step 12 and independently
reloaded its EMA checkpoint with zero prediction difference on the full
gene axis. The verifier rejected an injected forbidden case in each fitted
role after the case file checksum was recomputed.

## Warm-start tenv comparison

`warm_tenv.py prepare` uses split seed 20261105 to assign each parent one
test dose and to balance the number of test parent–dose combinations
across doses. It assigns approximately 15% of the remaining combinations
to validation, up to one per parent, and keeps a training combination for
each parent. All aliases and repeat wells stay with their parent–dose
combination. The manifest records the complete assignment. Prediction and
scoring retain the frozen drug–dose groups, case order, matched controls,
Morgan fingerprints and 18,164-gene axis.

Run the following commands from the repository root through the dedicated
scheduling subagent. Use the `arc` interpreter for prepare and score. Use
the isolated PerturbDiff interpreter for fit and predict. The following
replays the bounded check with 12 optimizer steps and seed 17.

```bash
step04=stage2/04_pdiff_hyperparameter
step05=stage2/05_pdiff_transfer
warm="$step05/execution/script/warm_tenv.py"
gpu=stage2/_envs/perturbdiff/bin/python
out="$step05/execution/output/warm-tenv-smoke"
python "$warm" prepare "$step04/execution/output/tenv-formal-ticket72" "$out" --smoke
digest=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$gpu" "$warm" fit "$out" 17 --plan-sha256 "$digest"
"$gpu" "$warm" predict "$out" 17 --plan-sha256 "$digest"
python "$warm" score "$out" "${out}-scores" --plan-sha256 "$digest"
```

For a formal run, prepare a fresh output without `--smoke` and complete the
fit and prediction phases for seeds 17, 42 and 2026 before scoring. The
three byte-identical Step 04 fold 0 recipes are frozen in both modes.
`verify` takes the same output and plan digest as the workers and checks the
recipe, canonical membership and predictor allowlist. Fit and prediction
seal separate outputs and reproduce the validation-selected EMA loss and
probe. Only the scorer reconstructs warm-start held-out treated expression,
after verifying the prediction and checkpoint evidence. Score outputs use
the Step 04 format and include all eight metrics, directional PR curves,
and the seedless `no_change` and training-group `delta_mean` references.

The [scheduled smoke record](execution/warm_record/replay.json) binds the
split, recipe, fit, predictions, scores and resource telemetry by SHA-256.
The split holds out 153 parent–dose combinations at 0.2 µM and 152 at
each of 1 and 10 µM. Train, validation and test contain 2,387, 417 and
1,402 wells. The seed 17 fit selected step 12 and reproduced the validation
loss and the full-axis prediction probe exactly on reload. Peak host RSS
was 3,155,924 KiB during fitting and 2,597,812 KiB during held-out prediction.

## Directly pooled Mix

Mix starts from random initialization and uses the matching CM recipe and
native group sampler. Training pools the two verified training partitions
with source-prefixed groups. Each sampled case keeps its original matched
control and expression delta. The dose vocabulary comes from the union of
the training doses. Validation and test cases come from tenv. Tahoe
validation and test files are outside the Mix predictor allowlist. In
cold-drug, Tahoe training comes from the filtered release above.

Run these commands from the repository root through the dedicated scheduler.
Prepare and score use the `arc` interpreter. Fit and predict use the pinned
PerturbDiff interpreter. Prediction independently reloads the selected EMA
checkpoint and reproduces its tenv validation loss and prediction probe.

```bash
step05=stage2/05_pdiff_transfer
cli="$step05/execution/script/mix.py"
pdiff=stage2/_envs/perturbdiff/bin/python
out="$step05/execution/output/mix-smoke"
python "$cli" prepare "$out" --smoke
digest=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" fit "$out" cold-drug 17 --plan-sha256 "$digest"
"$pdiff" "$cli" predict "$out" cold-drug 17 --plan-sha256 "$digest"
python "$cli" score "$out" cold-drug "${out}-scores" --plan-sha256 "$digest"
```

`plan.json` records all nine formal slot identities, source and case
inventories, dose categories and CM recipe hashes. For a formal run, prepare a
fresh directory without `--smoke`, run fit and predict for seeds 17, 42 and
2026 in each regime, then score each regime separately. `verify` checks
all three pools against the frozen releases. A score seals all eight metrics,
both directional PR curves and hashes for the held-out evidence.

The cold-drug qualification used seed 17 and 12 optimizer steps. It selected
step 12 by native EMA loss on the 138 tenv validation groups and reproduced
the validation probe exactly on reload. It predicted all 1,401 held-out tenv
cases across 462 groups. Peak host RSS was 3,300,484 KiB for fitting and
2,648,996 KiB for held-out prediction.
The [qualification record](execution/mix_record/qualification.json) binds the
plan, checkpoint, predictions, score outputs and scheduler evidence.
Independent scoring reproduced every metric, both directional PR curves and
the evaluation manifest byte-for-byte.

## Tahoe zero-shot on tenv

Tahoe reuses the validation-selected Step 04 Tahoe fold 0 EMA checkpoint
for each seed. The source verifier checks the historical producer revision,
the exact source requests and recipes, the train and validation memberships,
and the original checkpoint and replay evidence. Each prediction records
parent-connectivity overlap between every Tahoe role and every tenv role.
For interpretation of source exposure, consult the Tahoe training and
validation rows separately.

Run the source check, prediction and score as separate scheduled jobs.
The example qualifies seed 17 on the full tenv gene axis for both regimes.

```bash
step04=stage2/04_pdiff_hyperparameter
step05=stage2/05_pdiff_transfer
cli="$step05/execution/script/tahoe_zero_shot.py"
pdiff=stage2/_envs/perturbdiff/bin/python
out="$step05/execution/output/tahoe-zero-shot"
python "$cli" verify-source "$step04/execution/output/tahoe-formal-ticket73" "$out/source.json"
source_sha=$(sha256sum "$out/source.json" | cut -d ' ' -f 1)

for regime in warm-start cold-cm-drug; do
  if [ "$regime" = warm-start ]; then
    recipient="$step05/execution/output/warm-tenv-ticket77"
    plan=plan.json
  else
    recipient="$step05/execution/output/cold-tenv-ticket76"
    plan=split.json
  fi
  recipient_sha=$(sha256sum "$recipient/$plan" | cut -d ' ' -f 1)
  "$pdiff" "$cli" predict "$out/source.json" "$recipient" "$regime" 17 "$out/$regime/seed-17" \
    --source-sha256 "$source_sha" --recipient-sha256 "$recipient_sha"
  python "$cli" score "$out/source.json" "$recipient" "$regime" "$out/$regime" "$out/$regime-scores" \
    --source-sha256 "$source_sha" --recipient-sha256 "$recipient_sha" --seed 17
done
```

For the complete score output, predict seeds 17, 42 and 2026 and omit the
scorer's `--seed` argument. The scorer verifies each prediction seal and
case order before opening the test truth. It emits all eight group-weighted
metrics and both directional PR curves.

Each prediction directory contains `bridge.json`, the initialization
record used by Tahoe and Transfer. The record pins the historical EMA
checkpoint and the verified recipient request. It preserves the zero-dose
row and learned non-dose weights, remaps positive doses in log(µM), and
rebuilds the chemical fingerprint lookup and dose indices from tenv.
The loader checks the complete mapped state hash after loading it.
The source checkpoint stays in its original location and the bridge
stores the mapping and hashes needed to reconstruct the initialization.

The [scheduled qualification record](execution/tahoe_zero_shot_record/qualification.json)
verifies all three source checkpoints, each selected at step 8,000. The seed
17 prediction preserved 1,402 warm-start wells in 460 groups and 1,401
cold-cm-drug wells in 462 groups. All eight metrics were defined for every
group in both runs. Nine warm-start and two cold-cm-drug test parent blocks
appeared in Tahoe training. Tahoe validation included three warm-start test
parent blocks and zero cold-cm-drug test parent blocks. Both predictions
completed concurrently with a workstation VRAM peak of 3.30 GiB. Maximum
host RSS was 2,445,856 KiB for warm-start and 2,446,240 KiB for cold-cm-drug.

## Transfer and globally cold Tahoe

The Transfer plan freezes the tenv request and matching CM recipe for each
regime and seed. Warm-start and cold-cm-drug pin their selected Step 04 Tahoe
checkpoints. Cold-drug pins the filtered Tahoe request, source recipe and
per-seed fit location. After that source fit, the native validation replay
and fit seal establish the exact selected source checkpoint.

Run the following phases under the scheduler. Within each phase, the seed
slots are independent. The example shows one formal seed.

```bash
step=stage2/05_pdiff_transfer
cli="$step/execution/script/transfer.py"
pdiff=stage2/_envs/perturbdiff/bin/python
out="$step/execution/output/transfer-formal"
python "$cli" prepare "$out"
plan_sha=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
python "$cli" verify "$out" --plan-sha256 "$plan_sha"

# Run only for the globally cold Tahoe source, for seeds 17, 42 and 2026.
"$pdiff" "$cli" source-fit "$out" 17 --plan-sha256 "$plan_sha"

# Repeat these phases for each regime and seed after its source is ready.
regime=cold-drug
"$pdiff" "$cli" tahoe-predict "$out" "$regime" 17 --plan-sha256 "$plan_sha"
"$pdiff" "$cli" fit "$out" "$regime" 17 --plan-sha256 "$plan_sha"
"$pdiff" "$cli" predict "$out" "$regime" 17 --plan-sha256 "$plan_sha"

# Score a formal regime after its three predictions per arm are sealed.
python "$cli" tahoe-score "$out" "$regime" "$out/$regime-tahoe-score" --plan-sha256 "$plan_sha"
python "$cli" score "$out" "$regime" "$out/$regime-score" --plan-sha256 "$plan_sha"
```

`prepare --smoke` activates seed 17 with a 12-step target recipe. The
globally cold smoke route uses the separately qualified filtered Tahoe
checkpoint from ticket #79. Formal plans require their own filtered source
fit with the matching formal Tahoe recipe. Both stages can reuse the sealed
ticket #78 Tahoe predictions for seed 17 in the two historically sourced
regimes.

Transfer reads the bridge record from the corresponding Tahoe prediction and
strict-loads its mapped model state before constructing the optimizer,
scheduler and EMA. The target run records the source and recipient request
hashes, source checkpoint and history hashes, target recipe hash, bridge and
mapped-state hashes, fresh training state, and a gradient on the mapped
positive dose rows. The source and target validation histories remain in
their separate runs. The Transfer worker verifies case and bridge metadata
before fitting or predicting. Scoring opens the tenv test truth only after
checking the selected predictions and their seals.

The formal cold source fits are under `source/seed-*/fit`. New Tahoe
predictions and their bridge records are under `tahoe/<regime>/seed-*`.
Transfer fits and predictions are under `slots/<regime>/seed-*`. Successful
fit and predict commands verify existing sealed outputs on retry. Failed
attempts stay in their attempt directories.

The [scheduled qualification](execution/transfer_record/qualification.json)
records the full-axis cold-cm-drug Transfer fit and replay, the globally
cold Tahoe zero-shot prediction, and the two seed-17 evaluations.

## Fixed formal matrix

The dedicated scheduler uses `formal_matrix.py` from the integrated checkout.
The plan freezes 36 learned-arm slots on fold 0, with 24 new fits, 30 new
predictions and six logical references to the three historical cold CM
predictions. It binds each case and group order, the worker plans and code,
recipes, predictor files and source lineages before dispatch. Both cold
regimes share one target membership. Warm CM's independently prepared files
must match the qualified Transfer and Mix target data byte-for-byte.

```bash
step=stage2/05_pdiff_transfer
cli="$step/execution/script/formal_matrix.py"
out="$step/execution/output/formal-ticket82"
python "$cli" prepare "$out"
digest=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
python "$cli" verify "$out" --plan-sha256 "$digest"
python "$cli" preflight "$out" --plan-sha256 "$digest"
python "$cli" run "$out" --plan-sha256 "$digest"
python "$cli" verify-inventory "$out" --plan-sha256 "$digest"
```

Use the `arc` interpreter for these commands. The coordinator launches workers
with the pinned PerturbDiff interpreter under the Step 04 scheduler. Start
with at least 30 GiB free disk and 12 GiB available RAM. Training is serial.
Up to two predictions may overlap when both resource reservations fit.

During a run, use `status` with the same root and digest. To relay a requested
pause or resume, use `pause` or `resume` with those arguments and a
`--reason` describing the instruction. Pause holds dispatch across phases and
lets active workers drain under the scheduler's pressure limits. After an
explicit resume, rerun `run` if the coordinator already exited.

`phases/` contains immutable ready queues, scheduler status, telemetry,
admission decisions and attempt logs. Each failed or interrupted worker may
retry once with the same command. A prediction retry reads the sealed fit.
Exhausted prerequisites block only dependent jobs. When every slot has a
terminal outcome, `prediction_inventory.json` binds each prediction and
checkpoint by path and hash, accounts for failed slots and references reused
CM evidence. The inventory is the input to publication and scoring.

## Publication replay

The publisher verifies the frozen plan, all terminal outcomes, selected
checkpoint evidence and held-out prediction files before reading test
truth. It computes the eight Step 04 metrics on equally weighted tenv
drug–dose groups. Cold CM and deterministic references are scored once
under `cold-shared`. The two cold regimes link to those rows and have their
own Tahoe, Transfer and Mix results. Each seed pair uses the difference
between the two arms’ equal means of positive and negative AP@100. The
regime mean requires all three pairs.

Use the `arc` interpreter from the repository root. Submit real scoring
and replay to the dedicated scheduler. All output directories must be fresh.

```bash
step=stage2/05_pdiff_transfer
cli="$step/execution/script/publish.py"
out="$step/execution/output"
formal="$out/formal-ticket82"
plan=$(sha256sum "$formal/plan.json" | cut -d ' ' -f 1)
inventory=$(sha256sum "$formal/prediction_inventory.json" | cut -d ' ' -f 1)

python "$cli" score "$formal" "$out/scores" --plan-sha256 "$plan" --inventory-sha256 "$inventory"
scores=$(sha256sum "$out/scores/scores.json" | cut -d ' ' -f 1)
python "$cli" render "$out/scores" "$out/publication" --scores-sha256 "$scores"
release=$(sha256sum "$out/publication/publication.json" | cut -d ' ' -f 1)

python "$cli" score "$formal" "$out/scores-replay" --plan-sha256 "$plan" --inventory-sha256 "$inventory"
python "$cli" verify-replay "$out/scores-replay" "$out/publication" --publication-sha256 "$release"
```

Use the published plan and inventory digests when verifying a formal run.
Replay verifies the score record and reproduces the publication checksums,
including plotted points and PNG bytes. The full scoring output includes
group-level values and reasons for undefined metrics. The compact records
include per-run coverage and reasons, both directional PR curves, training
and validation points, and checkpoint hashes and paths. The checkpoint index
links each logical slot to its source and target fits by digest. Tahoe and
Transfer link to the same source fit for a given regime and seed. The two
historically sourced regimes also share those source records.

To verify the tracked release from a checkout, run the command below. It
uses the publication checksum and path map in
[run.json](execution/formal_record/run.json). The map includes the sealed
[generated report](execution/formal_record/results.md) under its original
logical name, `README.md`.

```bash
step=stage2/05_pdiff_transfer
python "$step/execution/script/publish.py" verify-release "$step"
```

Add `--scores "$step/execution/output/ticket83-scores"` to render a
temporary publication from the original local score record and compare
the resulting seal with the tracked release. An independently rescored
directory can be supplied the same way.
