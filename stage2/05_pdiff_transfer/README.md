# Step 05 — PerturbDiff transfer from Tahoe to tenv

## Retrieval evaluation

The active retrieval statistics are **AUPRC**, **AUPRC_positive**, and
**AUPRC_negative**, with corresponding **PR curve** panels. All use full-CM
DESeq2 significance and direction matching, under [ADR 0059](../../docs/adr/0059-stage2-significance-retrieval.md).
The [shared evaluator](../evaluation/README.md) defines the denominator,
undefined groups, pooling, and bounded resource policy.

- [PR curve](output/retrieval/PR_curve.png)
- [AUPRC and continuous metrics](output/retrieval/AUPRC.png)
- [Run statistics](output/retrieval/per_run.csv), [dose statistics](output/retrieval/per_dose.csv), and [summary](output/retrieval/summary.csv)
- [Replay manifest](output/retrieval/manifest.json)

From the repository root, after the relevant significance release is ready:

```bash
/home/chang/miniconda3/envs/arc/bin/python stage2/evaluation/publish.py --step 5
/home/chang/miniconda3/envs/arc/bin/python stage2/evaluation/evaluate.py replay --step 5
```

Previous retrieval scores, metric figures, and historical copies were removed;
Git retains their version history. Existing checkpoints and their original
validation-loss selections are unchanged. Other sections below describe input
and training provenance, rather than the retired retrieval publication.

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

Current results are linked above. Historical paired score differences remain only in Git.

Counts below are shared parent-connectivity blocks in the realized Tahoe
and tenv memberships. Mix uses Tahoe's training membership. Tahoe and
Transfer use the source training and validation memberships.

| Regime | Tahoe train ∩ tenv validation | Tahoe train ∩ tenv test | Tahoe validation ∩ tenv validation | Tahoe validation ∩ tenv test |
| --- | ---: | ---: | ---: | ---: |
| warm-start | 5 | 9 | 0 | 3 |
| cold-cm-drug | 1 | 2 | 0 | 0 |
| cold-drug | 0 | 0 | 0 | 0 |

## Shared cold tenv comparison

The original Step 04 tenv checkpoint and ordered prediction reuse are preserved.
The previous standalone cold replay score publication is retired. Current cold
regimes are included in the shared tenv inventory and replayed by the commands
above.

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

Each prediction directory contains `bridge.json`, the initialization
record used by Tahoe and Transfer. The record pins the historical EMA
checkpoint and the verified recipient request. It preserves the zero-dose
row and learned non-dose weights, remaps positive doses in log(µM), and
rebuilds the chemical fingerprint lookup and dose indices from tenv.
The loader checks the complete mapped state hash after loading it.
The source checkpoint stays in its original location and the bridge
stores the mapping and hashes needed to reconstruct the initialization.

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

Use the shared publication and replay commands above. Historical score releases
and their renderer arguments are retired; existing prediction, checkpoint,
training and exposure provenance remain available.
