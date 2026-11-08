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
