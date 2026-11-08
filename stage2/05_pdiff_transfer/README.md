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
