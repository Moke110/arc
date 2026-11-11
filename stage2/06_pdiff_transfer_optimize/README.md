# Step 06 — PerturbDiff transfer diagnostics

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
/home/chang/miniconda3/envs/arc/bin/python stage2/evaluation/evaluate.py score --step 6
/home/chang/miniconda3/envs/arc/bin/python stage2/evaluation/evaluate.py replay --step 6
```

Previous retrieval scores, metric figures, and historical copies were removed;
Git retains their version history. Existing checkpoints and their original
validation-loss selections are unchanged. Other sections below describe input
and training provenance, rather than the retired retrieval publication.

Step 05 observed a Tahoe Transfer deficit against the tenv-only CM arm.
Step 06 establishes a shared validation baseline and tests positive-dose
initialization and target base learning rate in separate paired arms.
[ADR 0049](../../docs/adr/0049-perturbdiff-transfer-optimization-comparisons.md)
extends the study to five optimizations on the Step 05 cold-cm-drug split.
[ADR 0048](../../docs/adr/0048-perturbdiff-transfer-diagnosis.md)
records the screening, confirmation and attribution rules.

## Paired optimization results

Run `publish` through the dedicated scheduler to verify all 18 cells
against their sealed artifacts and reproduce the publication. From the
runtime checkout, use the arc interpreter.

```bash
step=stage2/06_pdiff_transfer_optimize
python "$step/execution/script/paired_publication.py" publish "$step/output/paired-transfer"
```

The tracked numeric record also supports a figure-only replay. With the
arc plotting environment, the command verifies the regenerated PNG
against its recorded SHA-256. Supply `--figure /path/to/replay.png` to
write a separate copy.



## Historical cold-drug results

The dose reset improved the 0.2 and 1 µM strata by +0.01155 and +0.01470
respectively, with a −0.00080 change at 10 µM. Neither intervention reached
the +0.010 pooled promotion threshold. The added screening seeds and the
fold 1 source fits, target fits and test scoring were skipped. The dose
reset reduced the seed 17 validation gap by about 54%. The major-contribution
rule requires independent fold 1 replication and a matched gap comparison
there, so attribution remains unresolved.

The separate source audit covers nine shared tenv training parents and
50 Tahoe cell lines. Seven of those parents share identical standardized
parent SMILES. Historical Tahoe training overlapped one tenv validation
parent. After global exclusion, both Tahoe fitting roles overlap zero tenv
validation or test parents. Dose brackets and the assay, time, normalization
and coverage limitations are documented in the [source audit](execution/source_audit.md).

## Validation baseline

Use the arc Python interpreter for verification and scoring. Use the pinned
PerturbDiff interpreter for prediction. Submit these phases to the dedicated
scheduler. The baseline reads sealed predecessor artifacts in place.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/baseline.py"
pdiff=stage2/_envs/perturbdiff/bin/python
out="$step/output/baseline/fold-0-seed-17"

python "$cli" verify "$out" --seed 17
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" predict "$out" --seed 17 --plan-sha256 "$plan"
python "$cli" score "$out" --seed 17 --plan-sha256 "$plan"
```

The plan pins the selected EMA states and the common target-free validation
request. Both arms use the same validation well and group order, same-plate
controls, released Morgan fingerprints, and inference seed and batch layout.
The score phase joins frozen validation Δ after verifying both prediction
seals. `scores/` contains per-group directions and undefined reasons,
per-dose and pooled summaries, and the paired CM minus Transfer gaps. The
same commands accept seeds 42 and 2026 with their own output directories.
Use `--runtime /absolute/path/to/arc` to locate predecessor artifacts and the
pinned interpreter from a separate checkout.

## Paired Transfer reference

[ADR 0049](../../docs/adr/0049-perturbdiff-transfer-optimization-comparisons.md)
uses the Step 05 cold-cm-drug fold 0 original Transfer as a common reference
at seeds 17, 42 and 2026. Verify, predict and score each seed as separate jobs
through the dedicated scheduler. From the runtime checkout, use the arc
interpreter for verification and scoring and the pinned PerturbDiff
interpreter for prediction.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/paired_reference.py"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
seed=17  # Repeat for 42 and 2026.
out="$step/output/paired-transfer/reference/seed-$seed"

python "$cli" verify "$out" --seed "$seed"
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" predict "$out" --seed "$seed" --plan-sha256 "$plan"
python "$cli" score "$out" --seed "$seed" --plan-sha256 "$plan"
```

## Random positive-dose reset

The `random` arm restores the three native tenv positive-dose embeddings
after loading the same-seed Tahoe bridge. It then uses the original tenv
training recipe and selects its own EMA by the author validation diffusion
loss. Run all three seeds through the dedicated scheduler.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/paired_target.py"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
seed=17  # Repeat for 42 and 2026.
reference="$step/output/paired-transfer/reference/seed-$seed"
out="$step/output/paired-transfer/random/seed-$seed"
reference_plan=$(sha256sum "$reference/plan.json" | cut -d ' ' -f 1)

python "$cli" verify "$out" --arm random --seed "$seed" --reference-sha256 "$reference_plan"
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" fit "$out" --arm random --seed "$seed" --plan-sha256 "$plan"
"$pdiff" "$cli" predict "$out" --arm random --seed "$seed" --plan-sha256 "$plan"
python "$cli" score "$out" --arm random --seed "$seed" --plan-sha256 "$plan"
```

The fit seal and selected checkpoint are under
`output/checkpoints/04-paired-random-fold-0-seed-$seed/fit/`.
The fit record includes the restored dose row digests, actual optimizer
steps and selected EMA step. The replay writes `prediction/`, the common
validation evaluator writes `scores/`, and `paired.json` records the
per-dose differences from original Transfer at the same seed.

## Continuous log dose

The `log` arm fits Tahoe on the original fold 0 source memberships and
recipe. The dose encoder applies a trainable affine map to
`log10(d / 0.05 µM)` and uses a separate trainable control vector at 0 µM.
Tenv starts from the same seed's loss-selected log Tahoe EMA and rebuilds
its own dose and Morgan lookups. The encoder evaluates 10 µM directly.

```bash
step=stage2/06_pdiff_transfer_optimize
source_cli="$step/execution/script/paired_log_source.py"
target_cli="$step/execution/script/paired_target.py"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
seed=17  # Repeat for 42 and 2026.
reference="$step/output/paired-transfer/reference/seed-$seed"
out="$step/output/paired-transfer/log/seed-$seed"
reference_plan=$(sha256sum "$reference/plan.json" | cut -d ' ' -f 1)

"$pdiff" "$source_cli" fit "$out" --seed "$seed" --reference "$reference" --plan-sha256 "$reference_plan"
python "$target_cli" verify "$out" --arm log --seed "$seed" --reference-sha256 "$reference_plan"
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$target_cli" fit "$out" --arm log --seed "$seed" --plan-sha256 "$plan"
"$pdiff" "$target_cli" predict "$out" --arm log --seed "$seed" --plan-sha256 "$plan"
python "$target_cli" score "$out" --arm log --seed "$seed" --plan-sha256 "$plan"
```

The dedicated scheduler runs the jobs in order. Source and target fit
evidence lives under `output/checkpoints/08-paired-log-tahoe-fold-0-seed-$seed/fit/`
and `output/checkpoints/05-paired-log-fold-0-seed-$seed/fit/`, respectively.

After all three seeds have sealed scores, run the evidence verifier through
the scheduler with the arc interpreter.

```bash
python "$step/execution/script/paired_evidence.py" "$step/output/paired-transfer"
```

The command writes `log/matrix.json` after verifying each seed's
loss-selected source and target EMAs, reference pairing, and prediction
lineage. It recalculates the group-level metrics from the sealed
validation truth and prediction arrays.

## All-half learning rate

The `all-half` arm loads the original Tahoe bridge and uses a base
learning rate of `1e-4` in every tenv optimizer group. It keeps the
original group membership, weight decay and warmup and cosine schedule.
Use the target commands above with `--arm all-half` and
`out="$step/output/paired-transfer/all-half/seed-$seed"`.
The selected EMA and its fit seal are under
`output/checkpoints/06-paired-all-half-fold-0-seed-$seed/fit/`.

## Main-half learning rate

The `main-half` arm loads the original Tahoe bridge and fine-tunes the
Cross-DiT denoiser at a base learning rate of `1e-4`. The external
conditioning encoder, including its drug and dose components, uses
`2e-4`. Each component follows the original warmup and cosine schedule
with the original weight decay rules. The fit record lists the parameters,
base learning rate and weight decay of each group.

Use the target commands above with `--arm main-half` and
`out="$step/output/paired-transfer/main-half/seed-$seed"`.
The selected EMA and its fit seal are under
`output/checkpoints/07-paired-main-half-fold-0-seed-$seed/fit/`.

## AdaLN-only fine-tuning

The `adaLN-only` arm loads the original Tahoe bridge and trains only the
treated and control adaLN modulation weights and biases in every DiT
block and both output heads. It uses the original `2e-4` base learning
rate and native schedule. The fit record names the trainable parameters
and their optimizer groups. Frozen parameters, including the dose
embeddings, retain their bridge values in the loss-selected EMA.

Use the target commands above with `--arm adaLN-only` and
`out="$step/output/paired-transfer/adaLN-only/seed-$seed"`.
The selected EMA and fit seal are under
`output/checkpoints/09-paired-adaLN-only-fold-0-seed-$seed/fit/`.

## Reconstruct the paired comparison through the scheduler

Run the commands below from the runtime checkout. The matrix uses
`<matrix-root>/<arm>/seed-<seed>` for each cell, with `reference` as the
original Transfer arm. The target verifier finds its same-seed reference
in this layout. Use `--reference` to point a standalone target verification
at an explicit reference cell. Fit checkpoints always live in the runtime
checkout's Step 06 `output/checkpoints/`.

The queue generator pins the worker, cell and available input plan digests
for each phase. The Step 04 scheduler records attempts and telemetry in a
`scheduler/` directory alongside each queue. The commands assume a new
matrix schedule. For a fresh evidence check of a completed matrix, use a
separate directory for the new queues and leave the historical ledger in place.

```bash
set -e
step=stage2/06_pdiff_transfer_optimize
base="$PWD/$step/output/paired-transfer"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
scheduler=stage2/04_pdiff_hyperparameter/execution/script/schedule.py

run_phase() {
    local arm=$1 phase=$2 queue
    queue="$base/schedule/$arm-$phase/queue.json"
    "$pdiff" "$step/execution/script/paired_queue.py" "$base" "$arm" "$phase" --queue "$queue"
    "$pdiff" "$scheduler" run "$queue"
}

for phase in verify predict score; do
    run_phase reference "$phase"
done
run_phase log-source fit
for arm in random log all-half main-half adaLN-only; do
    for phase in verify fit predict score matrix-verify; do
        run_phase "$arm" "$phase"
    done
done
```

Each matrix verification writes `<matrix-root>/<arm>/matrix.json`. After
all five have verified, submit the [publication command](#paired-optimization-results)
as a CPU job through the same scheduler. The publisher checks the queued
phase and input digests against the current cells, then verifies each
successful attempt's log against the corresponding sealed output. For a
fresh run, use `--record` and `--figure` to write the newly generated
resource ledger and PNG to separate paths for comparison.
