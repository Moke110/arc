# Step 10 CarPS-VAE dose inputs

Compare nine CarPS-VAE dose inputs and a separate TranSiGen reference on the
[frozen Denoise-DEG 5k release](../03_data_train/03.7_tenv_denoise-5k/README.md).
The [experiment design](../../docs/adr/0103-carps-vae-dose-inputs.md) fixes the
conditions, model interfaces, objectives and primary metrics.

Test the intersection of held drug parents and held dose on each fold. Use
other parents at the other two doses for training and validation. The two
single-holdout strips are excluded. Each held dose is tied to a distinct
parent group, so fold-level differences conflate dose with the test parent
population.

| Fold | Held dose (µM) | Train | Validation | Test | Excluded |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.2 | 552 | 62 | 154 | 615 |
| 1 | 1 | 554 | 62 | 153 | 614 |
| 2 | 10 | 552 | 62 | 154 | 615 |

The CarPS-VAE arms use frozen UniMol and the selected shared Δ decoder and
objective. TranSiGen uses frozen KPGT, the original XX-XX architecture and
author objective, and measured same-plate DMSO without dose. Every arm is scored
on predicted treated expression minus measured matched control. Dose contrasts
pair CarPS-VAE scores with no-dose within each fold and seed. The TranSiGen
absolute scores combine architecture, representation and objective differences.

The primary held-out scores are full-axis Δ PCC and Δ RMSE. Calibrate MMD with
the corresponding fold's complete training partition on the frozen axis and
share the bandwidth across arms and seeds. Aggregate defined coordinate
values with equal weight and report sample SD. The learning curves report
each model's own objective.

## Bounded CUDA verification

[Ticket #377](https://github.com/Moke110/carps/issues/377) exercised all ten arms
on fold 0 and seed 17. The [acceptance receipt](log/377-bounded.json) records the
registered interfaces, minimum-validation-objective checkpoints, exact replay
of the validation objective and validation and test predictions, independently
derived scores, CSV-only figure reproduction and resource measurements. The
bounded fit uses the first 64 training, 24 validation and 24 test conditions
of the registered fold for two epochs.

## Reproduce the CUDA matrix and score it

Use the pinned `arc` environment and the [resource envelope](../../docs/agents/execution-environment.md).
Create separate scratch directories for bounded and formal runs. In a linked
worktree, pass the materialized `--step03`, `--representations` and `--source`
roots to the fit command and `--step03` to the report. Start the Stage 2
coordinator sampler in a separate shell before fitting. Keep it running until
its first post-fit sample is recorded. The formal invocation omits `--bounded`
and requires a matching accepted bounded receipt.

```bash
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/step10-bounded.XXXXXX")
mkdir -p "$work/scratch" "$work/matplotlib"
export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
stage=stage2_CarPS-VAE
step=$stage/10_cpsv_dose

# Separate shell with the same work path.
python "$stage/06_tsg_shape/06.1_round1/execution/sample_resources.py" \
  --output "$work/coordinator.jsonl"

systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  /usr/bin/time -v -o "$work/command.resources" \
  python "$step/execution/dose_matrix.py" --bounded --output "$work/matrix" --logs "$work/matrix-audit"

systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  /usr/bin/time -v -o "$work/report.resources" \
  python "$step/execution/dose_report.py" --bounded --matrix "$work/matrix" --output "$work/report" \
  --figure-data "$work/figure/data"

systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  /usr/bin/time -v -o "$work/figure.resources" \
  python "$step/execution/dose_figures.py" build --data "$work/figure/data" --output "$work/figure" \
  --scratch "$work/scratch" > "$work/figure-replay.json"

python "$step/execution/dose_gate.py" --bundle "$work" --output "$work/377-bounded.json"
```
