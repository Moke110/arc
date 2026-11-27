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

## Formal results

All 90 fits were accepted across three folds and training seeds 17, 42 and 2026.
Each checkpoint exactly replayed its validation objective, validation predictions
and held-out predictions in a fresh process. Values below are equal-coordinate
mean ± sample SD. Paired columns show the arm score minus no-dose within
the same fold and seed.

| Arm | Δ PCC ↑ | Δ RMSE ↓ | Paired Δ PCC | Paired Δ RMSE |
| --- | ---: | ---: | ---: | ---: |
| TranSiGen | 0.3640 ± 0.0166 | 0.2767 ± 0.0037 | — | — |
| no-dose | 0.3913 ± 0.0235 | 0.2472 ± 0.0046 | +0.0000 ± 0.0000 | +0.0000 ± 0.0000 |
| Xmol-scaling | 0.3421 ± 0.1135 | 0.2736 ± 0.0479 | -0.0491 ± 0.0979 | +0.0264 ± 0.0503 |
| Xmol-FiLM | 0.3342 ± 0.0618 | 0.2569 ± 0.0068 | -0.0570 ± 0.0593 | +0.0097 ± 0.0050 |
| Zmol-scaling | 0.3392 ± 0.1476 | 0.2679 ± 0.0424 | -0.0521 ± 0.1303 | +0.0208 ± 0.0448 |
| Zmol-FiLM | 0.3418 ± 0.0414 | 0.2543 ± 0.0053 | -0.0494 ± 0.0395 | +0.0071 ± 0.0053 |
| Xdose-num | 0.3190 ± 0.0708 | 0.2606 ± 0.0080 | -0.0722 ± 0.0668 | +0.0135 ± 0.0086 |
| Xdose-MLP | 0.3313 ± 0.0635 | 0.2567 ± 0.0078 | -0.0599 ± 0.0653 | +0.0096 ± 0.0074 |
| Zdose-num | 0.3512 ± 0.0544 | 0.2561 ± 0.0081 | -0.0400 ± 0.0425 | +0.0089 ± 0.0108 |
| Zdose-MLP | 0.3700 ± 0.0324 | 0.2511 ± 0.0041 | -0.0212 ± 0.0254 | +0.0039 ± 0.0056 |

No-dose had the highest mean Δ PCC and lowest mean Δ RMSE across all arms.
Zdose-MLP led the dose-input arms on both metrics. Its paired mean changes
relative to no-dose were −0.0212 in Δ PCC and +0.0039 in Δ RMSE.

Large transient objective spikes compress some learning-curve panels.
The [learning-curve CSV](figure/data/learning-curve.csv) contains the full traces.

[Fold and seed scores](log/fold-scores.csv) · [Foldwise summaries](log/foldwise-summary.csv) ·
[All metrics](log/summary.csv) · [Fit lineage](log/runs.csv) · [MMD calibration](log/bandwidths.csv)

[Learning curves](figure/learning-curve.png) · [Absolute metrics](figure/metrics.png) ·
[Paired metrics](figure/metrics-paired.png) · [CSV-only figure receipt](log/figure-render.json) ·
[Formal acceptance](log/378-formal.json)

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

## Verify and install the formal report

Run the fit command above without `--bounded` in a new scratch directory.
Keep the ten-second sampler running through the final fit measurement.
The formal publisher starts fresh scoring and independent audit processes.
The audit reconstructs held-out Δ from the frozen plate controls, checks
metric eligibility and fold–seed statistics, and compares plotted values
with the audited reports and sealed learning traces. Numeric audit tolerances
are 10⁻¹² absolute and relative. Prediction and figure replays are exact.

After the formal matrix finishes, use the same `work` and `step` paths.

```bash
formal="$step/execution/dose_formal.py"
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  python "$formal" build --bundle "$work" --publication "$work/publication" --scratch "$work/scratch"

systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  python "$formal" verify --bundle "$work" --publication "$work/publication" --scratch "$work/scratch"

systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  python "$formal" install --bundle "$work" --publication "$work/publication"
```

Pass `--step03` to `build` and `verify` when working in a linked checkout.
The install command checks the passing receipts and atomically materializes
the fits, phase audit evidence and report tables at `output/formal`.
It hardlinks the immutable fit files and copies the audit evidence and
report tables. Keep the bundle and destination on the same filesystem.
The published PNGs, plotted CSVs and compact report tables live in this
step's `figure/` and `log/`. To verify the installed study,
run `verify --scratch "$work/scratch"` with the default bundle and publication.
