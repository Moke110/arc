# Step 06.3 TranSiGen architecture combinations

The numbered commands use the [shared TranSiGen training and evaluation APIs](../../util/README.md#transigen-training-and-evaluation).

[ADR 0097](../../../docs/adr/0097-transigen-shape-round3.md) registers six
architectures on the frozen Step 06.2 denoise-5k protocol. The three tenv
cold-drug folds use measured same-plate DMSO, frozen KPGT and 32 repeats
of the numeric dose standardized on the training fold.

| Architecture | Registered configuration |
| --- | --- |
| `0_original` | Original Step 06.1 expression architecture and prediction head. |
| `1_output-delta` | Independent signed Δ decoder. |
| `2_output-delta-no-x2` | Treated decoder emits signed Δ and omits the treated reconstruction. |
| `3_z-interact_FiLM` | FiLM with LayerNorm before the Gaussian projection. |
| `4_1+3` | Independent signed Δ decoder and normalized FiLM head. |
| `5_2+3` | Treated decoder emits signed Δ with a normalized FiLM head. |

The four direct Δ arms add signed Δ to the measured control to predict
treated expression. Arms 2 and 5 preserve the measured treated encoder
and all three native KL terms. They use
`0.1 K + MSE_X1 + 2 MSE_Δ + 1.75 L1_Δ`. The remaining arms use the
Step 05.2 objective, `0.1 K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 L1_Δ`.

Arms 3, 4 and 5 apply `LayerNorm(100, elementwise_affine=False)` immediately
after FiLM modulation and before the Gaussian projection. The conditioner
starts with identity modulation. The normalization has no trainable parameters.
The [first formal attempt](log/356-blocked.json) recorded eight overflows in
the prediction KL during training and 52 exact independent replays.

The registered training seeds are 17, 42, 2026 and 43. They repeat fits
within each of the three held-out splits. The revised 72-coordinate matrix
trains all 36 fits of arms 3–5 from scratch. Arms 0–2 reuse nine accepted
seed-42 fits from Steps 06.1 and 06.2 and 27 independently replayed fits from
the first Round 3 attempt. The source manifests and registered populations
are checked before rescoring.

## Formal results

On the 12 paired fold and training-seed coordinates, `4_1+3` increased
mean Δ PCC by 0.004550 and reduced mean Δ RMSE by 0.009452 relative to
`0_original`. Relative to its signed Δ readout alone, its mean Δ PCC was
0.000485 lower and its mean Δ RMSE was 0.000973 higher. `5_2+3` reduced
mean Δ RMSE relative to the original and had lower mean Δ PCC. Both
combinations increased mean Δ PCC and reduced mean Δ RMSE relative to
normalized FiLM alone.

Each score below gives equal weight to the 12 fold and training-seed
coordinates. All 12 scores and differences were defined for each metric.

| Architecture | Δ PCC | Δ RMSE | JAC_top200 | JAC_pos100 | JAC_neg100 | Δ Spearman | Δ E-distance | Δ MMD | Δ R² |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `0_original` | 0.424154 | 0.251615 | 0.149628 | 0.126681 | 0.180589 | 0.407390 | 0.026717 | 0.158379 | -0.054433 |
| `1_output-delta` | 0.429189 | 0.241190 | 0.162041 | 0.145929 | 0.185989 | 0.406710 | 0.026999 | 0.139171 | 0.095294 |
| `2_output-delta-no-x2` | 0.421896 | 0.241106 | 0.160724 | 0.145181 | 0.184052 | 0.398976 | 0.029478 | 0.141678 | 0.092696 |
| `3_z-interact_FiLM` | 0.396089 | 0.259210 | 0.134883 | 0.123063 | 0.153621 | 0.383520 | 0.022962 | 0.148617 | -0.139372 |
| `4_1+3` | 0.428704 | 0.242163 | 0.163303 | 0.150076 | 0.183734 | 0.406540 | 0.024876 | 0.135851 | 0.084055 |
| `5_2+3` | 0.400866 | 0.243244 | 0.150913 | 0.138536 | 0.170243 | 0.375642 | 0.027724 | 0.137394 | 0.079498 |

The signed differences below subtract the original fit on the same fold
and training seed. Values are mean ± sample SD.

| Architecture | Paired Δ PCC difference | Paired Δ RMSE difference |
| --- | --- | --- |
| `0_original` | +0.000000 ± 0.000000 | +0.000000 ± 0.000000 |
| `1_output-delta` | +0.005035 ± 0.016649 | -0.010425 ± 0.004631 |
| `2_output-delta-no-x2` | -0.002258 ± 0.018312 | -0.010509 ± 0.004267 |
| `3_z-interact_FiLM` | -0.028065 ± 0.017481 | +0.007595 ± 0.005446 |
| `4_1+3` | +0.004550 ± 0.013717 | -0.009452 ± 0.004250 |
| `5_2+3` | -0.023288 ± 0.016215 | -0.008371 ± 0.004231 |

The component and interaction contrasts use the same paired coordinates.
Higher Δ PCC and lower Δ RMSE are favorable. The last two rows describe
departures from additive differences. The four seeds repeat training within
the same three held-out folds.

| Contrast | Δ PCC difference | Δ RMSE difference |
| --- | --- | --- |
| 4−1 | -0.000485 ± 0.019832 | +0.000973 ± 0.002256 |
| 4−3 | +0.032615 ± 0.019659 | -0.017047 ± 0.003908 |
| 5−2 | -0.021030 ± 0.021985 | +0.002138 ± 0.002250 |
| 5−3 | +0.004777 ± 0.018678 | -0.015966 ± 0.003005 |
| 4−1−3+0 | +0.027580 ± 0.032154 | -0.006622 ± 0.006550 |
| 5−2−3+0 | +0.007035 ± 0.030273 | -0.005457 ± 0.005543 |

[Held-out scores](figure/metrics.png) ·
[Paired differences](figure/metrics-paired.png) ·
[Original primary contrast figure](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.3_round3/figure/primary-contrasts.png) ·
[Learning curves](figure/learning-curve.png) ·
[Fold and seed scores](figure/data/fold-scores.csv) ·
[Aggregates and sample SD](figure/data/delta-metrics.csv) ·
[Coordinate contrasts](figure/data/primary-contrasts.csv) ·
[Contrast summaries](figure/data/primary-contrast-summary.csv) ·
[Epoch objectives](figure/data/epoch-objectives.csv)

All 36 normalized fits passed exact independent validation-objective,
validation-prediction, and held-out prediction replay. They include all
eight coordinates that failed in the first attempt. The 36 reused fits
passed source and artifact verification. The rescoring reproduced all
81 earlier seed-42 metric values and eligibility records for arms 0–2.
The [formal receipt](log/356-formal.json) records the independent rescore
and the [fit ledger](log/runs.csv) records the source of each fit. The
[table audit](log/356-formal-table-audit.json) independently recomputed
the paired differences, means, sample SDs, and component contrasts.
The [matrix samples](log/356-resources.jsonl) and
[coordinator samples](log/356-coordinator-resources.jsonl) record resource use.

## Reproduce the bounded CUDA smoke

Run from the materialized checkout in the `arc` environment using the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime).
Keep the [coordinator sampler](../06.1_round1/execution/sample_resources.py)
active at ten-second intervals.

The smoke uses all six arms on folds 0 and 1, training seeds 17 and 42,
and two epochs. Each fold supplies its first 64 training, 24 validation and
24 held-out registered conditions. Validation and test sampling seeds
remain 22617 and 100042.

The [original smoke receipt](log/355-smoke.json) records exact independent replays
for all 24 fits, 216 held-out scores, 48 coordinate contrasts and
12 contrast summaries. The [objective audit](log/355-objective-audit.json)
recomputed 96 objective values and verified the seed behavior. The
[figure receipt](log/355-figure-render.json) records two matching
CSV-only renders of four smoke figures. The resource records cover
the [matrix scope](log/355-resources.jsonl) and
[coordinator samples](log/355-coordinator-resources.jsonl).

Verify the original sealed smoke using `report_composition.py verify --smoke`
at revision `7ce9e6259c6554072a7193b2830ce9fd03c54248` with the original
matrix, matrix audit, report and publication paths. That revision contains
the verifier files whose hashes are recorded in the receipt. The receipt
records `29c7e1680e79b32d40de18950260ff57ad5615bd` as the matrix execution
revision. The commands below run a fresh smoke using the current heads.

```bash
set -e
STEP="$PWD/cpst/06_tsg_shape/06.3_round3"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-composition-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib" "$WORK/plot-csv"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"

capped_step04 python "$STEP/execution/matrix_composition.py" --smoke \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
python "$STEP/execution/audit_smoke_objectives.py" --matrix "$WORK/matrix"
python "$STEP/execution/prior_composition.py" > "$WORK/source-lineage.json"
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_composition.py" "$phase" --smoke \
    --matrix "$WORK/matrix" --output "$WORK/report" --scratch "$WORK/scratch"
done
cp "$WORK/report/"{delta-metrics,fold-scores,epoch-objectives,primary-contrasts,primary-contrast-summary}.csv "$WORK/plot-csv"
for output in figures figures-replay; do
  python "$STEP/execution/render_composition.py" --smoke \
    --data "$WORK/plot-csv" --output "$WORK/$output"
done
diff -r "$WORK/figures" "$WORK/figures-replay"
```

## Revised CUDA gate

The [revised bounded study](log/356-normalized-smoke.json) passed all 24
exact fit replays and independent rescoring. The
[objective audit](log/356-normalized-smoke-objective-audit.json) checked
96 weighted objective values and the
[figure receipt](log/356-normalized-smoke-figure-render.json) records
two matching CSV-only renders of four figures.

The first revised formal scope started at 15:34:33 Asia/Shanghai on
2026-11-26, before the revised bounded check. It paused at a sealed
coordinate boundary after 13 fits. The bounded scope started at 15:42:41
and its report passed at 15:48:43. Formal processing resumed at 15:48:47
after checking the 13 sealed registrations and replay manifests against
their pre-gate hashes. Both formal scopes used identical matrix, runner and
model source files.

## Reproduce the revised formal run

The [original release reproduction](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.3_round3/README.md#reproduce-the-revised-formal-run) pins the renderer and publication verifier for the historical results. Use that revision for the formal release verification below. The current comparison figures have their own CSV-only replay.

Complete the bounded CUDA check above first. Run from the materialized checkout in the `arc` environment. Use the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
and keep the coordinator sampler active for the entire matrix, including
at least one sample before launch and one after completion. Keep the
candidate and report scratch on the same filesystem.

```bash
set -e
STEP="$PWD/cpst/06_tsg_shape/06.3_round3"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-composition-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
REVISION=$(git rev-parse HEAD)
SAMPLES="$WORK/coordinator-resources.jsonl"

python "$STEP/../06.1_round1/execution/sample_resources.py" --output "$SAMPLES" &
SAMPLER_PID=$!
stop_sampler() {
  kill -TERM "$SAMPLER_PID" 2>/dev/null || true
  wait "$SAMPLER_PID"
}
trap 'stop_sampler || true' EXIT
while [ ! -s "$SAMPLES" ]; do
  kill -0 "$SAMPLER_PID" 2>/dev/null || exit 1
  sleep 1
done
capped_step04 python "$STEP/execution/matrix_composition.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
sleep 10
stop_sampler
trap - EXIT

capped_step04 python "$STEP/execution/report_composition.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --coordinator-samples "$SAMPLES" --code-revision "$REVISION"
```

The matrix driver seals separate-process replay of the validation objective,
validation predictions and held-out predictions. The report checks the 36
retained source fits and independently rescores and renders the full study.
If verified source products have moved, supply `--round-one-matrix`,
`--round-two-matrix`, or `--prior-unnormalized` to the report command.

To verify an installed release, run `report_composition.py verify` with
`output/formal/{matrix,matrix-audit,report}` and the step directory as
`--publication`. The installed verification reads the coordinator samples
and execution revision from the publication. Use fresh `--activity` and
`--scratch` directories below the checkout's `tmp/`.

## Replay the published figures

The three current figures use the [epoch traces](figure/data/learning-curve.csv),
[held-out scores](figure/data/metrics.csv) and [paired differences](figure/data/metrics-paired.csv).
Differences subtract `0_original` within each fold and seed. The
[replay receipt](log/figure-render.json) records the accepted fit evidence,
plotted CSV and PNG digests, and renderer versions.

```bash
STEP="$PWD/cpst/06_tsg_shape/06.3_round3"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-figures.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
python "$STEP/figure/script/render_comparison.py" \
  --data "$STEP/figure/data" --output "$WORK/figures"
for figure in learning-curve.png metrics.png metrics-paired.png; do
  cmp "$STEP/figure/$figure" "$WORK/figures/$figure"
done
```
