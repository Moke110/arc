# Step 06.3 TranSiGen architecture combinations

[ADR 0097](../../../docs/adr/0097-transigen-shape-round3.md) registers six
architectures on the frozen Step 06.2 denoise-5k protocol. The three tenv
cold-drug folds use measured same-plate DMSO, frozen KPGT and 32 repeats
of the numeric dose standardized on the training fold.

| Architecture | Registered configuration |
| --- | --- |
| `0_original` | Original Step 06.1 expression architecture and prediction head. |
| `1_output-delta` | Independent signed Δ decoder. |
| `2_output-delta-no-x2` | Treated decoder emits signed Δ and omits the treated reconstruction. |
| `3_z-interact_FiLM` | FiLM prediction head with identity modulation at initialization. |
| `4_1+3` | Independent signed Δ decoder and FiLM prediction head. |
| `5_2+3` | Treated decoder emits signed Δ with a FiLM prediction head. |

The four direct Δ arms add signed Δ to the measured control to predict
treated expression. Arms 2 and 5 preserve the measured treated encoder
and all three native KL terms. They use
`0.1 K + MSE_X1 + 2 MSE_Δ + 1.75 L1_Δ`. The remaining arms use the
Step 05.2 objective, `0.1 K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 L1_Δ`.

The registered training seeds are 17, 42, 2026 and 43. Training seeds
repeat fits within each held-out split. Twelve seed-42 fits from arms 0–3
reuse accepted Step 06.2 evidence, with the three originals traced to the
sealed Step 06.1 source.

## Reproduce the bounded CUDA smoke

Run from the materialized checkout in the `arc` environment using the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime).
Keep the [coordinator sampler](../06.1_round1/execution/sample_resources.py)
active at ten-second intervals.

The smoke uses all six arms on folds 0 and 1, training seeds 17 and 42,
and two epochs. Each fold supplies its first 64 training, 24 validation and
24 held-out registered conditions. Validation and test sampling seeds
remain 22617 and 100042.

```bash
set -e
STEP="$PWD/stage2/06_tsg_shape/06.3_round3"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-composition-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib" "$WORK/plot-csv"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"

capped_step04 python "$STEP/execution/matrix_composition.py" --smoke \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
python "$STEP/execution/prior_composition.py" > "$WORK/source-lineage.json"
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_composition.py" "$phase" --smoke \
    --matrix "$WORK/matrix" --output "$WORK/report" --scratch "$WORK/scratch"
done
cp "$WORK/report/"{delta-metrics,fold-scores,epoch-objectives,primary-contrasts,primary-contrast-summary}.csv "$WORK/plot-csv"
for output in figures figures-replay; do
  python "$STEP/figure/script/render_composition.py" --smoke \
    --data "$WORK/plot-csv" --output "$WORK/$output"
done
diff -r "$WORK/figures" "$WORK/figures-replay"
```
