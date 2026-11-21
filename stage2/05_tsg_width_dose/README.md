# TranSiGen width and dose encoding

The experiment crosses four dose encodings (no dose, latent numeric, latent
embedding, and FiLM) with four TranSiGen widths (original, 1.5×, 2×, and 2.5×).
The 48 fits reuse the [Step 04 tenv KPGT freeze](../04_tsg_dose/README.md)
and its three cold-drug-dose folds. The held concentrations are 0.2, 1, and
10 µM. The test set for each fold is the intersection of its held parents
and concentration. Both off-intersection strips stay excluded.

The architecture scales internal vector dimensions and retains the frozen
KPGT input and 5,000-gene expression interfaces. Every coordinate trains
from scratch with the same seed and selects its checkpoint by validation
treated-expression MSE. Test Δ metrics use the matched original vehicles.

The transfer experiment and its recovery artifacts are now in
[Step 06](../06_tsg_transfer/README.md).

## Reproduction

From the repository root with the `arc` environment active:

```bash
STEP=stage2/05_tsg_width_dose
INPUTS=stage2/04_tsg_dose/execution/checkpoints/tenv-kpgt-v1/inputs
SOURCE=pert_model/TranSiGen/source
WORK="$STEP/execution/checkpoints/width-dose-v1"
CAP=stage2/06_tsg_transfer/execution/capped_phase.sh

python "$STEP/execution/run.py" plan
"$CAP" "$WORK" width-dose-formal python "$STEP/execution/run.py" run \
  --inputs "$INPUTS" --source "$SOURCE" --candidate "$WORK" --device cuda
python "$STEP/execution/run.py" verify --candidate "$WORK"
python "$STEP/figure/script/prepare.py" --candidate "$WORK" --output "$WORK/report/data"
python "$STEP/figure/script/render.py" --data "$WORK/report/data" --output "$WORK/report"
```

For a bounded GPU check, use a separate candidate and pass
`--width 2.5x --method FiLM --fold 0 --smoke-rows 8 --max-epochs 1` to
`run`. The formal driver resumes completed, verified coordinates when
rerun with the same candidate identity.

## Figures

- `learning_curves.png` shows the Step 04 training loss and validation MSE
  for the three folds at each width and dose encoding.
- `metrics_balloon.png` shows one panel per Δ metric. Color and text report
  the equal-fold mean. Diameter increases with the mean's z-score against
  all 48 individual fold scores for that metric.
- `metrics_bar.png` shows the same means as bars, with individual fold dots
  and sample standard deviation error bars.
- `metrics_per_dose_ballone.png` shows sample-equal means at each held
  concentration. Each held concentration occurs in one fold, so the
  concentration and fold effects cannot be disentangled by these panels.
