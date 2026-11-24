# TranSiGen width and dose encoding

See the [historical evidence index](../03_data_train/03.1_hvg-deg-5k/historical/README.md) for sealed records.

The experiment crosses four dose encodings (no dose, latent numeric, latent
embedding, and FiLM) with four TranSiGen widths (original, 1.5×, 2×, and 2.5×).
The 48 fits reuse the [Step 04 tenv KPGT freeze](../07_tsg_dose/README.md)
and its three cold-drug-dose folds. The held concentrations are 0.2, 1, and
10 µM. The test set for each fold is the intersection of its held parents
and concentration. Both off-intersection strips stay excluded.

The architecture scales internal vector dimensions and retains the frozen
KPGT input and 5,000-gene expression interfaces. Every coordinate trains
from scratch with the same seed and selects its checkpoint by validation
treated-expression MSE. Test Δ metrics use the matched original vehicles.

The transfer experiment and its recovery artifacts are now in
[Step 06](../09_tsg_transfer/README.md).

## Results

All 48 fits were verified against the frozen inputs. Fold 0 scored 469 test
profiles at 0.2 µM, fold 1 scored 468 at 1 µM, and fold 2 scored 465 at 10 µM.
All twelve original-width fits reproduced their accepted Step 04 learning
curves, validation predictions and per-profile test scores byte for byte.

| Metric | Leading width and encoding | Equal-fold mean | Fold SD |
|---|---|---:|---:|
| Δ Pearson r | 2×, no dose | 0.27879 | 0.01057 |
| Δ Spearman ρ | Original, no dose | 0.28963 | 0.00868 |
| Δ RMSE | Original, no dose | 0.41786 | 0.00818 |
| Δ R² | Original, no dose | −0.05178 | 0.04466 |

The 2× no-dose variant led mean Δ Pearson r by about 0.001 over original-width
no dose. The three other metrics favored original-width no dose. Several
wider models had large transient training-loss spikes. The linear training
axes retain those spikes and compress the lower-loss epochs in the affected
panels. Validation MSE is plotted separately and determines checkpoint
selection. The results describe these three folds. Step 06 uses its own
model-selection protocol.

The [acceptance receipt](log/281-acceptance/acceptance.json) binds the input
freeze, 48 fit manifests, resource logs, plotted CSVs and four PNGs. The
[scientific audit](log/281-acceptance/scientific-acceptance.json) records exact
held-row coverage and the Step 04 comparison. The first capped attempt sealed
two fits before interruption. The resumed attempt reused those fits and
completed the other 46 in 45 min 41 s with peak process RSS of 1,876,744 KiB.
Its [10-second resource log](log/281-acceptance/resumed-resources-10s.txt)
contains 252 samples, with a minimum available WSL memory of 20,252,128 KiB,
no swap use, and a maximum total GPU memory use of 9,232 MiB.

## Raw-first comparison

Use the [pinned checkpoint replay](../04_tsg_axis/README.md#replay-a-sealed-runtime) for accepted pre-cutover checkpoints.

`execution/raw_first.py` prepares the versioned tenv KPGT inputs, evaluates
fixed predictions for a registered width, dose encoding and cold-drug-dose
fold, and renders the resulting summary. Preparation defaults to the
archived well-level raw-first full view. It can also reuse the verified tenv input from
the [dose comparison](../07_tsg_dose/README.md) with `prepare --tenv`.

```bash
STEP=stage2/08_tsg_width_dose
CANDIDATE=$(mktemp -d)
python "$STEP/execution/raw_first.py" prepare --output "$CANDIDATE/inputs"
python "$STEP/execution/raw_first.py" evaluate --inputs "$CANDIDATE/inputs" \
  --width original --method no_dose --fold 0 --prediction /path/to/prediction.npy \
  --prediction-ids /path/to/profile_ids.csv --output "$CANDIDATE/report"
python "$STEP/execution/raw_first.py" render --summary "$CANDIDATE/report/summary.csv" \
  --output "$CANDIDATE/metrics.png"
```

Supply a treated-expression prediction matrix in held profile and gene-axis
order. The identity CSV contains one `profile_id` column. The report records
the KPGT input, representation, corpus, axis and view digests. The renderer
uses only the summary CSV and Stage 2 plotting and palette APIs.

## Figures

- [learning_curves.png](figure/learning_curves.png) shows training loss and validation MSE from the Step 04 recipe
  for the three folds at each width and dose encoding.
- [metrics_heatmap.png](figure/metrics_heatmap.png) shows the equal-fold
  mean for each width and dose encoding in one panel per Δ metric.
- [metrics_bar.png](figure/metrics_bar.png) shows the same means as bars, with individual fold dots
  and sample standard deviation error bars.
- [metrics_per_dose_heatmap.png](figure/metrics_per_dose_heatmap.png) shows sample-equal means at each held
  concentration. Each panel row is one concentration and each panel's vertical
  axis follows the width order in the overall heatmap. Each held concentration occurs in one fold, so the
  concentration and fold effects cannot be disentangled by these panels.

Each heatmap cell displays its mean to three decimal places. Blue indicates
higher PCC, Spearman, and R² or lower RMSE. Red indicates poorer scores. White
marks the midpoint of the displayed metric range. The held-dose heatmap shares
one scale per metric across all widths and concentrations.

Reproduce the current heatmaps from the published metric CSVs:

```bash
python stage2/08_tsg_width_dose/figure/heatmaps/render_heatmaps.py --output "$(mktemp -d)"
```

The [heatmap acceptance record](log/heatmaps/acceptance.json) maps the sealed
balloon figures to their successors and records the CSV-only replay hashes.
