# Stage 2 benchmark results

Produced by `stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/evaluate.py` (called by `run.sh`).

## Files

| file | content |
|---|---|
| `metrics.csv` | long metrics: group/drug/fold/headline, all and atc-moa subsets; `k=0` rows are regression metrics |
| `per_drug.csv` | per model/fold/drug mean over doses; AUPRC is k=200 macro |
| `leaderboard.csv` | one row per model; rank by headline k=200 macro AUPRC |
| `bootstrap.npz` | stratified 92-drug bootstrap, seed 20261022, 2,000 replicates |
| `run_cards.csv` | per model/fold wall-clock, peak GPU (total device memory), final epoch/step, best epoch/step |
| `fidelity.csv` | one row per declared deviation from every `models/<model>/FIDELITY.md` |
| `../../README.md` | the substep report: leaderboard, conclusions, limitations |

## Conventions

- Headline = mean of the available fold-level macro AUPRCs at k=200.  A
  model needs all five folds to be `winner_eligible`; a dry run with one
  fold is reported but not eligible.
- AUPRC is always direction-macro and group-macro; there is no micro
  AUPRC pooling.  `auprc_up`/`auprc_down` are audit rows only.
- Regression rows use `k=0` and `direction=macro`; R² uses the true mean
  and RMSE is in raw `log1p_cp10k` Δ units.
- A constant prediction has NaN PCC/Spearman, counted in `nan_counts` and
  labelled `constant_prediction`; it is not imputed as 0.
- `subset=atc-moa` is descriptive; it never creates a second leaderboard.
