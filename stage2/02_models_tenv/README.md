# stage2/02_models_tenv — the tenv model comparison

The active Stage 2 model-comparison step (ADR 0063). Its first and only
entrant is PRnet, trained from random initialization on the indexed `tenv`
training corpus with the shared Morgan-2048 chemical adapter and the author's
Gaussian NLL objective; the deterministic `Mean` reference follows in a later
ticket. The step owns the trained entrants, their predictions and the report
figures; the shared metrics, palette, fold/seed identities and figure APIs stay
in [`../01_bases/`](../01_bases/README.md).

## Status (this branch, #163)

The bounded real-data admission run is recorded in
[`log/163-bounded-run.json`](log/163-bounded-run.json): fold row 0 / seed row 0
(fold 0, seed 17), 3 epochs, 512 training / 128 validation / 64 test wells, a
reloaded checkpoint, and a finite 64 × 18,164 prediction. The measured peak was
about 1.3 GiB RSS on CPU; the reasoning and the no-GPU deviation are recorded in
the same file and in [`log/resources/`](log/resources/). This qualifies
execution, not convergence or formal comparison performance. The full
2 folds × 2 seeds matrix belongs to #165; `metrics.png` and
`learning_curve.png` belong to #164.

## What the entrant is

- Inputs: baseline expression (same-plate DMSO mean, `layers/control`),
  binary Morgan-2048 fingerprint of the corpus's normalized parent structure,
  and the explicit dose in µM.
- Seam: a jointly trained `Linear(2048, 1024)` Morgan adapter multiplied by
  `log10(dose_um + 1)` feeds the unchanged upstream `models.PRnet.PGM`;
  `pert_model/PRnet/source` is imported, never edited.
- Objective: Gaussian NLL on the expression target (`delta + baseline`) with
  the native positive mean and `softplus` variance.
- Training wells only drive the fit; validation wells drive epoch-based early
  stopping. The selected lowest-validation-loss checkpoint is reloaded from disk
  before prediction, so every published value comes from an artifact.
- Early stopping (ADR 0063): 5 consecutive epochs without a strict decrease,
  independently 20 consecutive epochs above the running minimum, at most 100
  epochs; the 5-epoch rule dominates and that dominance is recorded.

## Products and paths

| path | content |
|---|---|
| `output/bounded/inputs/` | bounded predictor-input release and its manifest (step-owned, untracked payload) |
| `output/bounded/runs/prnet/fold_0/seed_17/` | bounded run: `prediction.npy`, `case_ids.json`, `best.pt`, `config.json`, `checkpoint.json`, `identity.json`, `terminal.json`, `progress.json` |
| `execution/script/` | input build, fit/predict, early stopping and artifact helpers |
| `execution/models/prnet/` | the fixed author recipe (`config.json`) and the bounded-epoch override (`bounded_run.json`) |
| `execution/tests/` | the seam tests described in [`AGENTS.md`](AGENTS.md) |
| `log/` | environment locks, bounded-run record, raw resource measurements |

## Limitations

- The bounded run is an execution/artifact admission, not a performance result:
  three epochs on a well subset cannot establish convergence.
- CPU is a documented ADR 0063 fallback; this boot had no GPU
  (`torch.cuda.is_available()` false, no `/dev/dxg`), so the bounded run is a
  CPU-only deviation and the GPU path is untested in this ticket.
- The 64 test wells are a deterministic bounded prefix, not the fold's full
  test set; no metric is reported by this step yet.
- Four run samples (once #165 runs them) are a spread indicator, not a
  confidence interval and not independent biological replicates.
