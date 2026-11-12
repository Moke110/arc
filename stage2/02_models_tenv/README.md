# stage2/02_models_tenv — the tenv model comparison

The active Stage 2 model-comparison step (ADR 0063). Its entrant is PRnet,
trained from random initialization on the indexed `tenv` training corpus with the
shared Morgan-2048 chemical adapter and the author's Gaussian NLL objective; its
deterministic reference is `Mean`. The step owns the trained entrants, their
predictions and the report figures; the shared metrics, palette, fold/seed
identities, evaluation driver and figure APIs stay in
[`../01_bases/`](../01_bases/README.md) and are never re-implemented here.

## Status (branch #165): accepted 2 folds × 2 seeds matrix

The declared matrix (ADR 0063) is 2 outer cold-drug folds × 2 training seeds,
giving four PRnet run samples per entry, plus the deterministic `Mean`
reference. All four PRnet runs and all four `Mean` fold copies reached terminal
status `complete` and were scored; both figures were rebuilt from their tracked
plot-value tables. The acceptance record is
[`log/165-acceptance.json`](log/165-acceptance.json): environment, declared
seeds, per-run identity hashes, per-run terminal status, wall time, peak RSS,
metric tables, file hashes and tolerances.

Runs are serialized, one heavy job at a time, each through
`execution/capped.sh` (10 GiB) under `/usr/bin/time -v`, 2 threads. This boot
had **no GPU** (`torch.cuda.is_available()` false), so every run used the
documented ADR 0063 CPU fallback.

| entry | fold | seed | epochs (best) | stop | wall s | peak RSS | PCC | Spearman | RMSE | R² |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| prnet | 0 | 17 | 49 (44) | no strict decrease | 87.0 | 1982 MiB | 0.1917 | 0.2864 | 0.2299 | −0.2940 |
| prnet | 0 | 42 | 36 (31) | no strict decrease | 66.1 | 2024 MiB | 0.1371 | 0.2715 | 0.2609 | −0.8487 |
| prnet | 1 | 17 | 57 (52) | no strict decrease | 104.1 | 2027 MiB | 0.2144 | 0.2729 | 0.2205 | −0.1733 |
| prnet | 1 | 42 | 41 (36) | no strict decrease | 75.0 | 1995 MiB | 0.1659 | 0.2755 | 0.2491 | −0.6674 |
| mean | 0 | 17/42 | — | deterministic | 3.0–6.6 | ≤1074 MiB | 0.1934 | 0.1717 | 0.2070 | 0.0321 |
| mean | 1 | 17/42 | — | deterministic | 3.1 | 1040 MiB | 0.1961 | 0.1754 | 0.2067 | 0.0308 |

Indirect metrics (equal-weight macro over defined test groups; 459 per fold):

| entry | fold | max_k@P80S40 | Recall@max_k |
|---|---:|---:|---:|
| prnet | 0 | undefined (both seeds) | undefined (both seeds) |
| prnet | 1 | 27.5 / 34.7 | 0.0190 / 0.0195 |
| mean | 0 | 13 | 0.1019 |
| mean | 1 | 31 | 0.0278 |

Reading: PRnet beats `Mean` on per-well `Spearman` in every run (≈0.27–0.29 vs
≈0.17) but is behind on `PCC`, `RMSE`, `R²` and the group-ranking family, where
`Mean` has a defined `max_k@P80S40` in both folds and PRnet has none in fold 0.
`Precision@k` and `Spearman@k` curves are in the figure data and in
`output/evaluation/per_group.csv`.

## Figures

| figure | PNG | complete plotted values | renderer |
|---|---|---|---|
| eight-panel metrics (2 × 4) | [`figure/metrics.png`](figure/metrics.png) | [`figure/data/metrics_values.csv`](figure/data/metrics_values.csv) | [`figure/script/render_figures.py`](figure/script/render_figures.py) |
| PRnet learning curve | [`figure/learning_curve.png`](figure/learning_curve.png) | [`figure/data/learning_curve_values.csv`](figure/data/learning_curve_values.csv) | same script |

`metrics.png`: direct panels (`PCC`, `Spearman`, `RMSE`, `R²`) and threshold
panels (`max_k@P80S40`, `Recall@max_k`) are bars of the four-run mean with one
dot per run; `Precision@k` and `Spearman@k` are `k`-curves with the four-run mean
line and min–max band. PRnet is the palette blue family, `Mean` the orange entry
slot; an entry's hue colours its line/edge/dots and the lighter shade fills its
band/bar. `learning_curve.png`: one PRnet subplot, epoch on the horizontal axis,
blue training and green validation, one shade per `(fold, seed)` run.

Both PNGs regenerate from the tracked CSVs alone, with the step's `output/`
unreadable; the committed replay test proves the drawn values equal the recorded
values byte for byte.

## Replay and rebuild

```bash
# Rebuild both PNGs from the tracked plot values alone (never opens output/):
python figure/script/render_figures.py --figure-data figure/data --out-dir figure

# Rebuild the whole accepted result from the canonical checkout:
# inputs -> four serialized PRnet runs -> Mean -> score -> figures -> record
python execution/run_matrix.py --all          # add --resume to reuse complete runs
python execution/run_matrix.py --plan-only    # just print the declared matrix
```

`run_matrix.py` writes a candidate tree under `execution/checkpoints/matrix-165/`
and touches no stable path until every run is `complete` and the whole set is
scored (fail closed). It publishes `output/{inputs,runs,evaluation,figure_data}/`,
copies the two plot-value tables into `figure/data/`, renders both PNGs and
re-verifies them with `output/` moved away.

## Products and paths

| path | content |
|---|---|
| `output/inputs/{fold_X}/{seed_Y}/` | complete predictor-input release per matrix cell and its manifest |
| `output/runs/prnet/{fold_X}/{seed_Y}/` | `prediction.npy`, `case_ids.json`, `best.pt`, `config.json`, `checkpoint.json`, `identity.json`, `terminal.json`, `progress.json`, `resources.json` |
| `output/runs/mean/{fold_X}/{seed_Y}/` | deterministic `Mean` run artifacts (same shape, empty training history) |
| `output/evaluation/` | `per_well.csv`, `per_group.csv`, `per_run.csv`, `evaluation.json` from the shared #164 driver |
| `output/figure_data/` | candidate copies of the two plot-value tables (identical to `figure/data/`) |
| `figure/data/*.csv`, `figure/*.png` | tracked plot values and report figures |
| `execution/run_matrix.py` | the one-command candidate → acceptance → publish/replay |
| `execution/script/mean_reference.py` | the deterministic reference (retired `delta_mean` arithmetic) |
| `log/165-acceptance.json`, `log/resources/{prnet,mean}-fold_*-seed_*.time` | acceptance record and raw `/usr/bin/time -v` measurements |
| `log/163-bounded-run.json`, `log/164-evaluation-driver.json` | earlier admission and driver records |

## What the entrant is

- Inputs: baseline expression (same-plate DMSO mean, `layers/control`), binary
  Morgan-2048 fingerprint of the corpus's normalized parent structure, and the
  explicit dose in µM.
- Seam: a jointly trained `Linear(2048, 1024)` Morgan adapter multiplied by
  `log10(dose_um + 1)` feeds the unchanged upstream `models.PRnet.PGM`;
  `pert_model/PRnet/source` is imported, never edited.
- Objective: Gaussian NLL on the expression target (`delta + baseline`).
- Training wells only drive the fit; validation wells drive epoch-based early
  stopping (5 consecutive epochs without a strict decrease, independently 20
  above the running minimum, at most 100 epochs, keep the lowest validation
  loss). In all four runs the 5-epoch rule stopped the run, as ADR 0063 predicts.
- The selected lowest-validation-loss checkpoint is reloaded from disk before
  prediction (recorded reload error 0 in the run identity evidence).

## The `Mean` reference

`Mean` implements the retired `delta_mean` semantics exactly: **one global
prediction vector**, the equal-weight mean over the training-partition drug–dose
groups of each group's mean measured Δ, copied to every test well. It does not
look up the test well's own group. Cold-drug test groups are absent from the
training partition by construction and receive the same global vector; there is
no fallback. The two seed rows for a fold share identical values, so `Mean`'s
four entries are four fold copies, not independent training results.

## Limitations

- Four run samples (2 folds × 2 seeds) are a **spread indicator, not a
  confidence interval** and not independent biological replicates; the figure
  shows their min–max band, not an uncertainty estimate.
- The truth gate (`padj < 0.05`, `|log2FC| >= 2`) leaves ≈65.5 gated genes per
  group on average; 283 of 1,377 groups have none. Per fold, 459 test groups have
  defined gates, and `max_k@P80S40` is bounded by each group's own gated count
  (a group with fewer than 13 gated genes cannot reach `Precision@10 = 0.8`).
  `Spearman@k` is undefined for `k < 3`; a PRnet fold-0 run with no qualifying
  `k` is reported as undefined, never as zero.
- Direct metrics are single-well statistics over the 18,164-gene axis and carry
  TempO-Seq well-level noise; they measure per-sample fidelity, not group-level
  accuracy.
- CPU-only acceptance: no GPU was available this boot, so the ADR 0063 CPU
  fallback is the measured path (wall 66–104 s per PRnet run, ≤2027 MiB peak
  RSS); GPU timing is untested.
- The entrant is PRnet alone; `no_change` is withdrawn and no second entrant is
  promised.

## Historical boundary

The superseded comparison implementation, its figures and records are not in the
working tree; they remain untouched and replayable from Git history
(`git show 1543781^:stage2/03_modes_tenv/...`). This step adds no dependency on
them and re-implements none of their code.
