# stage2/02_models_tenv — the tenv model comparison

The active Stage 2 model-comparison step (ADR 0064). Its entrant is PRnet,
trained from random initialization on the indexed `tenv` training corpus with the
shared Morgan-2048 chemical adapter and the author's Gaussian NLL objective; its
deterministic reference is `Mean`. The step owns the trained entrants, their
predictions and the report figures; the shared metrics, palette, fold/seed
identities, evaluation driver and figure APIs stay in
[`../01_bases/`](../01_bases/README.md) and are never re-implemented here.

## Status: accepted 2 folds × 2 seeds matrix

The declared matrix is 2 outer cold-drug folds × 2 training seeds, giving four
PRnet run samples per entry, plus the deterministic `Mean` reference. All four
PRnet runs and all four `Mean` fold copies reached terminal status `complete`.
The accepted results are recorded in
[`log/166-acceptance.json`](log/166-acceptance.json), the superseding
`stage2-tenv-acceptance-v2` record: environment, declared seeds and folds,
per-run identity hashes and terminal status, per-run `M`, `d*` and
`Reliability@d*`, the null/changed split, group counts, resource measurements,
every table and figure hash, the replay proof and the tolerances. The earlier
[`log/165-acceptance.json`](log/165-acceptance.json) remains as history and
carries a `superseded_by` pointer.

Runs were serialized, one heavy job at a time, each through
`execution/capped.sh` (10 GiB) under `/usr/bin/time -v`, 2 threads. This boot
had **no GPU** (`nvidia-smi` exits 255, "GPU access blocked by the operating
system"), so every run used the documented CPU fallback.

## Label rule and evaluated population

A **measured effect** is a gene with `padj < 0.05`; its direction is
`sign(log2FC)`. The earlier magnitude gate is gone: it removed 10 non-empty
groups from an already small label set, so the non-empty population rises from
**1,094 to 1,104**. The corpus holds **1,377 parent-dose groups**, of which
**273 are genuine null results** (no measured effect, 19.8% of the corpus). A
null group rewards a model that emits no signature and penalises one that
hallucinates a table; no group is dropped for an undefined metric.

The declared 2-fold matrix does not score the whole corpus. It evaluates the
**918 test groups** of fold 0 and fold 1 (459 each); **169 are null results**
(90 in fold 0, 79 in fold 1). Fold 2 has no runs and re-inference is out of
scope, so the accepted record documents the 918/169 population as a deliberate
deviation from the corpus-wide 1,377/273.

## The `@d` reliability family

For a run with largest test-group predicted magnitude `M`, the scan is 100
points `d = k · 0.01 · M` (`k = 1..100`, plotted as `d' = k/100`). A group's
**threshold window** is `W(d) = {i : |Δ̂_mean| ≥ d}`: the model's own changed-gene
table at that threshold. `H(d)` is the window genes that are measured effects
with matching direction. `Precision@d = |H|/|W|`, `Recall@d = |H|/|N|`,
`F1@d` their harmonic mean, and `Spearman@d` the predicted-versus-measured
magnitude rank correlation inside `H(d)`. `Reliability@d = F1@d · mult`, with
`mult = (1 + Spearman@d)/2` for three or more hits and `1` otherwise. The four
empty/non-empty grid cells take declared values, so every test group enters every
macro.

A run's `d*` is the argmax of its macro `Reliability@d`; its **delivered table**
is `W(d*)`. No signature table is published as a stable product — `d*` and
`Reliability@d*` are the recorded values a consumer applies to that run's
checkpoint.

## Accepted per-run results

`M`, `d*` and `Reliability@d*` are macro values over the fold's 459 test groups.
`null` / `changed` are the group counts behind the split. The two `Mean` seed
rows per fold are identical copies of the deterministic reference.

| entry | fold | seed | M | k* | d* | Reliability@d* | null groups | changed groups | null Reliability | changed Reliability |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PRnet | 0 | 17 | 5.9263 | 100 | 5.9263 | 0.19608 | 90 | 369 | 1.00000 | 0.00000 |
| PRnet | 0 | 42 | 5.9263 | 100 | 5.9263 | 0.19608 | 90 | 369 | 1.00000 | 0.00000 |
| PRnet | 1 | 17 | 5.9263 | 100 | 5.9263 | 0.16558 | 79 | 380 | 0.96203 | 0.00000 |
| PRnet | 1 | 42 | 5.9263 | 100 | 5.9263 | 0.16558 | 79 | 380 | 0.96203 | 0.00000 |
| Mean | 0 | 17/42 | 0.42570 | 32 | 0.13622 | 0.00705 | 90 | 369 | 0.00000 | 0.00877 |
| Mean | 1 | 17/42 | 0.45769 | 31 | 0.14188 | 0.00883 | 79 | 380 | 0.00000 | 0.01067 |

Direct per-well metrics (equal-weight mean over the fold's test wells) are
unchanged by the rescoring; they measure per-sample fidelity, not group-level
accuracy.

| entry | fold | seed | PCC | Spearman | RMSE | R² |
|---|---:|---:|---:|---:|---:|---:|
| PRnet | 0 | 17 | 0.1917 | 0.2864 | 0.2299 | −0.2940 |
| PRnet | 0 | 42 | 0.1371 | 0.2715 | 0.2609 | −0.8487 |
| PRnet | 1 | 17 | 0.2144 | 0.2729 | 0.2205 | −0.1733 |
| PRnet | 1 | 42 | 0.1659 | 0.2755 | 0.2491 | −0.6674 |
| Mean | 0 | 17/42 | 0.1934 | 0.1717 | 0.2070 | 0.0321 |
| Mean | 1 | 17/42 | 0.1961 | 0.1754 | 0.2067 | 0.0308 |

## Reading the result

PRnet's `Reliability@d` peaks at the **strictest grid point** `k = 100`
(`d' = 1`) in all four runs. At that edge the mean window size is **0.03 genes**
— the peak is a nearly-empty-table solution rather than a recovered gene set.
The value it reports is the share of the fold's test groups that emit an empty
window on a null group: fold 0 is exactly `90/459 = 0.196078`; fold 1 is
`76/459 = 0.165577`, because 76 of its 79 null groups are empty at `k = 100`
while three still emit one gene. PRnet's changed-group Reliability at `d*` is
zero in every run.

The consequence is load-bearing for downstream work: the **delivered table
`W(d*)` is near-empty**, and different runs deliver different table sizes. A
consumer of a PRnet checkpoint must treat table size as part of the intervention,
not as noise, and must not read `Reliability@d*` as gene-set recovery. `Mean`
peaks at small `d'` (31–32), where its window recovers a small part of the
measured set, but its macro Reliability is two orders of magnitude smaller.

Peak selection is **in-sample and optimistic**: the same group data chooses `d*`
and reports the value. The full curve ships beside the number so the optimism is
visible, and the four per-run `d*` values are the stability evidence. PRnet's
four identical `M = 5.926281` values are a genuine data property of the shared
cold-drug test wells, not a copied record.

## Figures

| figure | layout | PNG | complete plotted values | renderer |
|---|---|---|---|---|
| direct metrics | 1 × 4 | [`figure/metrics.png`](figure/metrics.png) | [`figure/data/metrics_values.csv`](figure/data/metrics_values.csv) | [`figure/script/render_figures.py`](figure/script/render_figures.py) |
| reliability threshold scan | 2 × 4 | [`figure/metrics_indirect.png`](figure/metrics_indirect.png) | [`figure/data/metrics_indirect_values.csv`](figure/data/metrics_indirect_values.csv) | same script |
| PRnet learning curve | 1 subplot | [`figure/learning_curve.png`](figure/learning_curve.png) | [`figure/data/learning_curve_values.csv`](figure/data/learning_curve_values.csv) | same script |

`metrics.png` has four bar panels (`PCC`, `Spearman`, `RMSE`, `R²`) of the
four-run mean with one dot per run. `metrics_indirect.png` has four curves over
`d'` (`Precision@d`, `Recall@d`, `Spearman@d`, `Reliability@d`) with the
four-run mean line, min–max band and `d*` marked; the Reliability panel also
carries `F1@d` as a thin reference line. Its second row has the
`Reliability@d*` and `d*` bars (four-run mean, one dot per run), a macro
three-class confusion heat map (`up` / `not significant` / `down`) at `d*`, and
the null-versus-changed Reliability split. `learning_curve.png` is one PRnet
subplot with blue training and green validation, one shade per `(fold, seed)`.

PRnet is the palette blue family and `Mean` the orange entry slot; an entry's
hue colours its line/edge/dots and the lighter shade fills its band/bar. The
reliability figure's confusion panel uses the shared sequential heat-map ramp.

All three PNGs regenerate from the tracked CSVs alone, with the step's
`output/` unreadable. The committed replay test proves each PNG rebuilds
byte-identically from `figure/data`.

## Replay and rebuild

```bash
# No-training acceptance/replay of the published matrix (never trains or re-infers):
python execution/run_matrix.py --accept-only \
  --python /home/chang/miniconda3/envs/arc/bin/python

# Rebuild the three PNGs from the tracked plot values alone (never opens output/):
python figure/script/render_figures.py --figure-data figure/data --out-dir figure

# Rebuild the whole accepted result from the canonical checkout (trains; long):
python execution/run_matrix.py --all          # add --resume to reuse complete runs
python execution/run_matrix.py --plan-only    # just print the declared matrix
```

`--accept-only` verifies the eight published runs and their identity digests,
rescores them twice with the shared driver (BLAS threads pinned to 2), compares
the candidate and published trees byte for byte, re-emits
`output/evaluation/` and `output/figure_data/`, and refuses to rewrite any
accepted figure or plot-value byte. It writes `log/166-acceptance.json`; a
failure leaves the accepted products untouched. The driver itself is
deterministic: the two rescoring passes must match byte for byte, and no
stochastic term enters the rescoring path.

`run_matrix.py --all` writes a candidate tree under
`execution/checkpoints/matrix-165/` and touches no stable path until every run
is `complete` and the whole set is scored (fail closed). It publishes
`output/{inputs,runs,evaluation,figure_data}/`, copies the three plot-value
tables into `figure/data/`, renders the three PNGs and re-verifies them.

## Products and paths

| path | content |
|---|---|
| `output/inputs/{fold_X}/{seed_Y}/` | complete predictor-input release per matrix cell and its manifest |
| `output/runs/prnet/{fold_X}/{seed_Y}/` | `prediction.npy`, `case_ids.json`, `best.pt`, `config.json`, `checkpoint.json`, `identity.json`, `terminal.json`, `progress.json`, `resources.json` |
| `output/runs/mean/{fold_X}/{seed_Y}/` | deterministic `Mean` run artifacts (same shape, empty training history) |
| `output/evaluation/` | `per_well.csv`, `per_group.csv`, `per_run.csv`, `per_run_curve.csv`, `confusion.csv` and `evaluation.json` from the shared driver |
| `output/figure_data/` | copies of the three plot-value tables (byte-identical to `figure/data/`) |
| `figure/data/*.csv`, `figure/*.png` | tracked plot values and report figures |
| `execution/run_matrix.py` | the one-command candidate → acceptance → publish, plus the no-training `--accept-only` replay |
| `execution/script/mean_reference.py` | the deterministic reference |
| `log/166-acceptance.json` | the current `stage2-tenv-acceptance-v2` record |
| `log/165-acceptance.json`, `log/resources/{prnet,mean}-fold_*-seed_*.time` | earlier acceptance record and raw `/usr/bin/time -v` measurements |
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
  loss). In all four runs the 5-epoch rule stopped the run.
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
  confidence interval** and not independent biological replicates; the figures
  show their min–max band, not an uncertainty estimate.
- The evaluated population is the declared matrix's 918 test groups (169 null),
  not the corpus-wide 1,377 groups (273 null). Fold 2 is not run.
- The label side is `padj < 0.05` with direction `sign(log2FC)`, so the median
  group has 11 measured effects and 19.8% of the corpus has none. PRnet's
  measured-effect recovery (`changed Reliability = 0` at `d*`) is correspondingly
  low; a peak at the strictest threshold is a finding about the entrant, not a
  failure of the scan.
- `d` is a predicted magnitude in ln-log1p CP10K units, not a fold change or a
  p-value. The predicted and measured criteria are asymmetric by construction:
  the truth side keeps a genuine FDR while the predicted side keeps a declared
  magnitude.
- Direct metrics are single-well statistics over the 18,164-gene axis and carry
  TempO-Seq well-level noise; they measure per-sample fidelity, not group-level
  accuracy.
- CPU-only acceptance: no GPU was available this boot; CPU wall was 68.07–105.38 s
  per PRnet run with ≤ 2000.8 MiB peak RSS (`Mean` 3.07–3.12 s, ≤ 1073.7 MiB).
  GPU timing is untested.
- The entrant is PRnet alone; no second entrant is promised.

## Retirement and historical boundary

The previous k-parameterized indirect family, its magnitude gate, its plot
values and its figure panels are retired by ADR 0064; they are not in the
working tree and Git history is their retention mechanism. The earlier sealed
records (`log/162`, `log/164`, `log/165`) keep their bytes and the superseded
[`log/165-acceptance.json`](log/165-acceptance.json) carries the forward
pointer to this record. The superseded comparison implementation and its
figures remain untouched and replayable from Git history
(`git show 1543781^:stage2/03_modes_tenv/...`); this step adds no dependency on
them and re-implements none of their code.
