# Rebuild the tenv model comparison with direct and indirect metrics

**Status:** accepted (2026-11-12). Supersedes the evaluation-metric and
figure portions of ADR 0042 and the reporting role of ADR 0059 for the rebuilt
comparison. Retains both historical implementations, their figures and their
records for replay. Introduces `stage2/02_models_tenv/` and the shared
evaluation, palette and figure products in `stage2/01_bases/`.

**Update (2026-11-12):** ADR 0064 replaces the indirect-metric and figure
portions of this record. The direct metrics, the fold/seed index tables, the
palette, the `Mean` reference and the PRnet entrant described below remain in
force; `Precision@k`, `Spearman@k`, `max_k@P80S40`, `Recall@max_k`, the
`|log2FC| >= 2` gate and the `metrics.png` indirect panels are retired.

The user removed every Stage 2 step except `01_bases` on 2026-11-11. That
foundation keeps the frozen gene axis, the training corpora, the measured tenv
truth and shared retrieval metrics, but no model comparison. This record
restores one comparison step on tenv and fixes what it reports.

## Boundary

`stage2/01_bases/` stays the foundation: it owns the shared metric computation,
the plotting APIs, the fold/seed index tables, the stage palette and the
evaluator truth. `stage2/02_models_tenv/` owns the trained entrants, their
predictions and the report figures, and instantiates nothing that duplicates a
`01_bases` implementation.

The first entrant is PRnet. The reference baseline is `Mean`. The former
`no_change` baseline is withdrawn. No other model is promised by this record.

## Evaluation target

Every metric is computed in **Δ space** (ADR 0061): the per-well target is the
well's observed `X − control`, where `control` is the mean `log1p CP10K`
expression of its matched same-plate DMSO wells. This is already the training
corpus `layers/delta`; it was verified exactly equal to `X − layers/control` on
all 4,206 wells × 18,164 genes with no missing or non-finite entry. No new
per-well Δ product is created.

A group (ADR 0061: one `(drug_key, dose)` pair) contributes its **Δ_mean**, the
mean of its wells' per-well Δ.

## Truth gate

A gene is a **gated truth gene** of a group when its DESeq2 contrast satisfies
`padj < 0.05` and `|log2FC| >= 2`. The direction of a gated truth gene is the
sign of its Δ_mean, which equals the sign of its measured log2 fold change. The
within-contrast BH convention, the matched same-plate control rule and the
fixed 18,164-gene family of ADR 0059 remain in force; only the fold-change
threshold is added.

The accepted `output/tenv-truth/` statistics are the sole truth source. The gate
is derived from them at evaluation time; it is not a new truth product. Truth is
evaluator-only and is restricted to the test groups of the fold being scored.

On tenv the gate leaves 65.5 gated genes per group on average, and 283 of the
1,377 groups have none. Because the median `|log2FC|` of significant genes is
about 6.2, the fold-change threshold removes 5.4 significant genes per group on
average and empties 10 further groups. It is a declared label rule, not a claim
about biology.

## Direct metrics

`PCC`, `Spearman`, `RMSE` and `R²` are **direct metrics**. Each is computed on
one test well: the predicted Δ̂ for that well against that well's Δ_truth, over
all 18,164 genes. No group aggregation is performed before reporting.

A run's direct-metric value is the equal-weight mean of its fold's test wells
(1,402 wells per fold). Per-well values are retained in the source tables.

A constant prediction or a constant truth makes PCC, Spearman and `R²`
undefined; the value is NaN with a reason, never zero.

## Indirect metrics

`Precision@k`, `Spearman@k`, `max_k@P80S40` and `Recall@max_k` are **indirect
metrics**, computed per group.

For a group, the predicted side first averages its several predicted wells into
a predicted `Δ̂_mean`. The group's **predicted absolute rank table** lists
`gene_index`, `rank`, `direction` and `|Δ̂_mean|`, ordered by `|Δ̂_mean|`
descending. `k` ranges over every integer from 1 to 100.

- `Precision@k` is the share of the top-`k` entries that are gated truth genes
  whose direction matches the prediction.
- `Spearman@k` re-orders the same top-`k` set by signed `Δ̂_mean` descending to
  give `rank_pred`, takes the true `Δ_mean` of those genes ordered descending to
  give `rank_gt`, and reports their Spearman correlation.
- `max_k@P80S40` is the largest `k` reaching both `Precision@k >= 0.8` and
  `Spearman@k >= 0.4`.
- `Recall@max_k` divides the direction-matched gated genes inside the top
  `max_k@P80S40` by all direction-matched gated genes of the group.

A group with no gated truth gene is undefined. `Spearman@k` is undefined for
`k < 3`. `max_k@P80S40` and its companion recall are undefined when no `k`
qualifies. Undefined values are NaN with a reason and an effective-group count;
none is filled with zero. Each indirect metric is the equal-weight macro mean
over the fold's defined test groups (459 test groups per fold).

## Run matrix and index tables

One run sample is one `(model, fold, seed)` result. The rebuilt comparison runs
**2 outer cold-drug folds × 2 training seeds**, giving four run samples per
entry. They are not independent biological replicates.

`01_bases/output/fold_n_seed/` holds the accepted assignment unchanged in
`cohorts.csv` and two index tables:

- `folds.csv` indexes the three frozen cold-drug folds in rows 0, 1, 2.
- `seeds.csv` indexes the training seeds `17`, `42`, `2026` in rows 0, 1, 2.

An experiment declares how many folds and seeds it uses and takes that many
rows from the top of each table. This comparison therefore uses folds 0 and 1
with seeds 17 and 42.

## Baselines

`no_change` is withdrawn from the rebuilt comparison. A constant-zero Δ leaves
the ranking family without a defined ordering — `Precision@k` is zero,
`Spearman@k` and `max_k@P80S40` are undefined — and leaves the correlation
family undefined. Reporting it as a bar would misrepresent an undefined entry as
a bad score.

`delta_mean` is renamed **Mean** with its definition unchanged: the mean Δ of
the same group in the training partition. It is a deterministic reference, not a
trained fit, and its fold copies are not independent training results.

## Stage palette and figure APIs

`01_bases/output/palette/` publishes the stage palette: six hue families (gray,
blue, green, orange, purple, magenta) crossed with nine lightness and saturation
shades, as a two-dimensional table of hex values together with its generation
recipe. Every Stage 2 figure selects colors from it, so a figure's color means
the same thing across steps.

Drawing APIs live in `01_bases/output/api/bases/plots.py`. Each step's
`figure/script/` reproduces its PNGs from `figure/data/*.csv` alone, with
`output/` unreadable, per ADR 0062.

## Figures

1. `metrics.png` — two rows by four columns, one subplot per metric. The two
   indirect panels plot `k` on the horizontal axis with a line through the
   four-run mean and a light band spanning the four-run minimum and maximum. The
   six direct and threshold panels are bar charts with one bar per entry, the
   bar height being the four-run mean and one dot per run on the bar. An entry's
   hue colors its line, bar edge and dots; the matching lighter shade fills its
   band and bar.
2. `learning_curve.png` — one subplot per learnable model that differs only by
   fold and seed (this study: PRnet alone). The horizontal axis is epoch and the
   vertical axis is the training and validation loss. The blue family is train
   and the green family is validation; one shade identifies one `(fold, seed)`
   run in both families.

## PRnet training

The entrant consumes upstream `models.PRnet.PGM` from `pert_model/PRnet/source`
without editing the source tree, through the same seam the retired adapter used
(a shared `Linear(2048, 1024)` Morgan-2048 adapter multiplied by `log10(dose+1)`
and a Gaussian NLL objective, with the adapter and predictor initialized and
trained jointly on training responses only). Baseline expression participates in
prediction and the actual dose is an explicit input.

Early stopping is epoch-based, not iteration-based:

- stop after **5 consecutive epochs** in which the validation loss does not
  strictly decrease;
- independently stop after **20 consecutive epochs** in which the validation
  loss is above the running minimum;
- never exceed **100 epochs**, and keep the checkpoint with the lowest
  validation loss.

The 5-epoch counter fires first in every trajectory that satisfies the 20-epoch
condition, so the second condition is declared and implemented but dominated.
It is kept because the user requires both conditions, and its dominance is
recorded rather than hidden.

The inherited ceilings remain: 14,400 s cumulative training, 1,800 s inference,
2 threads and the 10 GiB memory cap. The GPU is available again on this
workstation; CPU is the documented fallback rather than the plan.

## Consequences and limits

- Four run samples per entry replace the retired nine. A minimum–maximum band
  over four values is a spread indicator, not a confidence interval.
- A group with fewer than 13 gated genes cannot reach `Precision@10 = 0.8` even
  with every top gene correct, so `max_k@P80S40` is bounded by each group's own
  gated count. Macro means over groups with different gated counts are not
  directly comparable across models' worst groups.
- Direct metrics are single-well statistics and therefore carry TempO-Seq
  well-level noise. They measure per-sample fidelity, not group-level prediction
  accuracy.
- Withdrawing `no_change` amends the two-baseline requirement of ADR 0042. The
  rebuilt comparison keeps one deterministic reference; adding a further
  reference requires a new decision rather than a revived undefined baseline.
- ADR 0042's eight-metric table and ADR 0059's three-AUPRC reporting remain the
  contracts of their historical records. They are not the reporting contract of
  this step.

## See

- [Tenv model comparison redesign](0042-tenv-model-benchmark-redesign.md)
- [Three full-CM significance retrieval metrics](0059-stage2-significance-retrieval.md)
- [Shared Stage 2 bases and indexed corpora](0061-stage2-shared-bases-and-indexed-corpora.md)
- [Stage 2 accepted products and reproducible figures](0062-stage2-accepted-products-and-reproducible-figures.md)
- [External data provenance and step layout](0030-external-data-provenance-and-step-layout.md)
- `stage2/AGENTS.md`, `stage2/01_bases/AGENTS.md`
- [Execution environment](../agents/execution-environment.md)
