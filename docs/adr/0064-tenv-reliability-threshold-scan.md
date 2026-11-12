# Score the tenv comparison with a predicted-effect threshold scan

**Status:** accepted (2026-11-12). Supersedes the indirect-metric and figure
portions of ADR 0063. Keeps ADR 0063's direct metrics, truth gate source, fold
and seed index tables, stage palette, `Mean` reference and PRnet entrant.

The user rejected the ADR 0063 indirect family as too complex and not aimed at
the actual goal. That goal is: emit the gene table a model would hand to the
downstream phenotype task, and score how closely that table matches the table a
real cardiomyocyte Perturb-seq experiment would produce. ADR 0063 instead
reported four benchmark statistics whose defined-group sets did not coincide.

## Why the @k family is retired

The accepted ADR 0063 records measured it. Every fold has 459 test groups, but
91 (fold 0) and 84 (fold 1) have an empty gate. `max_k@P80S40` and
`Recall@max_k` were defined for 0–4 groups per run, and the PRnet and `Mean`
defined-group sets were disjoint in every fold/seed, so those two macro numbers
compared different handfuls of groups. A metric whose population is chosen by the
model being scored cannot rank models.

Two further defects were found while redesigning: `Precision@k` was NaN on an
empty gate (ADR 0063's prose said zero), and `reasons['Spearman@k']` recorded
`k_below_3` even when the `k >= 3` values were defined.

## Goal restated

A Perturb-seq contrast yields, for each measured gene, a direction
(`up` / `not significant` / `down`), and its downstream consumers take a ranked
up/down table. The model must therefore reproduce **both** the set of changed
genes and their order. The prediction side has no valid p-values, so the two
sides can never use the same criterion; ADR 0064 accepts that asymmetry and
scores the predicted side on its own effect magnitude.

## Boundary

`01_bases` owns the metric computation, the drawing API and the panel set.
`02_models_tenv` owns the runs, their scoring, the plot-value tables and the
figures. The model, its training and its per-well inference are unchanged: the
redesign rescores already-published prediction artifacts and re-renders figures.
No differential-expression tool runs in the prediction path.

## Predicted side

A group's predicted side is its **predicted absolute rank table**: all 18,164
genes ordered by `|Δ̂_mean|` descending, with direction `sign(Δ̂_mean)`. This is
the existing product of ADR 0063 and is reused unchanged.

## Label side

A **measured effect** is a gene whose accepted `tenv-truth` DESeq2 contrast has
`padj < 0.05`, with direction `sign(log2FC)`. The `|log2FC| >= 2` gate of ADR 0063
is dropped: it existed only to define a "final table" that this record no longer
needs, and it removed 5.4 genes and 10 non-empty groups from an already small
label set (median 11 measured effects per contrast, q25 = 1, max = 7,457).

The evaluated population is **every one of the 1,377 parent-dose groups**, not
only the 1,094 with a non-empty label set. The 273 groups with no measured effect
are real null results — 19.8% of the corpus — and a model that hallucinates a
signature on them must be penalised. ADR 0063 excluded them, which is what made
the weak-condition population shrink.

## Threshold scan

A group's **threshold window** at threshold `d` is `W(d) = {i : |Δ̂_i| >= d}`.
The window is the model's own changed-gene table, so its size follows the
group's predicted effect strength instead of a declared table size. A fixed table
size `k` was rejected by the user because it forces a weak condition to emit as
many genes as a strong one, and because it caps a weak condition's attainable
`F1` at `2/(k+1)` regardless of model quality.

The scan is 100 points: `d = k · 0.01 · M` for `k = 1..100`, where `M` is the
run's largest `|Δ̂_mean|` over its test groups and the whole gene axis, recorded
per run. The plotted axis is the scale-free `d' = d/M = k/100`, so the four runs
of an entry share one grid despite their different magnitudes. Normalising each
run by its own `M` is deliberate: it removes absolute magnitude calibration from
the comparison, so a model is judged on the best table it can produce on its own
scale rather than on whether its Δ happens to be in the right units.

The comparison is `|Δ̂| >= d`, not `> d`, so the strictest grid point can still
select each group's strongest gene instead of emptying every window.

## Totality

The window and the label set can each be empty. Rather than the NaN-plus-reason
machinery of ADR 0063, the four grid cells are given values, so every group
contributes to every macro average and no metric has an undefined-group set:

| `W(d)` | measured effects | `Precision@d` | `Recall@d` | `F1@d` |
|---|---|---|---|---|
| empty | none | 1 | 1 | 1 |
| empty | present | 0 | 0 | 0 |
| non-empty | none | 0 | 0 | 0 |
| non-empty | present | `|H|/|W|` | `|H|/|M|` | harmonic mean |

The empty/none cell is what lets a null group reward a model that correctly emits
no signature. The non-empty/none cell penalises hallucination on the same group.

## Metrics

A **direction-correct hit** `H(d)` is a gene of `W(d)` that is a measured effect
with matching direction.

- `Precision@d = |H|/|W|` — of the genes emitted, how many are real and correctly
  directed. A direction-flipped gene is both a false positive and a missed
  effect; direction is therefore carried by `Precision@d` rather than by a
  separate metric.
- `Recall@d = |H|/|measured effects|` — of the measured effects, how many were
  emitted. The denominator is direction-agnostic so that a correct table can
  reach 1.
- `F1@d` — the harmonic mean of the two. Harmonic averaging is load-bearing: an
  arithmetic mean of `Precision@d` and `Recall@d` peaks at the smallest window
  (`Precision = 1`, `Recall = 1/|effects|`), which would make the selection below
  choose a one-gene table.
- `Spearman@d` — inside `H(d)` only, the Spearman correlation between predicted
  magnitude rank (`|Δ̂_mean|`) and measured magnitude rank (`|log2FC|`). It asks
  the one question `F1@d` cannot: given that these genes are right and correctly
  directed, did the model order them by effect size?
- `Reliability@d = F1@d · mult`, with `mult = (1 + Spearman@d)/2` when
  `|H| >= 3` and `mult = 1` otherwise. Range `[0, 1]`.

`Spearman@d` is conditioned on the model's own hit set, so it is easier for a
conservative model and is not a standalone measure of ranking ability. The
neutral `mult = 1` below three hits is required for correctness, not politeness:
scoring it `0.5` would cap every weak condition at half its `F1` and reintroduce
the weak-condition penalty this record exists to remove.

## Selection and the delivered table

A run's macro `Reliability@d` is the equal-weight mean over its fold's 459 test
groups. The user selects models by the **peak** `Reliability@d`, and each run
carries its own `d*`, so a consumer of a given run's checkpoint uses that run's
`d*` and its delivered table `W(d*)`. No signature table is published as a
stable product in this decision; `d*` and `Reliability@d*` are recorded values.

## Figures

`figure/metrics.png` retains the four direct-metric panels, one row by four
columns. A new `figure/metrics_indirect.png` is two rows by four columns:

1. `Precision@d`, `Recall@d`, `Spearman@d` and `Reliability@d` over `d'`, four-run
   mean line and four-run min–max band, with `F1@d` drawn as a thin reference
   line on the Reliability panel and `d*` marked;
2. `Reliability@d*` and `d*` bars of the four-run mean with one dot per run, a
   macro three-class confusion heat map (`up` / `not significant` / `down`) at
   `d*`, and the Reliability split between the null and the non-null groups.

`AUPRC` reporting is withdrawn. It cannot be defined on the 273 null groups, so
keeping it would reintroduce a second, model-independent group set that this
record exists to remove; the threshold-free comparison it offered is available
from the `Reliability@d` curve itself.

## Consequences and limits

- Peak `Reliability@d` is an in-sample maximum: the same group data selects `d*`
  and reports the value, so the headline number is optimistic. The whole curve is
  reported beside it so the optimism is visible, and the four per-run `d*` values
  are recorded as stability evidence.
- Different runs deliver different table sizes, so the downstream phenotype
  comparison must treat table size as part of the intervention, not as noise.
- `d` is a predicted magnitude in ln-log1p CP10K units. It carries no statistical
  test, so the predicted and measured criteria are asymmetric by construction;
  the truth side keeps a genuine FDR while the predicted side keeps a declared
  magnitude.
- The label side inherits `n = 3` treated wells, so it inherits the corpus's low
  power: 19.8% of contrasts have no measured effect and the median has 11.
- The three-class confusion panel is a diagnostic and is not scored.

## Retired vocabulary

`Precision@k`, `Spearman@k`, `max_k@P80S40`, `Recall@max_k`, `gated truth gene`,
the `|log2FC| >= 2` gate, the NaN-plus-reason undefined machinery, the batched
`Cpred`/`Ceval` design and the `b`-prefixed metric names are all retired by this
record. Git history is their retention mechanism.

## See

- [Direct and indirect metrics for the tenv comparison](0063-tenv-comparison-direct-and-indirect-metrics.md)
- [Shared Stage 2 bases and indexed corpora](0061-stage2-shared-bases-and-indexed-corpora.md)
- [Stage 2 accepted products and reproducible figures](0062-stage2-accepted-products-and-reproducible-figures.md)
- `stage2/AGENTS.md`, `stage2/01_bases/AGENTS.md`, `stage2/02_models_tenv/AGENTS.md`
- `GLOSSARY.md` — the @d reliability family entries
- [Execution environment](../agents/execution-environment.md)
