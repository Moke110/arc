# Fusion amendment: non-informative sources, unmeasured `rank_norm`, and two removals

**Status:** accepted (2026-11-13). Amends ADR 0066's *Candidate table and the two-level
rank-sum fusion*, *Rank normalization* and *Consequences* without rewriting that record.
ADR 0066 remains the contract and is read together with this amendment; its text is not
edited, because accepted records are historical. Reference: spec #185 User Story 22
("a dataset that found almost nothing is not inflated to the weight of one that moved half
its transcriptome"), ADR 0066 *Candidate table*, and the post-code-review finding #194.

## Why this record exists

The two-axis review of the Stage 2 rebuild (#194) found that one declared fusion source
carried no ordering evidence. `output/disease-deg/af/gse255612__af` measures 18,092 genes
but has `padj ≡ 1` for every one of them (minimum raw p 1.4e-4; BH = 18,092·p ≈ 2.58), so
every effect score `|log2FC|·−log10(padj)` is exactly 0. Its within-dataset rank therefore
degenerated to the `gene_index` tie-break (Spearman 1.0 against gene index), yet it counted
as one of the seven equal-weight sources. That is the exact failure User Story 22 forbids.

The same review found a second, smaller deviation: `deg_rank.normalized_rank` clamped an
unmeasured gene's value to exactly `1.0`, while ADR 0066 *Rank normalization* defines
`rank_norm = rank / n_measured` with unmeasured `rank = n_measured + 1`, i.e.
`(n_measured + 1) / n_measured`.

This record amends the fusion for both, records two consequential removals, and resolves an
ADR 0062-vs-0066 removal-order conflict.

## Decision

### (a) A degenerate dataset is non-informative and is not a fusion source

A dataset is **informative** when its measured genes carry more than one distinct continuous
effect score. A dataset that measures genes but whose effect scores are **all tied** has no
informative ordering: its published rank is only the declared `gene_index` tie-break. Such a
dataset is:

- published **unchanged** — its DEG and rank tables stay in `output/disease-deg/`, because
  they are a faithful record of what the source measured;
- marked `informative: false` in its category `catalog.json` and release `manifest.json`,
  with an `informativeness` block recording the rule, the number of distinct effects and a
  human-readable reason;
- excluded from level-1 category fusion;
- and excluded from level-2 `rank_sum`, `pN` and `n_sources_measured`.

A **category whose datasets are all non-informative** is excluded from level-2 fusion and
from `pN` entirely. The released gene-axis manifest records `fusion.n_sources` (the count
actually fused), `fusion.sources`, `fusion.excluded_categories`, the per-category
`informative` flag and the exclusion reason. The candidate table keeps its contracted
`<category>_rank` column for a dropped category but leaves it empty (`NaN`), because the
category has no level-1 rank; the column is not silently repurposed.

**`pN` bounds follow the fused count.** With the current roster one category (`af`) is
dropped, so five categories plus tenv-HVG remain: `pN ∈ 0..6`, `rank_sum ∈ (0, 6]`, and the
rank-sum histogram is bounded by 6. The declared maximum remains seven sources; a future
category that becomes informative restores its weight and the bounds move with it.

**Rule, not dataset list.** The exclusion is data-driven and re-evaluated on every build from
the published effect tables. No dataset id is hard-coded as degenerate. The rule is narrow on
purpose: two distinct effect scores already count as informative; only a fully tied table is
dropped.

**Measured impact of this amendment (accepted rebuild, #194).**

| quantity | pre-amendment (ADR 0066) | post-amendment |
|---|---|---|
| fused sources | 7 | 6 (`af` dropped; five categories + tenv-HVG) |
| `pN` bounds | 0..7 | 0..6 |
| `G_2000 ∩ G_5000` | 1,919 (Jaccard 0.3777) | 1,972 (Jaccard 0.3922) |
| `G_2000` membership | — | 186 genes left, 186 entered, 1,814 shared |
| `G_5000` membership | — | 495 genes left, 495 entered, 4,505 shared |
| candidate `rank_sum` | 0.3845–6.5649 (median 3.5145) | 0.0237–5.8103 (median 3.0118) |
| candidate `p2000`/`p5000` max | 6 / 7 | 6 / 6 |
| `af` axis composition | 1,998 / 4,994 | 0 / 0 (excluded) |

The dropped source was not merely noise: its top set supplied 21.5% of `G_2000` and 39.0% of
`G_5000` (chance 11.0%/27.5%), so the amendment changes both axes materially. Both axes were
rebuilt from the amended candidate table, both views were rebuilt from the new axes, and all
replays are byte-identical.

### (b) Unmeasured `rank_norm` follows the ADR formula

`rank_norm = rank / n_measured` for **every** gene, measured or not. An unmeasured gene's raw
rank is `n_measured + 1` (unchanged, ADR 0066), so its normalized value is
`(n_measured + 1) / n_measured` — strictly worse than every measured gene and never clamped
to exactly `1.0`. The clamp let an unmeasured gene tie the worst measured gene and contradicted
the declared formula. All disease-DEG rank tables, the candidate table, both axes and both
views were rebuilt; the change touches only unmeasured rows.

### (c) Consequential removal: `bases/hvg_qa.py`

ADR 0066 *Removals* deleted `bases/metrics.py` as part of the evaluation chain. The retained
`bases/hvg_qa.py` imported `score_groups`/`summarize_runs` from that module to implement the
already-retired `verify_hvg_metrics` report, so it could not import after the deletion. It was
removed in commit `6078e9f` and its population-isolation assertion was dropped from
`test_hvg.py` (the retained `bases.hvg` release test keeps its two tests). This is a
**consequential** removal, not a new science decision: the module had no surviving purpose
once the metrics it reported were retired. Git retains the code.

### (d) ADR 0062-vs-0066 removal-order conflict for `02_models_tenv`

ADR 0062's acceptance rule is "candidate, then replace, then remove the superseded
products"; ADR 0066's *Removals* deletes the pre-0066 `training-corpora/` tree, which
`02_models_tenv` still read, while ADR 0066's *Deferred* says that step's re-point is out of
scope. Read literally, the two rules disagree about whether the removal may precede the
consumer's re-point.

**Resolution: ADR 0066's tree-specific removals govern.** The `01_bases` tree is removed on
its own acceptance, and `02_models_tenv`'s re-point to `view-2k/`/`view-5k/` remains deferred
and tracked as issue #193. ADR 0062's general ordering does not override a step-specific
removal decision; the consequence is a temporarily broken downstream reader, recorded rather
than hidden.

### (e) ADR 0067's released-axis view shares are superseded

ADR 0067 recorded released-axis CP10K shares of 12.9299987% / 37.2063189% for
`G_2000`/`G_5000`. Those were measured on the pre-amendment axes. After (a) the released
tenv shares are **15.0548953% / 37.4273748%** (treated-only 15.0892505% / 37.3846128%), as
recorded in both view manifests. The denominator argument of ADR 0066 is unchanged; only the
measured figures moved. ADR 0067's text is not edited.

## Consequences and limits

- The fusion is no longer "seven equal-weight sources". It is "the sources that carry
  ordering evidence", which is what User Story 22 asks for; the count and the bounds are
  published per release rather than assumed.
- A fully tied dataset can still be scientifically meaningful as a **negative** result. It
  stays published and marked; it just cannot weight the axis.
- The informativeness rule is deliberately coarse (two distinct values suffice). A dataset
  with near-degenerate ordering is not detected; that would need a declared effect-size
  criterion this record does not invent.
- ADR 0066's `pN ∈ 0..7` sentence and its "seven sources" consequence are amended for the
  current roster to `0..6`; the declared maximum of seven is retained for a future roster in
  which every category is informative.
- The AMENDMENT changes axis membership, so results computed on the pre-amendment `G_2000`/
  `G_5000` must be re-run against the new axes; the views carry the new gene sets and no old
  axis is retained.

## See

- [0066](0066-stage2-corpus-view-and-gene-axis.md) — the amended contract (text unchanged)
- [0067](0067-errata-to-0066-gse270788-arms-and-view-shares.md) — the errata whose view
  shares (e) supersedes
- [0062](0062-stage2-accepted-products-and-reproducible-figures.md) — candidate-then-replace
  acceptance, amended here for the removal-order conflict
- `stage2/01_bases/log/acceptance.md` — the measured replay/digest evidence
- `stage2/01_bases/output/gene-axis/manifest.json` — `fusion.n_sources`,
  `fusion.excluded_categories`, `release.pN_bounds`
- #194 (this fix), #185 (spec), #193 (deferred `02_models_tenv` re-point)
