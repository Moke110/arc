# Gene-axis v2: tiered source collection, DEG-only rank fusion and the f2000 fill rule

**Status:** proposed (2026-11-15). Supersedes ADR 0066's *Source eligibility*,
*Priority*, *Candidate table and the two-level rank-sum fusion*; keeps ADR 0068's
non-informative-source rule; re-admits `AF`, superseding ADR 0069 and ADR 0070.
Reference: the user decision of 2026-11-15 ("DEG 用 padj<0.05 & FC>=4；按 tier 顺序搜集，
tier1/2 累计 ≥2,000 DEG 就停；只让 DEG 且 disease 内前 2000 进下游；轴用 f2000 而不是 pN").

## Why this record exists

The released axis (ADR 0066, 0068) builds each source's rank over **all measured genes**:
every measured gene gets a position from the effect score `|log2FC|·−log10(padj)`, so a
dataset with almost no significant signal still supplies a complete `1..n` ordering. The
released `GSE270788` arms illustrate the failure: their within-dataset ranks are a near-noise
ordering (7 / 2 / 1 genes at BH<0.05; the PICM arm's top-5,000 contains 3 significant genes),
yet the source still carries 1/6 of `rank_sum` and the source-coverage figure credits the
`ICM` source with 2,350 of the 5,000 axis genes. Measured on the released products, 568 of
those 2,350 genes (11.4% of the axis) are there only because the axis was selected using the
same rank; the figure's bars are therefore not a signal measure, and the docx's "CM DEG count"
column (1) and the figure's count (2,350) describe different objects.

The user's decision is to make the disease side **DEG-only**: a source contributes only genes
that pass a declared differential-expression threshold, and only its top-2,000 DEGs enter
downstream fusion.

## Decision

1. **DEG definition.** A gene is a DEG in a dataset when `padj < 0.05` **and**
   `|log2FC| ≥ 2` (FC ≥ 4) in at least one of the dataset's contrasts; the multi-contrast
   collapse keeps the maximum-effect contrast (ADR 0066's collapse rule, unchanged). `NA`
   `padj` is non-significant. Direction is carried, not filtered on.
2. **Tiered source collection.** Sources are collected per disease in tier order:
   - **tier 1** — human atria/ventricle plus cardiomyocyte specificity (single-cell/single-nucleus
     CM pseudobulk, or CM-nuclear bulk such as FACS-sorted nuclei);
   - **tier 2** — hiPSC/hESC-CM cells, organoids or engineered cardiac tissue (bulk). For the
     clinical categories (`nHF`, `ICM`, `AF`) tier 2 admits only models of the clinical,
     non-genetic form: mutation-modelled iPSC-CM belongs to `genetic_cardiomyopathy` under
     ADR 0066's isogenic exception, not to a clinical category;
   - **tier 3** — human atria/ventricle bulk tissue, not cardiomyocyte-specific.

   After finishing a tier, the disease's DEG union (resolved to the frozen axis) is measured;
   when it reaches **2,000** the collection stops and lower tiers are not entered. The tier,
   the per-tier union size and the stop/descend decision are recorded per disease.
3. **Per-source DEG rank.** Inside a dataset, DEGs rank `1..n_DEG` by descending effect
   `|log2FC|·−log10(padj)` with `gene_index` as the deterministic tie-break; a gene that is
   not a DEG in that dataset gets `n_DEG + 1`; `rank_norm = rank / n_DEG`. The denominator is
   the DEG count, never the measured-gene count.
4. **Per-disease fusion (level 1).** `disease_rank_sum = Σ_source rank_norm`; the disease rank
   is the ascending position with `gene_index` tie-break. `DEG2000_d` = the genes that are a
   DEG in at least one of the disease's sources **and** whose disease rank is `≤ 2,000` (all of
   them when the disease has fewer than 2,000 DEGs).
5. **Cross-disease fusion (level 2).** `rank_sum = Σ_disease rank_d / n_d` where `n_d` is the
   disease's DEG count and a gene outside `DEG2000_d` takes `(n_d + 1) / n_d`. `f2000` = the
   number of diseases whose `DEG2000_d` contains the gene. `f2000` replaces `pN`.
6. **Axis fill.** For each `N ∈ {2000, 5000}`: the first `0.5·N` genes are the tenv-HVG top
   `0.5·N` (unchanged); the second `0.5·N` are the remaining candidates ordered by **`f2000`
   descending, `rank_sum` ascending, `gene_index` ascending**. When fewer than `0.5·N`
   candidates have `f2000 ≥ 1`, the remaining slots are filled with tenv-HVG genes not yet
   selected, in tenv-HVG rank order — never with genes that are not a DEG in any disease. The
   number of such fill genes is recorded in the release manifest.
7. **AF returns.** `AF` is declared again with its CM-specific sources (`GSE255612`,
   `GSE138252`, `GSE238242`) as tier 1 — their DEG union is 2 genes — and therefore descends to
   tier 3 bulk atrial tissue (`GSE269882` as a new item; `GSE128188` already owned). The
   category is published with a tier-3 evidence type on those sources.
8. **Products.** Each published `<dataset>.deg.csv` carries a `deg` flag under the new
   definition and the dataset's DEG rank; the candidate table carries `f2000` (the `p2000`/
   `p5000` columns are removed); `axis_composition.csv`, `source_topn_coverage.csv`, the axes
   and the rule document are rebuilt; the rule document opens with the new build rule and its
   second table lists, per disease, the number of source datasets, the total DEG count and the
   `DEG2000` count. The source-coverage figure is rebuilt on `DEG2000` (see 0066's figure
   contract for the two-panel layout).
9. **Removals.** The `pN`-based candidate columns, the old `source_topn_coverage.csv` schema
   and any measured-gene rank product are removed with their readers and tests (ADR 0062/0066
   deletion convention); git retains the code.

## Rejected alternatives

- **Raise the significance floor and keep the measured-gene rank.** A 3-DEG source would still
  supply a complete 18,164-gene ordering. The failure is the rank's universe, not the threshold.
- **Keep `pN` and add `f2000` beside it.** Two selection signals would have to be reconciled in
  every consumer; the user chose one (`f2000`).
- **Fill the second axis half with `f2000 = 0` genes ordered by the old `rank_sum`.** Rejected
  in favour of the tenv-HVG fallback (decision of 2026-11-15), so no gene enters the axis on
  disease evidence it does not have.
- **Keep `AF` out (ADR 0069/0070).** Superseded: the new rule tolerates small DEG sets, and its
  tier-3 descent gives AF a real bulk atrial contrast; the non-informative-source rule of
  ADR 0068 remains implemented for genuinely tied datasets.

## Consequences

- The disease half of both axes changes; `G_2000` / `G_5000` digests change, and both views are
  re-published on the new axes. Re-pointing `02_models_tenv` stays deferred (ADR 0066).
- The category rosters change: `AF` returns; `nHF` and `AF` gain tier-2/tier-3 sources; `ICM`
  gains the Kuppe MI `IZ vs CTRL` contrast (3,075 DEGs, 1,360 on the axis) and stops at tier 1.
- ADR 0066's sentence "whole-heart or whole-left-ventricle tissue bulk is **not** eligible" is
  superseded: such sources are tier 3 and are entered only when tiers 1–2 fall short of 2,000
  DEGs, with their evidence type recorded.
- The rule document's "CM DEG count" footnote is rewritten: the count is the new DEG definition,
  and the per-disease DEG2000 column is new.

## See

- [0066](0066-stage2-corpus-view-and-gene-axis.md) — the amended contract (text unchanged)
- [0068](0068-fusion-amendment-degenerate-sources-and-rank-norm.md) — the non-informative rule kept here
- [0069](0069-remove-af-category-from-01-bases.md), [0070](0070-af-stays-out-of-the-gene-axis.md) — the AF removals superseded here
- `stage2/01_bases/README.md` — products and contracts
- `docs/research/stage2-af-source-audit.md` — the AF source audit reused here
