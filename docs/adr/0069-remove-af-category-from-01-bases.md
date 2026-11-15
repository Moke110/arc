# Remove the AF category from Stage 2 `01_bases`

**Status:** accepted (2026-11-14). Amends ADR 0066's six-category roster and ADR 0068's
`AF` example without rewriting either record; accepted records are historical. Reference:
the user decision of 2026-11-14 ("given that AF has no source, remove it from 01"),
ADR 0066 *Disease categories*, ADR 0068 *Why this record exists*, and ADR 0030 for the
external item layout.

## Why this record exists

ADR 0066 declared six cardiomyocyte-derived disease categories, including `AF`
(`disease_data/GSE255612`, 18 AF / 16 control left-atrial CM donor pseudobulk). ADR 0068
then established that its only dataset, `gse255612__af`, carries **no ordering evidence**:
all 18,092 measured genes have `padj ≡ 1`, so every effect score `|log2FC|·−log10(padj)`
is exactly 0 and the within-dataset rank degenerates to the `gene_index` tie-break. ADR 0068
excluded the category from the fusion but kept its DEG/rank table published as evidence.

The result was a declared source that contributes nothing to any axis but still occupies
release surface in 01: an always-empty `af_rank` candidate column, an empty `%af_gN`
coverage cell, zero composition rows, `fusion.excluded_categories: ["af"]`, and a published
3.1 MB disease-DEG product that no axis, view or model reads. The user judged that an empty
source should not be carried inside 01 at all.

## Decision

Remove `AF` from `01_bases` entirely:

1. **Disease-DEG release.** Delete `output/disease-deg/af/` (its `catalog.json`,
   `manifest.json`, `gse255612__af.deg.csv` and `gse255612__af.rank.csv`) and drop `AF` from
   the DEG source declaration (`execution/disease_deg_sources.py`), so
   `build_disease_deg.py publish|verify --all` covers the three remaining #187/#188
   categories (`nHF`, `ICM`, `genetic_cardiomyopathy`).
2. **Gene axis.** Drop `af` from the `bases.gene_axis.CATEGORIES` roster: the candidate
   table loses its contracted `af_rank` column (12 columns), `axis_composition.csv` loses
   its `af` rows (40 rows), `source_topn_coverage.csv` loses `%af_gN` (7 columns), and the
   release manifest has five categories, `fusion.n_sources: 6`,
   `fusion.excluded_categories: []` and `release.pN_bounds: [0, 6]`.
3. **Rule kept.** The non-informative-source rule of ADR 0068 remains implemented and
   tested (hermetic fixtures still exercise a non-informative dataset and a fully
   non-informative category); it simply has no live instance in the released release.
4. **External item untouched.** `disease_data/GSE255612/` is not modified or deleted
   (ADR 0030): 01 stops declaring and publishing its DEGs; the item, its `source/` tree and
   its tracked summary remain for provenance and re-admission.
5. **Downstream re-pin.** Both views pin the gene-axis release manifest as an input, so they
   are re-published and re-verified; their payload digests are **unchanged** because the
   released axes and rank sums are unchanged.

## Why the axes do not move

`AF` contributed a constant worst-or-tied rank to every candidate, so removing it changes
only the schema and the metadata, not the ordering: `gene_axis_2000.csv` and
`gene_axis_5000.csv` keep their digests (`cb07535f…` / `15e5c49e…`), the intersection stays
1,972, `rank_sum` stays 0.0237–5.8103 (median 3.0118) and every `pN ∈ 0..6` distribution is
unchanged. The candidate table, composition table, coverage table and manifest change.

## Rejected alternatives

- **Keep the AF DEG table as evidence** (the ADR 0068 state): the user rejected carrying an
  empty source in 01; the non-informative evidence is preserved in ADR 0068's text and in
  the item's own tracked summary.
- **Keep `af` declared but always excluded**: identical empty release surface (empty column,
  empty coverage cell, zero composition rows) with no scientific gain.
- **Delete `disease_data/GSE255612`**: out of scope — external items are ADR 0030 read-only
  provenance, not 01 products, and a future decision may re-admit AF with a different
  contrast or dataset.

## Consequences

- `01_bases` no longer publishes or consumes AF; `build_gene_axis.py` requires exactly the
  five declared categories, and `output/disease-deg/` has twelve datasets in five
  directories.
- The released gene-axis schema changes (candidate columns 13 → 12, composition rows
  42 → 40, coverage columns 8 → 7); consumers must read the columns from the manifest, not
  hard-code `af_rank`.
- ADR 0068's `AF` example becomes historical; the rule it states is unchanged.
- The 2026-11-14 gene-axis evidence also replaces the ambiguous manifest phrase
  `"ascending is better"` with `"smaller is better (a descending value is better)"`; that is
  a wording clarification of the same rule, not a second decision.

## See

- [0066](0066-stage2-corpus-view-and-gene-axis.md) — the amended six-category contract
  (text unchanged)
- [0068](0068-fusion-amendment-degenerate-sources-and-rank-norm.md) — the non-informative
  rule kept here, whose `AF` example this record retires
- [0030](0030-external-data-provenance-and-step-layout.md) — why the item tree is untouched
- `stage2/01_bases/README.md` — the current five-category contract, products and digests
- `stage2/01_bases/log/acceptance.md` — the replay/digest evidence of this removal
- `stage2/01_bases/output/gene-axis/manifest.json` — categories, `fusion`, coverage columns
