# Publish only 5k expression views with axis-local normalization

**Status:** accepted design (2026-11-18); data implementation and acceptance pending.

The user needs only tenv-5k and tahoe_bulk-5k expression inputs, normalized within the genes each model receives. Retain the exact ADR 0073 nested 2k/3k/4k/5k axes, including original frozen `gene_index`, `symbol` and candidate/composition/coverage evidence, but retire both sources' 2k/3k/4k expression views. This amends ADR 0066's full-axis normalization choice and ADR 0073's four-view obligation; their original measured releases and sealed evidence remain history. Until coherent successor acceptance, the existing release remains authoritative.

For each treated or vehicle corpus row, first select the unchanged G_5000 integer counts, then calculate `log1p(10000 * counts_5k / sum(counts_5k))` using the natural logarithm. Normalize each matched vehicle separately on its own 5k total, average its matched group to obtain `control`, and calculate `delta = X - control`. Keep sample identities, roles, ordering, matched-control relationships and float32 expression publication. A non-positive selected total fails with source/profile identity; never drop the row or introduce a pseudocount. Build from corpus integers, not old log expression.

The trade-off is deliberate: model inputs use a selected-axis denominator, at the cost of losing the previous cross-axis normalization invariance and comparability to old expression/RMSE scales. Publish a new schema/recipe generation with explicit denominator, original axis identity and per-profile totals; reject mixed generations. Use existing bounded readers and deterministic candidate/replay publication. Validate both candidate sources before replacing the live 5k release and removing small-view products. Active builders, inventories and human-document population lookups must use retained products, and must not recreate retired expression views. Preserve upstream scientific bytes and all four named axes; no new axis selection is authorized.

This change is confined to Stage 2 01. Existing Stage 2 02/03 frozen inputs, weights, results and acceptance evidence retain their original full-axis normalization and remain historical evidence for that generation. New input preparation must identify the new generation explicitly; retraining or rescoring requires separate work. Live metadata and documentation must distinguish those generations without rewriting sealed records.

Acceptance uses analytical small-corpus checks through the public reader/build boundary, non-selected-count invariance, separately normalized matched vehicles, non-positive-total failures, schema/digest rejection, exact axis/population preservation and independent replay of both real 5k sources. Run repository health checks under the established resource limits and retain a new coherent receipt. This planning publication performs no data removal, payload rebuild or training.

## See

- [ADR 0066](0066-stage2-corpus-view-and-gene-axis.md)
- [ADR 0073](0073-full-source-tier-first-deg1000-and-four-axes.md)
- [Glossary](../../GLOSSARY.md)
- [Bases](../../stage2/01_bases/README.md)
- [Execution constraints](../agents/execution-environment.md)

- [Specification #243](https://github.com/Moke110/arc/issues/243)
- [Candidate implementation #244](https://github.com/Moke110/arc/issues/244)
- [Publication and retirement #245](https://github.com/Moke110/arc/issues/245), blocked by #244
