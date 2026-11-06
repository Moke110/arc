# Select a source before designing each Stage 3 phenotype corpus

**Status: accepted strategy (2026-11-06); source selection research in progress**

The user reset the corpus-building approach after the universal representation discussion became too complicated. For each desired `(condition, phenotype)`, select the best source, build a corpus around that source, and only then find similar sources that can be converted and added to that existing corpus. Do not require a universal multi-source schema before selecting the anchor source. This replaces the universal-corpus-first approach in ADRs 0039–0043; earlier sample representation proposals are not implementation requirements for the new approach.

The initial search targets cardiomyocytes under chemical treatment plus baseline, hypoxia or hypoxia–reperfusion/reoxygenation, with `cell_viability` measured using CCK8 or MTT. Chemicals must be resolvable to canonical SMILES. Do not substitute ATP/CellTiter-Glo, beat count or imaging viability for these requested assays, or chemical-toxicity injury for an oxygen-condition match.

Rank eligible sources by chemical breadth first, dose coverage second, 24 h timing over other timings third, and cell-model preference fourth: human primary CM or human iPSC-CM, then animal primary CM, then subcultured cells. Record the actual treatment/oxygen schedule so 24 h pretreatment is distinguishable from total assay duration. Screen-library size is not equivalent to publicly recoverable chemical/readout coverage. Assess those facts before recommending an anchor.

The previous requirement to derive significance tests and endpoint labels remains removed. Source interpretation and provenance remain necessary; further schema choices should follow the selected source rather than anticipate every possible experiment.

See [ADR 0043](0043-stage3-source-measurement-corpus.md), [source feasibility research](../research/stage3-viability-sources.md), and [GLOSSARY.md](../../GLOSSARY.md).
