# Select a source before designing each Stage 3 phenotype corpus

**Status: accepted strategy (2026-11-06); amended below to include CTG**

The user reset the corpus-building approach after the universal representation discussion became too complicated. For each desired `(condition, phenotype)`, select the best source, build a corpus around that source, and only then find similar sources that can be converted and added to that existing corpus. Do not require a universal multi-source schema before selecting the anchor source. This replaces the universal-corpus-first approach in ADRs 0039–0043; earlier sample representation proposals are not implementation requirements for the new approach.

The initial search targets cardiomyocytes under chemical treatment plus baseline, hypoxia or hypoxia–reperfusion/reoxygenation, with `cell_viability` measured using CCK8 or MTT. Chemicals must be resolvable to canonical SMILES. Do not substitute ATP/CellTiter-Glo, beat count or imaging viability for these requested assays, or chemical-toxicity injury for an oxygen-condition match.

Rank eligible sources by chemical breadth first, dose coverage second, 24 h timing over other timings third, and cell-model preference fourth: human primary CM or human iPSC-CM, then animal primary CM, then subcultured cells. Record the actual treatment/oxygen schedule so 24 h pretreatment is distinguishable from total assay duration. Screen-library size is not equivalent to publicly recoverable chemical/readout coverage. Assess those facts before recommending an anchor.

The previous requirement to derive significance tests and endpoint labels remains removed. Source interpretation and provenance remain necessary; further schema choices should follow the selected source rather than anticipate every possible experiment.

See [ADR 0043](0043-stage3-source-measurement-corpus.md), [source feasibility research](../research/stage3-viability-sources.md), and [GLOSSARY.md](../../GLOSSARY.md).

## Assay scope amendment (2026-11-06)

The user subsequently included CellTiter-Glo (CTG) as another assay for `cell_viability`, alongside CCK8 and MTT. This replaces the CTG exclusion above. Preserve assay method, source scale and normalization context; sharing a phenotype does not imply direct numerical interchangeability of assay values.

Under the existing chemical-first ranking, Heliste 2021 becomes the leading currently assessed anchor candidate for both baseline and hypoxia–reoxygenation: 689 source chemical labels at four doses (0.04, 0.2, 1 and 5 µM) with CTG readouts and an archived numerical supplementary workbook. Chemical labels are not yet a count of unique SMILES-resolved structures. Its schedule is 24 h pretreatment, 48 h hypoxia and 24 h reoxygenation, with a duration-matched normoxia arm; it must not be described as a 24 h endpoint. The corresponding provenance assessment places total elapsed time at about 96 h. This update does not establish a pure-hypoxia anchor, complete chemical resolution or approve cross-source conversion.

Source: [Heliste 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8371130/), [archived source description](../../pheno_data/README.md).

## Imaging viability scope amendment (2026-11-06)

The user also included imaging-based viability in `cell_viability` and requested acquisition/feasibility assessment of Huang 2022 as a baseline anchor. This replaces the earlier imaging exclusion. Whole-well Hoechst/propidium-iodide live/dead quantification is therefore eligible; retain its measurement definition and do not substitute functional imaging endpoints such as calcium or beat rate. Huang's reported 2,374/2,375-compound, 13-donor, 24 h screen is a leading design-fit candidate, but its full numerical data must be obtained before corpus construction can be claimed. See the [access audit](../research/stage3-huang-2022-access-audit.md).

## Heliste anchor selection (2026-11-06)

The user selected Heliste 2021 to build separate baseline and hypoxia–reoxygenation cardiomyocyte chemical-perturbation `cell_viability` corpora now. The requested Stage 3 step layout is `stage3/01_corpus/output/baseline/` and `stage3/01_corpus/output/H-R/`. This replaces the earlier `01_data_corpus/output/corpus/<release-id>/` layout for the new build; migration of the historical step is still to be designed. Broader sources remain potential later additions rather than prerequisites for this anchor.

The user accepted replacement of the existing `01_data_corpus` step by `01_corpus`, retaining historical work through Git. Each condition directory contains `corpus.parquet`, a readable `corpus.csv`, `source_annotation.json` and a chemical-resolution report. One record represents one supplied chemical–dose–condition value, with `rep: null` and `readout: {"source_reported": value}` because the workbook does not establish physical replicate identity or a group-mean interpretation. Preserve the author-normalized scale and unresolved normalization details in the annotation. Treatments require evidence-backed canonical SMILES; unresolved identities remain in the report pending resolution, without invented structures or control rows. Existing endpoint-gated code and the historical Heliste `empty-after-QC` assessment do not establish admission rules for this source-first build.

## Desalted parent amendment (2026-11-06)

The user requested desalted forms for all Heliste chemicals, including the six
previously pending labels. This replaces the first implementation's preservation
of reference salts in corpus canonical SMILES. Apply a deterministic, versioned
parent transformation to every resolved reference: remove counterions and
solvates, normalize removable protonation, retain stereochemistry and permanent
charges. Active metal coordination centers and their ligands remain part of the
parent even when a reference SMILES encodes them as disconnected components;
a largest-fragment rule must not turn a coordination drug into an isolated ligand.
Separate genuine metal counterions require supported parent evidence.
Keep original labels, reference structures and resolution evidence so
the administered form remains auditable. Supplied molar doses remain unchanged;
parent mapping does not authorize an inferred dose or potency conversion.

Reference candidates may be admitted when evidence supports the supplied
chemical identity and all relevant candidates converge to one desalted parent.
Desalting alone does not resolve stereoisomer mixtures, contradictory identities
or a polymer without a defensible finite molecular representation. Keep these
explicitly pending rather than select an arbitrary structure. Preserve separate
source rows even when parents coincide, and report source-label and unique-parent
counts separately. Freeze new identity evidence as a new declared local addition;
retain the previous snapshot and immutable source bytes. Rebuild both condition
exports through the existing offline build/verify seam and retain complete
5,512-cell accounting.
