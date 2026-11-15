# AF stays out of the gene axis: the clinical CM-specific audit

**Status:** accepted (2026-11-15). Amends ADR 0069's rationale (dataset-level, not only the
informativeness rule) and ADR 0068 (informativeness is necessary but not sufficient) without
rewriting either record. Reference: the re-admission request of 2026-11-14, ADR 0066 *Source
eligibility*, the user rule that only data able to resolve cardiomyocytes specifically may
qualify, and the audit note `docs/research/stage2-af-source-audit.md`.

## Why this record exists

ADR 0069 removed `AF` because its only dataset carried no ordering evidence. The user then asked
to re-admit AF with `GSE255612` confirmed as a clinical left-atrial AF snRNA source after the
"AF data search" session. Re-admitting on provenance alone would repeat ADR 0068's mistake in the
opposite direction: the dataset exists and is clinical, but a source must carry ordering evidence.

So the whole AF corpus found by that search was audited before deciding: every candidate was
cross-checked against the tracked repo, filtered by the CM-specificity rule, its processed payload
downloaded, its actual donor-level DEG count computed, and its DE ordering tested against a
**label-permutation null** (100 permutations of the sample group labels, same model refit, top-`N`
overlap with the released axes). A uniform-random gene baseline was not used: it flatters any
dataset whose measurable genes happen to be large-effect genes.

## What the audit found

- 117 verified GEO series; **only 4 had any tracked repo record** (`GSE255612` item, `GSE128188`
  item + AF-catalog rejection, `GSE234141` used for `genetic_cardiomyopathy`, `GSE79768` in a
  research note only). The rest were never evaluated; under the CM rule almost all are ineligible
  (whole-tissue bulk, non-CM cell types, blood, microarray, controlled access).
- **CM-specific clinical candidates, all negative:**
  - `GSE255612` — 18 AF / 16 control, author CM label, 32,438 CM nuclei. This repo's pipeline: 0
    DEGs (min padj 0.7359; all 18,092 axis-measured effect scores 0). The paper's own CM table: 2
    DEGs (LINC01257, ATRNL1, adj p 0.0174), authors state power was insufficient. A count-based
    re-derivation gives 0/2/6 hits (voom/edgeR-QL/DESeq2) but the union is not credible AF biology
    and the canonical CM markers stay null. Label-permutation null: within-axis p=0.634 (top-2000)
    and p=0.347 (top-5000); genome-wide p=0.515/0.277 → the apparent axis enrichment is an
    expression-magnitude artifact.
  - `GSE138252` — CM-nuclear bulk (FANS PCM1+), LA 5 AF+HF vs 5 non-failing. 2 (limma-voom) / 6
    (edgeR QL) DEGs, but the hits track the HF arm (log2FC r=0.577 with HF-only; sex-adjusted 0),
    and the permutation null gives p=0.574 (top-2000) / 0.901 (top-5000, observed below null). The
    AF+end-stage-HF confound is not de-confoundable (HF-only arm is n=3 and batch-separated).
  - `GSE238242` — 4 SR / 3 AF donors, author CM label, 3,049 CM nuclei: 0 DEGs (min padj 0.217).
  - `GSE271229` — 3/3/3 donors: the deposit ships no cell-type annotation at all, so a CM gate
    cannot be built; ineligible.
  - dbGaP phs003644 (left-atrial CM RNA-seq): controlled access.
- **Only tier-2 in-vitro models carry large DEG lists:** `GSE175944` (PITX2−/− atrial EHT, 586),
  `GSE255991` (ATRNL1 siRNA/OE, 1,963/7,783), `GSE307021` (TBX5-degron, 302). They are perturbation
  models of AF GWAS genes, not clinical AF.

## Decision

1. **AF stays out of `01_bases`.** ADR 0069 stands; the audit gives the dataset-level reason ADR
   0069 lacked. `AF` is a declared gap like `HHD`/`PPCM`, documented by
   `docs/research/stage2-af-source-audit.md`.
2. **The audit becomes visible in the rule document.** `gene_axis_fill_rule.docx` gains a **source
   dataset information table** at its start: one row per released dataset (category, accession,
   sample type, groups, sample sizes, sequencing technology, published CM DEG count), followed by
   a block listing the audited-and-excluded AF sources with their actual DEG counts. The DEG count
   is the published `catalog.json` `n_significant` (one row per axis gene, best contrast,
   padj<0.05); item-level row counts are a different quantity and are never mixed in.
3. **Accession and sequencing technology are verified at the primary source.** Each row's
   accession and technology are looked up against the primary record (NCBI GEO series and platform
   records, cellxgene collection) and kept in a tracked per-dataset metadata mapping read by the
   docx generator; the release catalogs are not re-published for it, so the release manifest and
   views are untouched by this table.
4. **ADR 0068 is amended.** Informativeness is necessary but not sufficient: a dataset or category
   newly admitted to the fusion must also beat a label-permutation null. The rule is not applied
   retroactively to the six sources already released; that limitation is recorded, not hidden.

## Rejected alternatives

- **Re-admit `GSE255612` as a fused source** (existing table, or a count-based re-derivation): no
  AF-specific ordering; the count-based table would pass ADR 0068 while carrying no AF signal,
  which is exactly the failure mode the rule exists to stop.
- **Re-admit it as an evidence-only source** (the ADR 0068 state): ADR 0069 explicitly rejected
  carrying a permanently empty source, and the audit strengthens that rejection.
- **Use `GSE138252` as the AF source:** AF+HF confound, sex-confounded RA arm, permutation null.
- **Use `GSE238242` / `GSE271229`:** null at n=4/3; no resolvable CM annotation.
- **Redefine AF as a hiPSC-CM model category** (`GSE175944`/`GSE255991`/`GSE307021`): credible
  perturbation models, but the category would describe AF genetics, not AF disease.
- **Wait for dbGaP phs003644:** controlled access; re-admission is a separate decision if access
  is ever granted.

## Consequences

- `01_bases` keeps five categories, six fusion sources and `pN ∈ 0..6`; the AF sources are
  documented, not published.
- A future AF re-admission needs a CM-specific AF dataset whose ordering beats the permutation
  null; the audit note records what that would take.
- The docx table adds a tracked metadata mapping and two catalog inputs to the docx replay test,
  but no release product changes.
- The "AF has no source" premise behind ADR 0069 is corrected to "AF has sources and none carries
  AF-specific CM ordering".

## See

- [0066](0066-stage2-corpus-view-and-gene-axis.md) — source eligibility and the category roster
- [0068](0068-fusion-amendment-degenerate-sources-and-rank-norm.md) — the non-informative rule this
  record amends
- [0069](0069-remove-af-category-from-01-bases.md) — the removal this record confirms
- `docs/research/stage2-af-source-audit.md` — the full audit (candidates, methods, numbers)
- `docs/research/stage2-gene-axis-cm-source-survey.md` — the earlier CM-source survey
- `stage2/01_bases/README.md` — the five-category contract and the docx's source table

## Erratum (2026-11-15)

The Consequences section's "two catalog inputs" was an early estimate. The implemented docx
replay reads the **five** `output/disease-deg/<category>/catalog.json` inputs plus the tracked
`stage2/01_bases/execution/gene_axis_source_metadata.py` module (issue #197), and the document
under test is `output/gene-axis/gene_axis_fill_rule.docx`. The decision above is unaffected; this
accepted record is not rewritten.
