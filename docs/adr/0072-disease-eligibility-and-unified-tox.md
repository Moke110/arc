# Disease eligibility and a unified TOX category

**Status:** accepted and implemented (2026-11-16). The coherent five-category
release passed fresh TOX and downstream source replays and was promoted in #219.
The user confirmed the complete shared understanding and authorized specification
and ticket publication, with implementation/deployment deferred to a new session. Specification: [#213](https://github.com/Moke110/arc/issues/213).
Amends ADR 0071's disease roster,
small-signature admission and cardiotoxicity source unit, and ADR 0066's
mechanism-class cardiotoxicity admission. Other DEG, tiering and axis-fill rules
remain in force.

AF's released axis-resolved DEG union contains one gene (SLC9A7), while
cardiotox_TOP2 and cardiotox_other contain 4,514 and 297 respectively (ADR 0071).
We distinguish a declared disease's complete evidence trail from its eligibility
to contribute to cross-disease fusion. Drug cardiotoxicity becomes one disease
from source discovery onward, rather than two signatures merged after analysis.

## Decisions

1. **A uniform disease gate.** Complete source discovery, verification,
   selection, per-source DEG analysis, within-disease fusion and DEG2000 for each
   declared disease. Then count the disease's distinct DEG union resolved to the
   frozen 18,164-gene axis: fewer than 100 excludes the disease from cross-disease
   rank-sum fusion, f2000 and disease-driven axis selection; exactly 100 passes.
   This is not a per-source threshold. Retain the upstream products, counts,
   DEG2000 and explicit exclusion reason even when the disease fails. An empty
   union is excluded without attempting a zero-denominator rank calculation.
2. **AF remains a declared disease.** Preserve its complete collection,
   validation, screening, DEG, fusion and DEG2000 process. Its current one-gene
   union fails the new gate; its exclusion is a measured outcome, not removal of
   the category or deletion of its provenance. Apply the gate to future rebuilt
   results rather than hard-coding AF as permanently excluded.
3. **TOX is one disease from the beginning.** Replace cardiotox_TOP2 and
   cardiotox_other as disease definitions with TOX. Within each dataset, pool the
   admitted TOX drug exposures against their matched vehicle controls, retaining
   the existing source-specific tests, donor blocking and applicable cell-line,
   plate and time strata. Collapse contrasts with the existing DEG/significance
   and maximum-effect rules; each dataset supplies one source rank to TOX's
   level-1 fusion. Do not fuse the six old mechanism-class tables as a substitute
   for recomputing dataset-level TOX DEGs. Mechanism remains annotation only.
4. **Drug admission does not depend on mechanism classes.** Remove the
   per-mechanism minimum of three drugs and the requirement for a unique ATC/MOA
   class. Keep positive small-molecule identity, traceable cardiotoxicity evidence,
   experiment quality and matched-control requirements. Missing ATC/MOA alone is
   not an exclusion; antibodies and peptides remain outside this small-molecule
   definition. An unresolved chemical identity does not become eligible merely
   because the mechanism restriction was removed.
5. **Cardiotoxicity is an independently supported drug property.** Use the
   authors' explicit positive cardiotoxicity labels or traceable external evidence,
   recording source, date, rationale and uncertainty for each admitted drug.
   Panel membership, a chemical structure, ATC/MOA, or the new DEG result alone
   cannot establish positivity. Explicit negative-control drugs do not enter the
   TOX treated group; insufficient or conflicting evidence stays pending review.
   These drug controls do not replace matched vehicle controls.
6. **Admission is by drug identity, not a toxicity-positive exposure condition.**
   Include all doses and times of admitted drugs that satisfy the experiment and
   matched-control criteria. A demonstrated toxicity phenotype at every dose/time
   is not required. TOX denotes the transcriptional response to known cardiotoxic
   drugs, not proof of toxicity in every exposure. Retain the existing stratum
   collapse and weighting rules; no new drug-balanced weighting or directional
   agreement requirement is introduced by this decision.
7. **Reassess discovery and exclusions under TOX.** Review the existing candidate
   sources, searches and rejected drugs, supplement searches as necessary, and
   reuse valid downloaded data and verification evidence. Reassess GSE242692 and
   GSE313746, whose old single-drug exclusion no longer applies; neither is admitted
   without validation. Continue the tier ladder and cumulative axis-resolved
   2,000-DEG stopping rule. A second platform measuring the same experiment does
   not become an independent source merely because class restrictions disappear.

## Consequences and verification boundaries

With current non-TOX results and provided rebuilt TOX passes the gate, the eligible
roster is nHF, ICM, genetic_cardiomyopathy and TOX, so f2000 ranges from 0 to 4.
This is an expectation, not a measured new release. TOX's union, DEG2000 and final
axis composition must be recomputed; the old union sizes cannot predict them.

Publish explicit eligibility separately from the declared disease catalog so that
retaining AF does not silently reintroduce its rank term. Rebuild disease-DEG,
disease-fusion, both gene axes, both model-facing views, coverage tables, figures,
rule document and affected reader/documentation contracts. Keep the existing
DEG threshold (padj < 0.05 and |log2FC| >= 2), 18,164-gene identity, tier policy,
rank formulas for admitted terms, tenv-HVG allocation, nested axes and full-axis
view normalization. Corpus/count payloads are not changed by this design.

Acceptance must establish 99/100 gate behavior, retained excluded-disease
evidence, removal of excluded diseases from every cross-disease contribution,
dataset-level rather than class-level TOX sources, evidence-backed drug
membership, exclusion of explicit negative drugs, correct matched-control and
stratum handling, and reproducible downstream publication. Existing releases
remain authoritative until this successor is implemented and verified.

## Evidence and alternatives

- ADR 0071 and the released catalogs supply the current DEG counts and six-term
  fusion contract. A per-AF exception was rejected in favor of a uniform gate.
- The current `bases/cardiotox.py` implements six class-level source tables across
  three datasets. Merely renaming their categories would preserve the old unit
  of evidence and was rejected.
- GSE244740's deposited series summary explicitly distinguishes 33 cardiotoxicants
  from 9 non-cardiotoxicants. GSE217421's summary describes a mixed drug panel,
  not a universal positive-toxicity label. Thus pooling every structurally resolved
  panel drug was rejected. These are local deposited metadata observations;
  a fresh external drug-evidence review remains implementation work.
- Source of decisions: user confirmations in the 2026-11-16 grill-with-docs
  discussion covering the six numbered questions; no new DEG analysis was run.

## Accepted measured release 2026 11 16

The independently replayed #217 candidate measures nHF 1,661/1,661, ICM
1,647/1,647, AF 1/1, genetic_cardiomyopathy 1,482/1,482 and TOX 9,589/2,000
(axis-resolved DEG union / DEG2000). Four diseases pass the gate; AF retains
SLC9A7 but contributes no cross-disease rank or f2000 term. TOX has twelve
independent dataset sources. The axes contain exactly 2,000 and 5,000 ascending
frozen indices, are nested, and need zero tenv-HVG fallback slots. The coherent #219 release preserves these measured results and replaces live
mechanism-class products. Fresh source and integrated replays are byte-identical;
see `stage2/01_bases/log/acceptance.md` for acceptance evidence. Earlier design
expectations and historical ADR 0071 counts above remain as decision context.
