# Stage 3 starts with drug–phenotype data, independent of transcriptome pairing

**Status: accepted (2026-11-04)**

Stage 3 will predict the direction of change of one cardiomyocyte phenotype readout
under one specified condition from a drug-induced perturbome. The intended output
has three classes: increase (+1), no change (0), and decrease (−1). The input
representation (fold change, rank, direction, or a combination) is undecided.

The first effort plans the experimental drug–phenotype database: feasible sources,
inclusion criteria, and curation rules, ready for the repository's subsequent spec
workflow. Only drug–phenotype pairing is required. Transcriptomes will be supplied
by Stage 2 model predictions; generating, matching, or evaluating those predictions
is outside this effort. Requiring jointly measured transcriptomes would unnecessarily
exclude phenotype experiments and couple this data search to a separate model task.

The requested search covers:

- Drug exposure for 24–48 h followed by MTT/CCK8 viability measurement.
- Drug exposure for 24–48 h followed by hypoxia/reoxygenation and MTT/CCK8 measurement.
- Other cardiomyocyte readouts after drug exposure for 24–48 h, optionally followed
  by hypoxia/reoxygenation, with at least 500 independent drug–readout observations
  for one readout under comparable conditions. This minimum applies to this third
  category only; it is not a minimum of 500 distinct drugs.

Search all human iPSC-derived, human primary, animal primary, and cardiomyocyte-like
cell-line models, retaining species and cell-model strata. Inclusion in a shared
training set remains undecided. Report distinct drugs, drug × condition × readout
combinations, and repeats separately; repeat wells do not increase the independent
observation count. Cross-study aggregation is allowed with provenance retained;
the exact condition-equivalence and deduplication rules remain decision tickets.

Data-source research may establish that the requested threshold cannot be verified.
It must report that result rather than pool unrelated endpoints or relax the
requirements silently. Control choice, three-class thresholds, and treatment of
missing uncertainty remain open; non-significance is not an agreed definition of
the no-change class.

The explicitly invoked wayfinder workflow governs this planning effort: unresolved
decision issues and its map remain open across sessions, and research notes are
committed on `dev`. This overrides the generic close-all-issues run rule and the
skill's throwaway research-branch convention for this effort. No implementation
tickets or production data ingestion are authorized by completion of the map alone.

See [GLOSSARY.md](../../GLOSSARY.md), [roadmap](../roadmap.md),
[tracker operations](../agents/issue-tracker.md), and ADR
[0030](0030-external-data-provenance-and-step-layout.md). The GitHub wayfinder map
is the canonical index of subsequent decisions:
[Stage 3 决策地图：心肌细胞药物表型数据库](https://github.com/Moke110/arc/issues/21).
