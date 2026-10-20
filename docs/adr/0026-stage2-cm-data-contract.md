# Stage 2 CM data — frozen corpus and filtered views

**Status:** accepted (2026-10-20)

Stage 2 needs cardiomyocyte perturbation data for model adaptation, but Stage 1's
`Union finetune data` is unusable for it: 978 landmarks, per-plate quantile
alignment to the L1000 pretrain, HG/`drug_idx` filtering, and a single collapsed
`cardiomyocyte` cell identity are all XPert-specific. `stage2/02_cm_data/` therefore
builds a model-agnostic corpus from the same four sources (`dtoxs_primary`,
`dtoxs_ipsc`, `gse244740`, `gse262419`) in **six decisions**:

1. **Gene space = four-source intersection, latest HGNC symbol.** The four
   sources use symbols (38,478), Ensembl IDs (44,154), and `SYMBOL_entrez` rows
   (22,537 → 19,703 symbols), so a shared space requires an identity rule.
   Genes are matched alias/previous-symbol-aware and Entrez-verified where
   available; unresolvable or ambiguous genes are dropped and reported. The
   intersection is capped by the TempO-Seq panel (~18.7k), so `gse262419` gates
   the universe. Chosen over per-file native spaces because a frozen shared
   space is what makes cross-view and cross-dataset comparison reproducible
   without re-deciding gene identity at every adapter.
2. **`X` = raw counts**, duplicate probes summed by symbol before any
   transform; one `log1p_cp10k` convenience layer. No quantile alignment, batch
   correction, or HVG selection — those are model-specific and would make the
   pharmacology-tier comparison depend on a preprocessing artifact.
3. **Controls are rows**, with `control_group = dataset::plate` and a `vehicle`
   column, plus per-group means in `uns`. CPA-style models need a control
   population; mean-based models need the mean; one representation cannot serve
   both without discarding information.
4. **Drug identity = desalted parent InChIKey, stereochemistry retained**
   (`drug_key`). Connectivity-only keys are unusable: doxorubicin and epirubicin
   are 4′-epimers and share the same InChIKey connectivity block. Salt and
   free-base entries stay separate rows but share identity, so `Tamoxifen` /
   `Tamoxifen citrate` and `Estradiol` / `17beta-Estradiol` cannot leak across a
   cold-drug split. Structure wins over name: the DToxS `AMI` arm carries a
   pentyl-benzofuran homolog (CID 16218975), not amiodarone, so it enters
   `all-smiles` but not `atc-moa`.
5. **Three pharmacology tiers for `gse262419`**, built as filtered projections of
   one corpus: `atc-moa` (113 drug keys / 1,047 samples), `atc-moa-mesh`
   (252 / 2,310), `all-smiles` (460 / 4,206). MeSH is a separate tier because it
   admits herbicides, insecticides, carcinogens and disinfectants — chemical
   taxonomy, not pharmacology. The two curated drug panels use `atc-moa` with a
   cited **manual supplement** for real drugs missing from ChEMBL's mechanism
   table (currently SB-431542); the `gse262419` tiers stay strictly automatic,
   because that boundary *is* the comparison.
6. **The corpus is frozen without matrices.** `corpus/v1/` stores the gene
   set, the drug set, the per-source resolved gene sets (the drift baseline)
   and the exclusion report; no expression matrix is stored. Bundles re-read
   the sources and apply the frozen universes.
   The build **hard-fails** on any gene or drug outside the corpus, so a given
   corpus version always yields the same rows. Promoting a corpus version is a
   deliberate edit. Loaders are copied from
   `stage1/xpert-fda-miri-v3/03_cm_finetune_data/_datasets.py` (commit `c921f7d`)
   rather than imported, because `stage2/AGENTS.md` treats Stage 1 as
   non-dependency and Stage 2 must diverge on naming, dose provenance, and gene
   identity.

**Consequences.** `datasets-1/` ships five views: `dtoxs_cm.h5ad` (957 + 118
rows, 46 drugs), `gse244740.h5ad` (2,520 + 552, 42), and three `gse262419` views
(1,047/2,310/4,206 perturbed rows, 113/252/460 drugs, each + 272 controls). The
intersection drops ~12k unmapped `dtoxs` identifiers and ~1.4k TempO-Seq genes
that are not shared; a future model needing them must rebuild with a new corpus
version, not silently widen a view. Stage 1 artifacts remain untouched and serve only as
the well-set regression baseline for QA.
