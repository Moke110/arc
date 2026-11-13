# AGENTS.md — stage2/01_bases

> **Superseded layout, rebuild pending (ADR 0066).** This step's `output/` tree is being
> rebuilt as `raw_counts/` → `corpus/` → `view-2k/`/`view-5k/`, with `disease-deg/` and the
> released `gene_axis_2000.csv` / `gene_axis_5000.csv` in `gene-axis/`. The measured truth,
> the successor envelope, the diagnostic PR figure and `log/` are retired. Everything below
> describes the pre-0066 tree until the rebuild lands.

## Rebuild foundation (#186, ADR 0066)

- Public implementation: `output/api/bases/raw_counts.py` (raw-count release,
  readers, bridge), `output/api/bases/corpus.py` (corpus build/replay, matched
  controls, reader), `output/api/bases/release.py` (candidate-then-replace and
  digest-replay harness) and `output/api/bases/deg_rank.py` (the DEG-table/rank
  contract tickets #187-#190 import: `SymbolResolver`/`align_table`,
  `deg_effect`, `per_gene_deg`, `dataset_deg_rank`, `rank_scores`,
  `complete_rank`, `normalized_rank`). `bases.gene_axis` re-exports the moved
  functions; do not create a second live copy.
- Drivers: `execution/build_raw_counts.py` and `execution/build_corpus.py`, both
  `publish|verify`. Builds stage candidates under `execution/checkpoints/`,
  verify digests, then replace the accepted product; `verify` replays from the
  pinned inputs and compares every digest. Products: `output/raw_counts/`
  (`tenv.h5ad` native 4,510 x 19,613; `tenv_probe_mapping.csv`;
  `tenv_axis_mapping.csv` frozen-index bridge; `tahoe.h5ad` 67,018 x 18,164 CSR;
  `token_mapping.csv`; `manifest.json`) and `output/corpus/` (`tenv.h5ad`
  4,206 treated + 272 vehicle; `tahoe_bulk.h5ad` 65,318 treated + 1,400
  `DMSO_TF`; `exclusions.csv`; `manifest.json`). Counts are integer only; the
  corpus has controls as rows and `matched_control_group` correspondence.
- Provenance: `raw_counts` default mode is `source-accepted-replay` — the payloads
  are rehomed, byte-checked against `log/records/native_tenv_verify.json` and
  `log/records/tahoe_counts_source_acceptance.json`, and the Tahoe token mapping
  is re-derived from the immutable source metadata. `--full-rebuild` rebuilds both
  payloads from the external `source/` trees under `pert_data/` through
  `native.py`/`tahoe.py` and replays them from source; it is the stronger mode if the rehomed bytes are in doubt.
- Tests: `execution/tests/test_raw_counts.py`, `test_corpus.py`,
  `test_deg_rank.py`, hermetic fixtures in
  `execution/tests/release_fixtures.py`; run each file individually. The pre-0066
  raw/training products stay present read-only until the removal ticket (#192).

## Disease DEGs: nHF / ICM / AF (#187, ADR 0066)

- Public implementation: `output/api/bases/disease_deg.py` — the contract layer
  (`align_source_table`, `build_dataset_products`, `genotype_split`,
  `contrast_selection`, `validate_table`, `validate_rank`, the
  `DiseaseDeg` reader, category catalog/release-manifest assembly) plus
  `output/api/bases/disease_deg.R`, the **one live DEG implementation** for the
  disease side (limma: voom for integer counts, `lmFit` on the log-scale deposit;
  the pinned arc R env has limma but no edgeR, so voom's own library-size
  normalisation is used and TMM is deliberately not re-implemented).
- Sources: `disease_data/GSE270788`, `disease_data/GSE255612` and
  `disease_data/Reichart-CM` (ADR 0030 items; CM gate + per-donor pseudobulk +
  item-local DEG in each item's `script/`). They import the public module through
  `execution/public_api.py`; do not copy the pseudobulk or DEG logic into them.
- Driver: `execution/build_disease_deg.py publish|verify --all` (also
  `--category nHF|ICM|AF`), with `--root/--input-root/--sources` as test seams.
  Candidates stage in `execution/checkpoints/`; verify replays from the pinned
  item-local DEGs and compares every digest. Products:
  `output/disease-deg/{nhf,icm,af}/` — per dataset `<dataset>.deg.csv` +
  `<dataset>.rank.csv` (18,164-row frozen-axis rank frame) plus per-category
  `catalog.json` and `manifest.json`. There is deliberately no top-level
  catalog/manifest: `#190` aggregates per category.
- Source declarations (eligibility, tier, `shared_control_group`,
  `genotype_status`, `time_window_basis`, contrasts, n_case/n_control) are data
  in `execution/disease_deg_sources.py`; the published catalog carries them
  verbatim. Local input paths there are **item-relative**: the driver composes
  them from `source_item` at runtime, and the catalog pins tracked files as
  `inputs` and item-local derived CSVs as `derived_inputs` (relative name +
  sha256), so no repo-relative literal points at a gitignored `local/` payload
  and `script/check_paths.sh` stays green in a merged checkout.
- Tests: `execution/tests/test_disease_deg_human_sc.py` + `disease_deg_fixtures.py`
  (hermetic; schema, rank/tie-break, NA/unmeasured, genotype split, joint-design
  arm collapse, determinism, no absolute paths). Run each file individually.
- GSE270788 erratum: the deposit's arm sizes are Donor 7 / NICM 6 / ICM 5 / AMI 5
  (both label sources), not the ticket's 6/7; the deposit is used unchanged and
  `n_case`/`n_control` are 6/7, 5/7, 5/7. Recorded in issue #187 and the item
  author log; #192 carries it into ADR 0066.

## Isogenic CM models and `genetic_cardiomyopathy` (#188, ADR 0066)

- Public implementation: `output/api/bases/bulk_deg.py` is the single **derivation
  adapter** for the six in-vitro bulk ADR 0030 items (GEO series-matrix parsing,
  fail-closed sample-token binding, feature naming, the shared integerisation
  helper `integerise_counts` with its count-scale guard, and
  `derive_item_deg`), so each item's `script/build_deg.py` only declares its
  design and contrast arms. The live DEG implementation stays
  `output/api/bases/disease_deg.R`; the shared count contract is not relaxed.
- Sources: six new ADR 0030 items `disease_data/GSE224129|GSE234141|GSE304575|
  GSE305648|GSE253223|GSE194297` (two payloads each: the deposited matrix +
  the GEO series matrix; item-local `local/deg/` matrix/samples/full DEG table
  with a tracked `_deg.json` summary) and `disease_data/Reichart-CM`'s new
  `script/build_genetic_deg.py` (Option A, additive: the two genetic contrasts
  are derived from #187's tracked pseudobulk + genotype export, read-only; the
  derived table and tracked summary live in that item's declared `local/deg/`).
- Driver: unchanged `execution/build_disease_deg.py publish|verify --category
  genetic_cardiomyopathy` (also in `--all`). Declarations live in
  `execution/disease_deg_sources.py`; the driver now carries an optional
  category-level `isogenic_exception` into the catalog and optional per-dataset
  `isogenic_exception` / `strata` records, and verification accepts a
  multi-contrast table (every produced label must be declared and every declared
  label must win at least one gene). Product:
  `output/disease-deg/genetic_cardiomyopathy/` — eight datasets + catalog +
  manifest, ~30 MB force-tracked (compact CSVs).
- Tests: `execution/tests/test_disease_deg_isogenic.py` (35 tests) — the
  integerisation contract and its guard, as-deposited count validation, series
  matrix/sample binding, the declared isogenic control arm of every source, the
  multi-arm strata decision recorded in item summary and catalog and applied in
  the published table, catalog/manifest schema and digest replay, rank semantics
  and multi-contrast collapse, determinism of a driver replay, no absolute
  paths, and the genotype split (positive/ARVC donors only in this category,
  never nHF). Run it individually under the `arc` env.

## Drug-cardiotoxicity class evidence (#189, ADR 0066)

- Public implementation: `output/api/bases/cardiotox.py` (pinned code→drug→
  mechanism crosswalk, ATC level-4 class rule, TOP2 mesh corroboration, the
  ≥3-small-molecule rule, pooled class DEGs, product writers, publish/replay and
  the `DiseaseDeg`/`CategoryDeg` readers). The R compute steps are
  `output/api/bases/cardiotox_edgeR.R` (DToxS pinned pipeline on pooled wells)
  and `cardiotox_deseq2.R` (GSE243674 pooled DESeq2); both run under
  `/home/chang/miniconda3/envs/biomni_e1/bin/Rscript` (recorded as
  `conda:biomni_e1`), because the `arc` env's R lacks edgeR/DESeq2. edgeR >= 4
  dropped `decideTestsDGE`, which the pinned DToxS pipeline calls; the R driver
  supplies the documented two-group `decideTests(method = "separate")` shim.
- Driver: `execution/build_disease_deg_cardiotox.py publish|verify`; builds all
  three products in one candidate under `execution/checkpoints/`, verifies each
  manifest, then replaces `output/disease-deg/{cardiotox,cardiotox_TOP2,
  cardiotox_other}/`. `verify` rebuilds from the pinned inputs and compares every
  digest; never publish from a stale worktree without re-running it.
- Products (per-category `catalog.json` + `manifest.json`, no shared top-level
  JSON so tickets #187/#188 cannot collide): `cardiotox/` holds the crosswalk,
  `classes.csv`, `class_members.csv`, `excluded.csv`, `eligibility.csv`,
  `rejected.csv`; each category holds `<source>.<class_slug>.deg.csv` /
  `.rank.csv` with the ticket-01 DEG/rank contract; `catalog.json` uses
  schema `stage2-disease-deg-catalog-v1` with a `datasets` list of full
  records, each record keyed `method`, and a top-level `catalog_sha256`
  (payload digest with the field removed; `bases.cardiotox.catalog_digest`
  verifies it, and the category reader checks it). Each record also carries the
  real pooled arm sizes (`n_case`/`n_control`, derived from the pooled wells,
  with `n_case_basis`/`n_control_basis` and per-stratum counts in `contrasts[]`)
  — never a bare null. The 15.7 MB of compact CSVs are force-tracked (`*.csv` is
  globally ignored); no large payload is kept.
- Tests: `execution/tests/test_drug_classes.py` (crosswalk resolution, the
  `(lookup_name, dataset_inchikey)` collision, the ≥3 rule, TOP2 mesh
  corroboration, antibody/unresolved exclusion, evidence-product determinism and
  relative-path manifests) and `test_disease_deg_cardiotox.py` (pooled-class
  arithmetic, multi-contrast collapse/rank, published-table schema); run each
  file individually. Inputs: `disease_data/GSE243674/**`,
  `pert_data/DToxS-iPSC/**`, `pert_data/GSE244740/**` and
  `pert_data/Annotation/**` are read-only — derived evidence is step-owned.

## Released gene axes `G_2000` / `G_5000` (#190, ADR 0066)

- Public implementation: `output/api/bases/gene_axis.py` — the pure two-level
  fusion (`winning_tier`, `category_rank`, `fusion_rank_sum`, `candidate_table`,
  `axis_fill`, `axis_intersection`, `axis_composition`, `pn_distribution`,
  `rank_sum_histogram`, `kneedle_knee`, `second_difference_peak`,
  `detect_inflection`, `inflection_stability`) plus the digest-checked
  `GeneAxis` reader (`manifest`, `axis_sizes`, `candidate_table`, `axis`,
  `axis_frame`, `composition`, `digests`). It keeps the #186 re-export contract:
  the `bases.deg_rank` functions (`SymbolResolver`, `complete_rank`,
  `deg_effect`, `per_gene_deg`, `rank_scores`, `ALIGNMENT_CATEGORIES`,
  `PADJ_FLOOR`, `PADJ_THRESHOLD`, `normalized_rank`, `AXIS_SIZE`) are the same
  objects here — a test asserts identity, so a rewrite cannot drop them.
- Driver: `execution/build_gene_axis.py publish|verify`. Reads only
  `output/disease-deg/**` (through `bases.disease_deg.DiseaseDeg` for the
  catalog and digest checks, with the rank CSVs re-validated against the
  declared `n_measured`), `output/tenv-hvg/descriptive.csv` (never recomputed)
  and `output/gene-identity/genes.csv`. Candidates stage under
  `execution/checkpoints/`; only this ticket's products are replaced (the
  retired pre-0066 `gene-axis/` file set stays byte-identically in place for
  #192's removal); `--verify` replays from the pinned inputs and compares the
  manifest and every digest.
- Contract decisions implemented, with ADR 0066 as the basis: **winning tier**
  is the category's minimum catalog `tier`, plus records whose
  `isogenic_exception.applies` is true (ADR 0066 *Priority* and *Isogenic
  exception*; `genetic_cardiomyopathy` therefore fuses 2 tier-1 human
  CM-pseudobulk datasets + 6 tier-2 isogenic models); the **seventh source** is
  the published `descriptive.csv` `rank` column normalized as
  `rank / frozen-axis size` (ADR 0066 "each normalized by the frozen-axis
  size"; GLOSSARY *Source rank*). Level 1 is `sum(rank_norm)` over the
  winning-tier datasets (ties by `gene_index`); level 2 is the seven-source
  `rank_sum`; `pN` is the count of sources with `rank <= N`; `n_sources_measured`
  counts sources that measure the gene.
- Products (`output/gene-axis/`, force-tracked compact CSVs): `candidate_rank_sum.csv`
  (18,164 rows: `symbol`, `gene_index`, `tenv_hvg_rank`, six `<category>_rank`,
  `rank_sum`, `n_sources_measured`, `p2000`, `p5000`), `gene_axis_2000.csv` /
  `gene_axis_5000.csv` (ascending original `gene_index` + `symbol`; never
  renumbered), `axis_composition.csv` (per-category and `pN` distributions of
  both axes, so `pN` is not read as a continuous score) and `manifest.json`
  (fill rule, sizes, intersection, per-category composition and winning-tier
  decisions, 33 input digests, 5 output digests, environment; no absolute
  paths).
- Evidence panel: `figure/gene_axis_rank_sum_histogram.png`, complete plotted
  values `figure/data/gene_axis_rank_sum_histogram.csv`, CSV-only renderer
  `figure/script/render_gene_axis_rank_sum_histogram.py` (palette colours from
  the public `bases.palette`; reproduces the PNG byte for byte with `output/`
  unreadable). The Kneedle + second-difference inflection is a diagnostic in
  the manifest (parameters, candidates, seven leave-one-source-out replays)
  and may honestly be `false`; it never chooses `N`.
- Tests: `execution/tests/test_gene_axis.py` (38 tests, hermetic 8-gene
  fixture) covers winning tier/isogenic exception, level-1 identity and
  level-1/2 sums with `gene_index` ties, `pN` 0/7 boundaries, fill rule
  `0.5*N + 0.5*N` and its ordering, non-nested axes + intersection, ascending
  frozen-space axes, histogram, inflection detection/stability, candidate-then-
  replace determinism and tamper detection, manifest relative-path hygiene,
  figure CSV-only replay, `deg_rank` re-export identity and retired-file-set
  preservation. Run it individually under the `arc` env.

## Model-facing views `view-2k/` / `view-5k/` (#191, ADR 0066)

- Public implementation: `output/api/bases/view.py` — the one live view
  definition (`library_totals`, `log1p_cp10k`, `build`, `verify`,
  `publish_release`/`replay_release`) plus the bounded public reader
  (`View`, `ViewSource`, `open_view`, `View.digests`, `.library_totals()`,
  `.control_groups()`, `.matched_control_rows()`). It imports the corpus reader
  and the digest-checked `GeneAxis` reader; it never re-derives the corpus axis,
  the axis fill rule or the matched-control correspondence.
- Definition: normalize every corpus profile on the **full 18,164-gene library
  total** (`CP10K`, natural `log1p`), subset to the released axis's genes in
  ascending original `gene_index`, `control` = mean of the row's matched
  vehicle rows normalized the same way, `Δ = X − control` exactly. The view **never**
  re-normalizes on its own gene subset — one normalization, two gene subsets
  (ADR 0066 *Corpus and view semantics*). No `{N}k_eval_truth/` is created.
- Measured basis (released axes, this build): tenv axis = 12.93% (`G_2000`) /
  37.21% (`G_5000`) of full-axis CP10K (pooled all rows; treated-only 12.96% /
  37.16%); `tahoe_bulk` 10.46% / 24.78%. ADR 0066's 12.9% / 36.3% are
  design-record figures from an interim pre-#190 axis; the 2,000 figure
  reproduces and the manifest records the measured released-axis values as
  authoritative (the ~2.88× denominator ratio leaves the full-axis decision
  unaffected). Per-source full-axis library totals live in
  `manifest.sources.<source>.full_axis_library_total` and in
  `library_totals.csv`.
- Driver: `execution/build_view.py publish|verify` (one view at a time,
  streaming 256-row corpus batches; candidates stage in
  `execution/checkpoints/`; `verify` replays both views from the pinned corpus
  and gene-axis manifests and compares every digest). Products:
  `output/view-2k/` and `output/view-5k/`, each `tenv.h5ad` + `tahoe_bulk.h5ad`
  (`X` float32, `layers/control`, `layers/delta`; `obs` carries corpus row
  identity + `library_total`; `var` carries ascending original `gene_index` +
  `symbol`), `library_totals.csv` (`source`, `profile_id`, `role`,
  `library_total`; 71,196 rows) and `manifest.json` (recipe, normalization,
  input/output digests, counts, basis, environment; no absolute paths). The
  manifests are force-tracked; the ~4.7 GB float32 h5ads and
  `library_totals.csv` are gitignored and registered in
  `docs/data-inventory.md`. Measured cost: publish 3 m 35 s / 445,512 KiB peak
  RSS, verify 3 m 29 s / 443,036 KiB.
- Consumer contract: `02_models_tenv` still reads the retired
  `output/training-corpora/tenv.h5ad`; re-pointing it and re-running its accepted
  matrix are deferred (ADR 0066 *Deferred*). Its target is
  `output/view-{N}k/tenv.h5ad` `X`/`control`/`delta`; the view keeps both corpus
  roles as rows (`obs.role`), so the model selects `role == 'treated'` (the
  retired training h5ad held treated rows only). The two axes share one scale,
  so an axis comparison must declare one hyperparameter protocol across axes or
  a per-axis tuning budget.
- Tests: `execution/tests/test_view.py` (17 tests) + `view_fixtures.py`
  (hermetic corpus h5ads and gene-axis release written directly): the fixture
  where view-local and full-axis denominators differ and the full-axis result is
  asserted, control = mean of the named group, `Δ = X − control` exactly,
  replay reproduces `Δ`, exactly 2,000/5,000 ascending frozen genes (no
  renumbering), control scope per group, no `{N}k_eval_truth/`, determinism,
  tamper detection, manifest relative-path hygiene. Run it individually under
  the `arc` env.

Routing for the standalone Stage 2 foundation; commands and contracts are in
[`README.md`](README.md). ADR 0062 fixes the layout; ADR 0061 defines the frozen
gene and evaluation semantics. Every path below is a stable semantic path.

- Public implementation: `output/api/bases/` is the **single live** public
  reader/metric implementation (ADR 0062, ticket #156 decision D1). Build, audit,
  figure-provenance and test commands bind it with `import public_api` from
  `execution/`; never copy a module out of `output/api` to create a second live
  implementation. One deliberate exception is the accepted successor envelope
  `output/tenv-truth-successor/deseq2_contrasts.R`: it is a byte-identical,
  digest-sealed **provenance copy** of `output/api/bases/deseq2_contrasts.R` that
  `truth_successor.py` re-hashes to prove which worker produced the envelope.
  Keep it — it is sealed evidence, not an alleged duplicate.
- Identity: `execution/build_identity.py` and `output/gene-identity/`
  (`genes.csv`, `manifest.json`, and the frozen `tables/`). Historical v1 mapping
  evidence is `log/historical/shared-v1-mappings/`.
- Raw preparation: `execution/build_native.py`, `output/api/bases/native.py`,
  `output/api/bases/raw.py`; products are `output/tenv-raw-counts/` and
  `output/tahoe-raw-counts/`. Read the repository execution-environment rules
  before data access. External `source/` inputs are immutable.
- Shared expression: `execution/build_shared.py` and `output/api/bases/raw.py`;
  raw files are bound independently from the identity tables and the product is
  `output/shared-expression/`.
- Training corpora: `execution/build_corpora.py`,
  `execution/verify_training_source.py` and `output/api/bases/corpora.py`;
  product is `output/training-corpora/`; bounded raw-count replay verifies all
  three layers. `output/historical-training-release/` is the historical v1
  truth-fit release, not a training destination.
- Evaluator-only truth: `execution/build_significance.py` and
  `output/api/bases/significance.py`; product is `output/tenv-truth/`. The exact
  successor binding is `output/tenv-truth-successor/`, re-issued after a path
  move by `execution/reissue_truth_successor.py` and never edited. Ordinary
  training readers open only the h5ad.
- HVG/metrics: `execution/hvg_metrics.md` routes the train-only selector,
  `build_hvg.py`, `evaluate_predictions.py` and `replay_pr.py`; product is
  `output/tenv-hvg/`.
- Gene-axis products: the **released** axes `G_2000`/`G_5000` are documented
  in *Released gene axes (#190, ADR 0066)* above. The pre-0066 #183 file set
  (`source_summary.csv`, `fused_rank.csv`, `overlap.csv`, `inflection.json`,
  `tenv_non_hvg_curve.csv`, `direction_summary.csv`,
  `mouse_human_orthologs.csv`) and its figure/adapters
  (`figure/gene_axis_source_summary.*`, `execution/gene_axis_sources.py`,
  `execution/build_mouse_orthologs.R`) are **retired**, retained byte-identically
  only until the removal ticket (#192) deletes them; nothing reads them.
  Superseded decision-preparation semantics (direction composition, unweighted
  quota-free fusion, animal half-weighting, tenv non-HVG curve, ortholog table)
  must not be resurrected.
- Figure: `figure/diagnostic_tenv_pr.png`, complete plotted values in
  `figure/data/diagnostic_tenv_pr.csv`, and the CSV-only renderer
  `figure/script/render_diagnostic_tenv_pr.py`. The renderer must keep working
  with `output/` unreadable.
- Records: acceptance/execution records in `log/records/`, repo-health records in
  `log/health/`, dissolved and historical audit evidence in `log/historical/` and
  `log/diagnostic/`. Disposable build state lives in `execution/checkpoints/`,
  never inside `output/`; it is now empty — #157 removed the superseded truth-fit
  build state and the Tahoe per-plate cache, and `log/157-removals.json` holds
  their digests and resolution proof.
- Migration: `execution/migrate_products.py` is the one-shot, idempotent rename
  migration that produced this layout; its record is
  `log/156-migration-record.json` and `--rollback RECORD` reverses it.
- Heavy work: run one job at a time through `execution/capped.sh` (systemd user
  scope when available, otherwise an equivalent 10 GiB RLIMIT_AS) with
  `/usr/bin/time -v`, and record the measured peak.
