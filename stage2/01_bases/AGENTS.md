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
  verifies it, and the category reader checks it). The 15.7 MB of compact CSVs
  are force-tracked (`*.csv` is globally ignored); no large payload is kept.
- Tests: `execution/tests/test_drug_classes.py` (crosswalk resolution, the
  `(lookup_name, dataset_inchikey)` collision, the ≥3 rule, TOP2 mesh
  corroboration, antibody/unresolved exclusion, evidence-product determinism and
  relative-path manifests) and `test_disease_deg_cardiotox.py` (pooled-class
  arithmetic, multi-contrast collapse/rank, published-table schema); run each
  file individually. Inputs: `disease_data/GSE243674/**`,
  `pert_data/DToxS-iPSC/**`, `pert_data/GSE244740/**` and
  `pert_data/Annotation/**` are read-only — derived evidence is step-owned.

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
- Gene-axis decision evidence (#183): public implementation
  `output/api/bases/gene_axis.py`; build/replay `execution/build_gene_axis.py`,
  `--verify` included; `execution/gene_axis_sources.py` holds the two
  `pert_data/` cardiotoxicity adapters and `execution/build_mouse_orthologs.R`
  regenerates the pinned `output/gene-axis/mouse_human_orthologs.csv` (pinned
  `arc` R entry point, homologene + org.Mm.eg.db); the R script writes a candidate
  under `execution/gene_axis/`, checks its sha256 against the accepted table and
  replaces it only when byte-identical or with `--args --accept` (ADR 0062
  candidate-then-replace). Product is
  `output/gene-axis/` (`source_summary.csv`, `fused_rank.csv`, `overlap.csv`,
  `inflection.json`, `tenv_non_hvg_curve.csv`, `direction_summary.csv`,
  `manifest.json` + the ortholog table; large per-source intermediates stay
  untracked under `execution/gene_axis/`). Figure:
  `figure/gene_axis_source_summary.png`, plotted values
  `figure/data/gene_axis_source_summary.csv`, CSV-only renderer
  `figure/script/render_gene_axis_source_summary.py` (must keep working with
  `output/` unreadable). Records: `log/gene-axis/{record.md,inputs.json}`.
  Decision preparation only — it chooses no `N`, quota, direction rule or animal
  weight (#173), builds no corpus/model and opens no truth or Heliste file; the
  only tenv input is the published descriptive HVG ranking.
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
