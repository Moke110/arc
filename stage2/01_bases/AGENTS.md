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
