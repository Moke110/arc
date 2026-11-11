# AGENTS.md — stage2/01_bases

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
