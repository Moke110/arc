# #155 product, consumer and record inventory (2026-11-12)

Complete disposition inventory for `stage2/01_bases` at `dev` `a53c9ff`, read-only
from the main checkout. The machine record is
[`155-inventory.json`](155-inventory.json) (every file: size, tracked flag, role,
disposition, proposed stable path); this page is the reviewable summary.
Nothing was moved, rebuilt or deleted by #155.

Totals: 12 product directories, 36.10 GB (21.7 MB of it tracked), 120 tracked
execution/acceptance records (1.44 MB), 1 published figure and its plotted CSV at
the step root. No expression payload was opened while building this inventory.

## 1. Products

`kind` uses ADR 0062 / `stage2/AGENTS.md` vocabulary: *accepted reusable product*
(reused by another step), *reconstruction intermediate* (needed to rebuild an
accepted product), *plot data*, *audit evidence*, *disposable intermediate*,
*rejected/superseded*.

| directory | kind | files / size | seal (sha256) | tracked | proposed stable path |
|---|---|---:|---|---|---|
| `identity-v1/` | accepted reusable product | 2 / 0.9 MB | `5f7517c3…` | 2 | `output/gene-identity/` |
| `native-v1/` | accepted reusable product (tenv raw counts) + frozen axis/crosswalk tables | 5 / 95 MB | `e0941b24…` | 4 | `output/tenv-raw-counts/` + `output/gene-identity/tables/` |
| `tahoe-counts-v1/` | accepted reusable product (projected Tahoe raw counts) + build state | 31 / 5.78 GB | `7ea21f45…` | 0 | `output/tahoe-raw-counts/` |
| `shared-v1/` | historical audit evidence (mapping CSVs only) | 2 / 2.8 MB | — | 2 | `log/historical/shared-v1-mappings/` |
| `shared-v2/` | reconstruction intermediate (verified) | 8 / 5.20 GB | `fb352f2b…` | 0 | `output/shared-expression/` |
| `training-v1/` | historical accepted truth-fit input + rejected Tahoe payload | 5 / 11.69 GB | `3543df76…` | 3 | `output/historical-training-release/` |
| `training-v2/` | accepted reusable product | 5 / 11.51 GB | `aac36d2d…` | 0 | `output/training-corpora/` |
| `truth-tenv-v1/` | accepted reusable product (measured tenv truth) | 1384 / 0.60 GB | `ab73b301…` | 0 | `output/tenv-truth/` |
| `truth-tenv-v1.building/` | disposable intermediate | 567 / 1.21 GB | — | 0 | none — delete in #157 |
| `truth-tenv-v2/` | accepted reusable product (successor envelope) | 2 / 0.9 MB | `24b65aee…` | 1 | `output/tenv-truth-successor/` (re-issued, not moved) |
| `hvg-tenv-v1/` | accepted reusable product (HVG rankings + partitions) | 6 / 6.2 MB | `953103a1…` | 0 | `output/tenv-hvg/` |
| `hvg-metrics-diagnostic-v1/` | plot data + audit evidence | 12 / 2.5 MB | — | 0 | dissolved into `figure/` and `log/` |

Step-root `diagnostic_tenv_pr.png` (132 026 B, `35d14266…`) and
`diagnostic_tenv_pr.csv` (484 988 B, `eaad870b…`) are byte-identical copies of
`hvg-metrics-diagnostic-v1/diagnostic_pr.{png,csv}`; they are the published
figure and its plotted values, not a second product.

### Payload roles worth naming

- `native-v1/tenv.h5ad` — tenv raw integer counts, 4510 wells × 19613 native genes;
  digest `efd07632…` is the truth-fit input seal.
- `tahoe-counts-v1/tahoe.h5ad` — projected 18164-gene Tahoe integer pseudobulk
  counts, 67018 profiles; `token_mapping.csv` keeps the original token evidence;
  `.checkpoints/*.npz` (28 files, 2.85 GB) are build resumption state only.
- `shared-v2/*.npy` — CP10K/log1p normalized projection, 67018 × 18164 and
  4478 × 18164 float32.
- `training-v2/*.h5ad` — the current treated/control/delta corpora (65318 Tahoe,
  4206 tenv samples) plus explicit exclusion tables.
- `truth-tenv-v1/` — 1377-contrast `groups.csv`, `contrast_samples.csv`, three
  statistic matrices and 1377 per-contrast directional list payloads.
- `truth-tenv-v2/` — the successor envelope (`manifest.json`, 892 KB of
  per-contrast proof) and its sealed copy of the R DESeq2 worker.
- `hvg-tenv-v1/` — descriptive and three train-only rankings, plus the
  deterministic parent-hash partition assignments.

## 2. Consumers (every live reference)

Method: `git grep` over tracked files in the ticket worktree (untracked and
ignored trees excluded by construction); the external `source/` trees were never
scanned. Hits inside `execution/build_record/**` are historical evidence, not
runtime consumers, and are listed separately.

### Runtime consumers (code defaults and imports)

| consumer | reference | binds |
|---|---|---|
| `execution/bases/raw.py:60` | `output / "native-v1"` | tenv raw release (legacy default) |
| `execution/bases/raw.py:67` | `output / "tahoe-counts-v1/tahoe.h5ad"` | Tahoe raw release |
| `execution/build_identity.py:60,62,63` | `native-v1`, `identity-v1`, `shared-v1` | axis/crosswalk, identity output, mapping evidence |
| `execution/build_native.py:9` | `native-v1` | tenv raw output |
| `execution/build_shared.py:8` | `shared-v2` | shared-expression output |
| `execution/build_corpora.py:15,23` | `shared-v2`, `training-v2` | training build inputs/output |
| `execution/prepare_tahoe.py:25–27` | `native-v1/tables/{genes,gene_crosswalk}.csv`, `tahoe-counts-v1` | axis, crosswalk, Tahoe raw output |
| `execution/build_hvg.py:13,14` | `training-v2/tenv.h5ad`, `hvg-tenv-v1` | HVG input/output |
| `execution/verify_hvg_metrics.py:12–14` | `training-v2`, `hvg-tenv-v1`, `truth-tenv-v2` | HVG/diagnostic verification |
| `execution/build_significance.py:11,14` | `training-v1`, `truth-tenv-v1` | **truth-fit input** (historical release) and truth output |
| `execution/bind_truth_successor.py:12–19` | `truth-tenv-v1`, `training-v1`, `training-v2`, `native-v1/tenv.h5ad` (×2), `truth-tenv-v2` | successor envelope create/verify |
| `execution/verify_training_source.py:11` | `training-v2` | training source replay |
| `execution/tests/test_corpora.py:114` | `tahoe-counts-v1` (sibling of `native`) | fixture binding |
| `stage2/01_bases/README.md`, `AGENTS.md`, `execution/hvg_metrics.md`, `execution/truth_successor.md` | all of the above | operator commands |
| `docs/data-inventory.md:118–126` | `native-v1`, `tahoe-counts-v1`, `identity-v1`, `shared-v1`, `shared-v2`, `training-v1`, `shared-v2`/`training-v2`, `truth-tenv-v1` | size/provenance table (already stale: no HVG, diagnostic or successor row) |
| `.gitignore:143,146` | `output/` ignore + published-CSV exception | tracking policy |

### Data-level consumers (digest bindings read by code or acceptance)

- `training-v2/manifest.json → inputs.shared_manifest` = `shared-v2/manifest.json`
  digest `fb352f2b…`.
- `truth-tenv-v1/manifest.json → provenance` binds `native-v1/tenv.h5ad`
  (`efd07632…`), `training-v1/tenv.h5ad` (`a95bba66…`) and
  `training-v1/manifest.json` (`3543df76…`).
- `truth-tenv-v2/manifest.json → bindings` records **absolute paths** for the
  truth manifest, both training manifests, both training tenv files, both native
  files, the accepted verification record and the Rscript — and
  `bases/truth_successor.py:check_successor()` reads them at runtime.
- `hvg-tenv-v1/manifest.json → training_sha256` binds training-v2's `tenv.h5ad`.
- `hvg_metrics_actual_acceptance.json` binds the HVG manifest, truth manifest,
  successor manifest and training `tenv.h5ad` digests.

### Historical references (evidence records, not runtime)

`execution/build_record/**` contains ~60 path literals to the old directories
(launch records, resource records, acceptance records, failure diagnoses). ADR
0062 keeps those bytes and their path strings; the old→new mapping lives in
`155-migration-map.json`.

## 3. Acceptance and execution records

`execution/build_record/` — 120 tracked files, 1.44 MB, all tracked:

| role | count | examples |
|---|---:|---|
| acceptance record | 17 | `v2_actual_acceptance.json`, `truth_complete_verify.json`, `truth_successor_acceptance.json`, `hvg_metrics_actual_acceptance.json`, `training_v2_verify.json`, `shared_v2_verify.json`, `native_tenv_verify.json`, `truth_partition_acceptance.json` |
| build/replay record | 15 | `shared_v2_build.json`, `training_v2_build.json`, `training_v2_{tenv,tahoe}_replay.json`, `truth_complete_build.json`, `truth_successor_build.json`, `hvg_actual_build.json`, `tahoe_counts_smoke.json` |
| compact failure evidence | 6 | `tahoe_training_v1_rejection.json`, `tahoe_historical_raw_diagnosis.json`, `training_v1_acceptance_blocked.json`, `tahoe_training_v1_replay_failed.log`, `falsepass_diagnosis.json`, `truth_rebuild_assessment.md` |
| launch/service record | 8 | `truth_launch.json`, `truth_successor_launch.json`, `*_service.properties` |
| resource measurement | 13 | `truth_full_fit_resources.txt`, `truth_successor_build_resources.txt`, `hvg_actual_*_resources.txt` |
| contract/preflight/planning | 7 | `hvg_metrics_contracts.json`, `pre_migration_inventory.md`, `tahoe_counts_preflight.md` |
| publication record | 1 | `diagnostic_publication.json` |
| tracked evidence table | 1 | `tenv_parent_identity.csv` (81 explicit source-key conversions) and `tenv_training_source_replay.json` (reused v1 tenv replay) |
| repo-health record | 51 | `final_health/**` — `check_paths`, `check_tracked`, `verify_all`, per-item rechecks, seal ledger |
| execution log | 1 | `tahoe_counts_source_verify.log` |

Six acceptance records are authoritative for the current release chain; the
audit in [`155-audit.md`](155-audit.md) re-checks the digest bindings between them
and their products.

## 4. Reconstruction dependencies

| dependency | owner | why it must stay available |
|---|---|---|
| `pert_data/GSE262419/source/` (16 plates, hash table, compound annotation) | external item | tenv raw counts; digests sealed in `native-v1/tenv_manifest.json` |
| `pert_data/Tahoe-100M/source/` (3388 expression shards, metadata, `gene_metadata.parquet`) | external item | Tahoe projected raw counts; 3388 shard digests sealed in `tahoe-counts-v1/manifest.json` |
| `pert_data/Annotation/source/hgnc_complete_set.txt` | external item | axis/crosswalk/identity; digest sealed in the identity manifest |
| drug metadata table (digest `7a04c7a6…`) | bases input | desalted parent identity for both corpora |
| R 4.5.3 / DESeq2 1.50.2 at `/home/chang/miniconda3/envs/biomni_e1/bin/Rscript` (`4b08b95a…`) | workstation runtime | truth fit and successor recipe; versions sealed in the truth and successor manifests |
| RDKit 2026.03.2 | `arc` environment | parent identity (`Cleanup`, `FragmentParent`, `Uncharger`) |
| `shared-v1/` mapping CSVs | bases historical | the identity seal binds their digests; the Tahoe mapping describes the deleted 62710-gene matrix, so it is evidence, not a rebuild input |
| `training-v1/tenv.h5ad` + `training-v1/manifest.json` | bases historical | original truth-fit input seal and successor-envelope binding |
| per-plate `tahoe-counts-v1/.checkpoints/` | bases build state | resumption only; the release build completed |

## 5. Known omissions in this inventory

- Sizes are `stat()` sizes at 2026-11-12; they change while builds run.
- The inventory does not re-hash expression payloads (`*.h5ad`, `*.npy`, `*.npz`);
  it records the digest each seal declares. Re-verification is a #156 gate.
- 34 GB of Tahoe source shards and all other external trees are outside this
  step and outside this inventory.
