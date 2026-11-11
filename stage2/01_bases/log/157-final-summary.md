# Stage 2 bases — final audit and accepted state (#157)

Ticket #157, parent [#154](https://github.com/Moke110/arc/issues/154). This is the
closing record of the Stage 2 bases migration: what the accepted science is, what
was verified for this ticket rather than reused, what was removed and why, and
what a reader must not conclude from these numbers.

Branch `codex/stage2-bases-157` from `dev` `7496743688a9872f8588ad122216b671c4bb45a8`.
Authorities: [ADR 0062](../../../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md),
[ADR 0061](../../../docs/adr/0061-stage2-shared-bases-and-indexed-corpora.md), the
[`stage2/AGENTS.md`](../../AGENTS.md) step contract and
[`155-migration-map.json`](155-migration-map.json).

## 1. The accepted science

**One frozen gene axis.** Every accepted product is on the same 18,164-gene axis
in one frozen row order (`gene_index` 0–18163), byte-identical in
`output/gene-identity/genes.csv`, `output/gene-identity/tables/genes.csv` and
`output/shared-expression/genes.csv` (`5877a3cd…`). It is an indexed axis, not a
symbol list: `GeneIdentity` resolves native Tahoe and tenv namespaces with
explicit ambiguity errors, `labels`/`native_ids` cross the human/native boundary,
and the full-axis HGNC round trip plus both source-specific FAM153B identities
are checked on every `--verify`. FAM153B is index 5091: tenv `ENSG00000182230`,
Tahoe `ENSG00000289731` / token 61482.

**Raw counts.** tenv: 4,510 wells (4,206 treated / 272 control / 32 technical
blank) over 19,613 native genes, digest `efd07632…`; every source probe of the 16
GSE262419 plates is replayed in independent verification. Tahoe: 67,018 projected
18,164-gene integer pseudobulk profiles over 14 plates with 1,400 controls;
3,388 source shards, 95,624,334 cells independently replayed. Neither release
carries normalized expression or significance.

**Shared expression and training corpora.** The bounded CP10K/log1p projection is
67,018 × 18,164 and 4,478 × 18,164 float32, replayed row by row with
`max_abs_error` 0.0. The current training release holds **65,318 Tahoe** and
**4,206 tenv** samples with individual treated wells, `layers/control` (mean of
independently normalized matched controls) and `layers/delta` (treated minus
control); tenv uses same-plate DMSO, excludes Media controls, and carries no truth
labels. Both source replays pass with `max_abs_error` 0.0.

**Measured truth.** `output/tenv-truth/` holds **1,377 eligible parent–dose
contrasts** over **4,206 treated wells** and **16,596 control memberships**
(mean ≈12.05 controls per contrast), fitted with a two-sided Wald DESeq2 test on
exactly the indexed axis, within-contrast BH at FDR 0.05, plate effects for the
**6 multi-plate** contrasts (1,371 others are `~ condition`). **273** contrasts are
a distinct *tested-empty* state, not "untestable". Three independent R refits were
accepted by `log/records/truth_complete_verify.json`; no refit ran for this
ticket. The exact successor envelope `output/tenv-truth-successor/` re-authorizes
the *current* training release against those statistics: integer contrast counts,
sample/control identity, the statistical recipe and bit-exact
`X`/`control`/`delta` layers all match, with no new fit.

**HVG and partitions.** `output/tenv-hvg/` publishes four complete 18,164-gene
rankings — one full-tenv descriptive ranking plus three **train-only** rankings —
and deterministic parent-hash cohorts (`bases-cold-parent-v1`, seed in the
manifest). All 4,206 wells appear in each of the three folds; per fold the split
is train **2,508 / 2,508 / 2,535**, validation 297 / 282 / 282 and test
1,401 / 1,416 / 1,389. Full-tenv descriptive ranks are **not** model inputs.

**The report panel.** `figure/diagnostic_tenv_pr.png` (`35d14266…`) is the single
figure. Its complete plotted values are `figure/data/diagnostic_tenv_pr.csv`
(`ff3a82e8…`, 7,222 curve rows, 9 summary rows, 3 panels, 3 series) and the
CSV-only renderer reproduces the PNG byte for byte. The displayed quantity is the
**exact equal-group macro AP** of the measured delta layer, not the area under
the sampled fixed-group pooled curve:

| fold | AUPRC (combined) | AUPRC positive | AUPRC negative | eligible groups (c / + / −) |
|---|---:|---:|---:|---|
| 0 | **0.07981072489453635** | **0.3349457015461418** | **0.13737174931309615** | 369 / 107 / 366 |
| 1 | 0.08780864619894788 | 0.3673685095005024 | 0.14057355398255195 | 380 / 110 / 373 |
| 2 | 0.09289711674818489 | 0.3522541054443466 | 0.1452336179429624 | 355 / 114 / 349 |

Fold 0 evaluates 459 groups; the evaluated groups digest is
`7a5431492d50f029bc2830913afa18a1dd3c55e1b52353b5873746db3b27ab77`. These are
**measurement-derived diagnostic QA**, not trained-model performance: no model
prediction survives in this scope, so the panel asks whether the measured delta
layer ranks its own held-out groups sensibly, and nothing more.

## 2. Accepted products and their stable paths

| path | role | binding |
|---|---|---|
| `output/api/bases/` | the single public reader/metric implementation | tracked code; every build/audit/figure/test command imports it (#156 D1) |
| `output/gene-identity/` | indexed identity + frozen `tables/` | manifest `5f7517c3…`, genes `5877a3cd…` |
| `output/tenv-raw-counts/` | tenv integer counts, probe evidence, seal | `efd07632…` |
| `output/tahoe-raw-counts/` | projected Tahoe integer counts, token evidence, seal | `7ea21f45…` |
| `output/shared-expression/` | verified CP10K/log1p projection (reconstruction intermediate) | `fb352f2b…` |
| `output/training-corpora/` | current separate Tahoe and tenv training releases | `aac36d2d…`, Tahoe `fa6e5a4c…`, tenv `bab4140c…` |
| `output/historical-training-release/` | historical v1 truth-fit release — tenv h5ad + seal + exclusion CSVs only after #157 | manifest `3543df76…`, tenv `a95bba66…` |
| `output/tenv-truth/` | accepted 1,377-contrast measured truth | `ab73b301…` |
| `output/tenv-truth-successor/` | accepted exact-equivalence envelope | `2798e3b8…`, worker `65db039a…` |
| `output/tenv-hvg/` | descriptive + three train-only rankings and partitions | `953103a1…`, cohorts `c1dcce79…` |
| `figure/` | panel PNG, complete plotted values, CSV-only renderer | `35d14266…` / `ff3a82e8…` |

`output/` now contains exactly these ten entries (nine products plus the tracked
`output/api`); nothing else is a formal product.

### The truth / successor split, stated plainly

Two files look like they answer the same question and do not.

- `output/tenv-truth/` is the **fit**: the statistics themselves, produced from
  the historical v1 tenv release and sealed against it.
- `output/tenv-truth-successor/` is a **proof, not a fit**: an exact-equivalence
  envelope that says "the labels in the current `output/training-corpora/tenv.h5ad`
  are the same experiment as the ones these statistics were computed from". It
  binds absolute paths and is re-issued after a path move, never edited.

A reader who wants truth opens `output/tenv-truth/` through the partition reader.
A reader who wants to know whether the current training release may use those
labels checks `output/tenv-truth-successor/`. Neither file authorizes a CLI
default; `--successor` does.

## 3. Limitations

1. **Historical identity vintage (deferred decision).** The sealed identity
   product keeps the historical v1 mapping strings. **229** Tahoe
   `native_symbol` values differ from the current shared-expression annotation
   with **0** token-id differences and **no axis change**. Re-sealing would change
   a public evidence interface for a cosmetic vintage, so #156 D2 deferred it to
   an explicit user decision; it is still deferred and still documented, not
   silently repaired.
2. **The successor envelope binds absolute paths.** It records canonical
   `/home/chang/arc/...` paths, so it is valid for this checkout and this path
   layout. Moving the payload to another root requires a re-issue
   (`execution/reissue_truth_successor.py`), which is the documented behaviour,
   not a bug.
3. **Metric tolerances.** The successor proof and the diagnostic replay are
   deterministic and are compared exactly; fold-level macro AP is reported as the
   exact rational value the code computes. **No stochastic tolerance is claimed**
   because no stochastic step is in scope. File-integrity hashes are recorded
   separately from metric values and never stand in for them.
4. **Retained historical evidence.** `log/historical/seals/**`,
   `log/historical/shared-v1-mappings/**`, `log/diagnostic/**`,
   `log/records/**` and `log/health/**` keep their original bytes and their
   original path strings. One dependency is a documented `EXPECTED_GAP`:
   `output/historical-training-release/manifest.json → inputs.shared_manifest`
   (`8ac09bcb…`) seals a shared-v1 manifest that was deleted before this ticket.
   It is a historical gap, and it is reported as one rather than papered over.
5. **The rejected historical v1 Tahoe release** is not repaired, not relabelled
   and not claimed equivalent to the current release. Its payload was removed by
   #157; the rejection evidence and the compact counterexample remain
   (see §6).
6. **Retained by choice.** `output/shared-expression/` stays because the training
   manifest seals it and gate D4 (a preserved seal copy) is not satisfied.
   `log/diagnostic/cohorts.csv` stays as accepted audit evidence even though it
   duplicates `output/tenv-hvg/cohorts.csv`.
7. **No model claim.** Nothing here measures a trained model. The panel is QA of
   the measured delta layer; the HVG rankings are inputs for later modelling.

## 4. Provenance and reconstruction inputs

The science is reproducible from immutable inputs, never from the accepted
payloads alone. External `source/` trees and Stage 1 are untouched
(`git diff --name-only a53c9ff codex/stage2-bases-157 -- pert_data stage1` → 0
files; `script/check_paths.sh` → 0 unresolved literals).

| product | reconstruct from |
|---|---|
| identity | frozen axis/crosswalk in the product + pinned HGNC `hgnc_complete_set.txt` (`12363916…`, external) + the surviving historical v1 mapping CSVs |
| tenv raw counts | 16 immutable GSE262419 plate files |
| Tahoe raw counts | immutable Tahoe source (3,388 shards) |
| shared expression | both accepted raw-count releases |
| training corpora | both raw releases + the shared release |
| truth | the historical v1 tenv h5ad + tenv raw counts + the sealed R DESeq2 worker |
| successor | truth statistics + current training release + the same worker |
| HVG/partitions | `output/training-corpora/tenv.h5ad` |
| figure | `figure/data/diagnostic_tenv_pr.csv` only |

Supported dependency order (each heavy command refused existing destinations and
ran one at a time under the repository caps):

```bash
python stage2/01_bases/execution/build_identity.py --verify
python stage2/01_bases/execution/build_native.py --verify
python stage2/01_bases/execution/prepare_tahoe.py verify
python stage2/01_bases/execution/build_shared.py --verify
python stage2/01_bases/execution/build_corpora.py --verify
python stage2/01_bases/execution/verify_training_source.py --source tenv
python stage2/01_bases/execution/verify_training_source.py --source tahoe
python stage2/01_bases/execution/reissue_truth_successor.py   # only if the payload root moves
python stage2/01_bases/figure/script/render_diagnostic_tenv_pr.py \
    --csv stage2/01_bases/figure/data/diagnostic_tenv_pr.csv \
    --png stage2/01_bases/figure/diagnostic_tenv_pr.png
```

## 5. Migration map

The old→new layout is one committed, idempotent rename:
[`migrate_products.py`](../execution/migrate_products.py), recorded with
per-file inode/size/mtime/sha256 evidence in
[`156-migration-record.json`](156-migration-record.json) (20 moves, 1 `rmdir`, 0
problems) and bound to the plan in
[`155-migration-map.json`](155-migration-map.json). Reversal:

```bash
python3 stage2/01_bases/execution/migrate_products.py --rollback \
    stage2/01_bases/log/156-migration-record.json
```

Historical records were **not** rewritten to the new paths; the mapping lives in
that record and in `output/api/bases/identity.py`'s callers.

## 6. The removal record

Three payloads were removed after the #155 gates D1/D2/D3/D6 were satisfied and
the proof in [`157-removals.json`](157-removals.json) was recorded. Total
**15,274,607,787 B (15.27 GB)**; disk 869 G → 855 G used.

| id | removed | bytes | digest | why it could go |
|---|---|---:|---|---|
| R1 | `output/historical-training-release/tahoe_bulk.h5ad` | 11,201,616,617 | `6e3ac2ec…` | **rejected and superseded**: independent raw replay failed (`X row48`); the old aggregator dropped marker entries then indexed compacted arrays with original offsets, misattributing and dropping counts. Every accepted binding was enumerated (successor envelope 11/11, truth manifest, training manifests, identity) and **none** names it; the truth-fit code opens only `tenv.h5ad` and `manifest.json` |
| R2 | `execution/checkpoints/truth-tenv-building` (567 files) | 1,206,080,986 | tree `34d7afc0…` | a **superseded half-built** truth candidate (56 `parts/`, `identity.json`, `planned_*`, three `.npy` matrices). The accepted complete product `output/tenv-truth/` replaced it; the only reconstruction path is a full refit through the sealed R worker, and no resumption is claimed |
| R3 | `execution/checkpoints/tahoe-raw-counts` (28 files) | 2,825,100,184 | tree `715dd384…` | **rebuildable** per-plate resumption state for a build that is documented complete and source-accepted; no accepted consumer |

What survives from R1 is deliberate and checkable: the historical
`manifest.json` still declares `6e3ac2ec…`, `tahoe_bulk_exclusions.csv`
(`667ae1a8…`), `tenv.h5ad` (`a95bba66…`), `tenv_exclusions.csv` (`6359c415…`),
and `log/records/tahoe_training_v1_rejection.json` (`f40b14b9…`) plus
`tahoe_historical_raw_diagnosis.json` (`b4ce2c38…`). The regenerated audit check
now requires that absence **and** the matching removal-record entry **and** both
rejection records, so an unrecorded disappearance fails.

No empty old-path directory needed removing — the migration had already removed
`output/native-v1`, and every other old name was a rename. The only empty
directory left is `execution/checkpoints/`, the documented disposable-state
boundary, which holds no bytes.

## 7. Verification scope for this ticket — reused vs new

The distinction matters: re-reading an accepted record is not re-running the
science.

**Newly executed by #157**

- the regenerated metadata/structure audit (44 binding + 21 structural checks,
  `non_pass: []`, one documented `EXPECTED_GAP`), run twice from the main checkout
  with single-checkout roots and byte-identical both times
  (`1b5642cb…`) — the stale #156 record is replaced;
- the figure reproduction from the CSV alone (byte-identical `35d14266…`);
- the identity product reproduction through the corrected defaults
  (`--verify`; manifest and genes byte-identical) — this also fixed finding F1;
- 13/13 individual test files (46 tests, 1 documented opt-in skip);
- `check_paths.sh` (once: PASS, 0 unresolved), `check_tracked.sh` (PASS) and
  `verify_all.sh` (once: PASS 34/34, 4:48.32, 79.4 MiB peak);
- the removals themselves, their digests and the declared-vs-absent scan;
- the Stage 1 / external-tree immutability diff (0 files).

**Reused (record re-read, digest/binding re-checked, no re-execution)**

- tenv and Tahoe raw-count acceptance, the shared-expression replay, both training
  corpora and their source replays, the 1,377-contrast truth fit and its three R
  refits, the four HVG rankings and the held-out isolation result.

Machine-readable per-claim ledger: [`157-evidence-ledger.json`](157-evidence-ledger.json).
Gate-level detail, measured resources and deviations: [`157-gates.md`](157-gates.md).
Health logs: [`health/157_*.log`](health/).

`verify_all.sh` hashes the immutable **external** item manifests (34 items across
`pert_data`, `pert_model`, `disease_data`, `drug_emb_model`, `pheno_data`). It does
**not** hash Stage 2 expression payloads and is not a bases-product integrity
check; that is the metadata audit plus the #156 E1 payload re-hash.

## 8. Final tree state

`stage2/01_bases/` contains only `output/`, `figure/`, `log/`, `execution/`,
`AGENTS.md` and `README.md`, with no step-root duplicates and no old-name
directories anywhere below it. The tracked public implementation is
`output/api/bases/**`; everything else under `output/` is ignored payload with one
stable path per product.

## 9. Records

| record | what it holds |
|---|---|
| [`157-removals.json`](157-removals.json) | each removal: old path, size, sha256/tree digest, resolution proof, consumer scan, retained evidence, post-removal checks |
| [`157-gates.md`](157-gates.md) | acceptance gates A1–A8, health checks C1–C3, tests, findings F1–F5, deviations |
| [`157-evidence-ledger.json`](157-evidence-ledger.json) | per-claim reused-vs-new ledger with record digests |
| [`health/157_*.log`](health/) | check_paths, check_tracked, verify_all, tests and figure-reproduction measurements |
| [`156-metadata-binding-audit.json`](156-metadata-binding-audit.json) | regenerated single-checkout audit (A8) |
| [`155-migration-map.json`](155-migration-map.json), [`156-migration-record.json`](156-migration-record.json) | the plan and the executed old→new mapping |
| [`156-gates.md`](156-gates.md), [`156-decisions.md`](156-decisions.md) | #156's gates (including the pre-removal G14 state) and decisions D1–D3 |

`156-gates.md` G14 ("rejected Tahoe still present") and #156's G9 old-path claim
describe the tree as #156 left it. #157 supersedes them with A4 and the R1
removal; their bytes are unchanged, as historical records should be.
