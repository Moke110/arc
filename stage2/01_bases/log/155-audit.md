# #155 audit: evidence checked against actually bound metadata

Scope: `stage2/01_bases` products, seals and existing acceptance records at `dev`
`a53c9ff`, read-only from the main checkout. Method: metadata and
non-expression products only. No `*.h5ad`, `*.npy`, `*.npz` or `*.bin` was opened;
no `backed="r"`; no external `source/` tree was scanned; `verify_all.sh` was not
run. Raw output: [`155-metadata-binding-audit.json`](155-metadata-binding-audit.json).
Re-run: `python3 stage2/01_bases/log/155-reproduce.py --step-root <main checkout>/stage2/01_bases --write /tmp/audit.json`.

Measured: **0.42 s wall, 130 MiB peak RSS, 0 expression payloads opened.**
Result: **33 binding checks (32 PASS, 1 documented `EXPECTED_GAP`) and 13/13
structural checks PASS.** No product was modified.

## What "checked against bound metadata" means here

Two different claims are kept apart:

- **Reused evidence** — an acceptance record produced by the original build or
  audit run. #155 read it and checked that its recorded seals still agree with
  the files on disk. It did **not** re-run the scientific computation.
- **New checks** — digest recomputation of non-expression files, digest-to-digest
  cross-binding between records, and structural/scientific-invariant checks on
  the small tables. These are #155's own work.

Nothing below claims a re-fit, a re-replay or a full rescan.

## New checks (this ticket)

Binding checks, all PASS:

1. `identity-v1/manifest.json` → `genes.csv` (`5877a3cd…`), native
   `tables/genes.csv` (`19ce44ac…`), `tables/gene_crosswalk.csv` (`bc19c6a8…`),
   HGNC annotation (`12363916…`, external item), `shared-v1` mappings
   (`73cf97fc…`, `f1c3abfb…`).
2. `native-v1/tenv_manifest.json` → `tenv_probe_mapping.csv` (`f3e3dfcd…`);
   `tenv.h5ad` digest is declared (`efd07632…`) but not re-hashed (payload).
3. `shared-v2/manifest.json` → all six non-expression files; `genes.csv` and
   `gene_crosswalk.csv` input digests; `tahoe.h5ad`/`tenv.h5ad` input digests are
   declared, not re-hashed.
4. `tahoe-counts-v1/manifest.json` → `mapping_sha256` matches
   `token_mapping.csv`; 3388 source-shard digests declared (not re-hashed — the
   immutable external source is 337.6 GB).
5. `training-v1`/`training-v2` manifests → both exclusion tables.
6. `truth-tenv-v1/manifest.json` → 1377 per-contrast payloads, `groups.csv`,
   `contrast_samples.csv`, `genes.csv`; statistic matrices declared only.
7. `truth-tenv-v2/manifest.json` → `source_truth_files` (the v1 payload digests),
   `worker_sha256` = `deseq2_contrasts.R`.
8. `hvg-tenv-v1/manifest.json` → cohorts and all four rankings.
9. `hvg-metrics-diagnostic-v1/acceptance.json` → all ten diagnostic files.

Cross-binding checks, all PASS:

10. `training-v2.inputs.shared_manifest` = `shared-v2/manifest.json` digest.
11. `truth-v1.provenance.training_manifest_sha256` = `training-v1/manifest.json`.
12. `truth-v1.provenance.native_counts_sha256` = native `tenv.h5ad` digest;
    `tenv_training_sha256` = `training-v1/manifest.json → files.tenv.h5ad`.
13. `truth-v1.provenance.worker_sha256` = `truth-tenv-v2/deseq2_contrasts.R`.
14. `truth_complete_verify.manifest_sha256` = `truth-tenv-v1/manifest.json`;
    `truth_successor_acceptance.manifest_sha256` = `truth-tenv-v2/manifest.json`.
15. `hvg_metrics_actual_acceptance.json` → HVG manifest, truth manifest,
    successor manifest and training `tenv.h5ad` digests.
16. Step-root published `diagnostic_tenv_pr.{png,csv}` = byte-identical to the
    internal accepted artifacts.
17. All eleven `truth-tenv-v2` binding paths exist, are canonical absolute paths,
    and every non-expression binding digest matches.

Structural / scientific-invariant checks, all PASS:

18. Identity axis: 18164 rows, `gene_index` 0…18163, unique symbols and HGNC ids,
    index 5091 = FAM153B.
19. Frozen axis order identical between `native-v1/tables/genes.csv` and
    `identity-v1/genes.csv`; native tables and `shared-v2/genes.csv` are the same
    bytes.
20. Truth: 1377 contrast records, all `eligible`, unique ids; 4206 treated well
    memberships and 16596 matched-control memberships; 273 tested-empty
    contrasts; 6 multi-plate contrasts; 1377 per-contrast payloads covering
    exactly the table.
21. HVG partitions: 12618 rows over 4206 wells, each well in exactly three folds
    and exactly one `test` role; test/train/validation sizes 1401/2508/297,
    1416/2508/282, 1389/2535/282.
22. All four HVG rankings are complete permutations of 18164 genes.
23. Diagnostic panel CSV: three fold labels × three metrics, 7222 curve rows.

## Findings that change the plan

**F1 — the successor envelope is path-bound at runtime, so it must be re-issued,
not moved.** `bases/truth_successor.py:check_successor()` resolves every
`bindings[*].path`, requires the canonical absolute path, hashes it, and derives
the truth directory from `bindings.truth_manifest.path`. Renaming
`truth-tenv-v1`, `training-v1`, `training-v2`, `native-v1/tenv.h5ad` or the
accepted verification record silently breaks `SuccessorTruth`,
`evaluate_predictions --successor` and `verify_hvg_metrics`. Measured cost of the
re-issue: create 4:20 wall / 867 MiB peak, verify 4:21 / 870 MiB (single core).

**F2 — `native_column` in the Tahoe evidence is not an off-by-three defect.**
`tahoe-counts-v1/token_mapping.csv` records `native_column` = the immutable
`gene_metadata` row while `native_token_id` = that row's `token_id` column; the
external `gene_metadata.parquet` row 0 has `token_id` 3, so the two columns are
uniformly three apart. `bases/raw.py:tahoe_axis_mapping()` documents this.
Checked and dismissed; do not "repair" it.

**F3 — `shared-v1` and `shared-v2` Tahoe mappings are different evidence, not
duplicates.** Tenv mappings are byte-identical (`f1c3abfb…`). Tahoe mappings keep
identical token ids but differ in `native_column` (deleted 62710-gene matrix
column vs projected 18164 axis column) and in `native_symbol` for **229 genes**
(annotation vintage; e.g. `ACKR5` was `GPR182`; `ACTMAP` was `C19orf54`).
Consequence: re-sealing `identity-v1` against `shared-v2` would change 229
native-symbol evidence strings and the manifest digest — the axis bytes would not
change. #155 keeps the sealed identity product as accepted and records the
alternative as a deferred option (see the migration map).

**F4 — `groups.csv` has 1377 CSV records over 1414 physical lines.** Some fields
(`size_factors`, `messages`) contain newlines, so `wc -l` overstates the table.
The audit first hit this and was corrected; the manifest digest matches, and the
record count equals `group_count` 1377.

**F5 — `v2_actual_acceptance.json` carries a stale `pending` list.** It lists four
gates that later records closed (Tahoe source acceptance, Tahoe training-v1
replay, truth successor envelope, HVG/metrics acceptance). Read it as a
point-in-time record; the later records are authoritative.

**F6 — the published panel cannot yet be reproduced from its CSV alone.** The
displayed macro-AP annotation values live in
`hvg-metrics-diagnostic-v1/diagnostic_pr.json`, not in `diagnostic_pr.csv`
(columns: `recall, precision, contributing_groups, eligible_groups, label,
metric`). ADR 0062 requires the plotted CSV to carry the complete plotted values,
so #156 must add the annotation values (and fixed-group definitions) to
`figure/data/diagnostic_tenv_pr.csv` before the CSV-only renderer can reproduce
the PNG bytes. Replay evidence shows byte-identical PNG rendering is achievable
from the same values.

**F7 — byte-identical duplicate artifacts are intentional but must not multiply.**
`diagnostic_tenv_pr.{png,csv}` at the step root duplicate the internal artifacts;
`hvg-tenv-v1/cohorts.csv` duplicates the diagnostic copy; `shared-v2/genes.csv`
duplicates `native-v1/tables/genes.csv`. The migration keeps one copy per role.

**F8 — one documented historical gap.** `training-v1.inputs.shared_manifest`
(`8ac09bcb…`) seals a `shared-v1/manifest.json` that no longer exists: the
historical release's expression payload and seal were deleted. `training-v1` is
therefore a historical seal, not a reconstruction input. Recorded as
`EXPECTED_GAP`, not a defect to repair.

## Explicitly not done

- No expression payload was opened, sliced, re-hashed or replayed.
- No `verify_all.sh`; no external item hashing (~52 GB); no Tahoe source rescan
  (337.6 GB).
- No product was rebuilt, moved or deleted; no manifest was rewritten.
- Health checks are recorded separately in [`155-health-check.md`](155-health-check.md);
  the repo-wide PASS records reused here come from the earlier audit
  (`execution/build_record/final_health/`, commits `b17049c`/`aaa48b5`). #156 must
  re-run `check_paths.sh` and `check_tracked.sh` for the new layout.

## Reused evidence reused verbatim

`native_tenv_verify.json` (4510 source rows replayed), `shared_v2_verify.json`
(67018 / 4478 rows, `max_abs_error` 0.0), `training_v2_verify.json` (65318 rows /
1400 controls, 4206 / 192, `max_abs_error` 0.0), `training_v2_{tenv,tahoe}_replay.json`,
`v2_actual_acceptance.json`, `truth_complete_verify.json` (1377 contrasts, 3
independent R refits, 6 multi-plate, 273 tested-empty),
`truth_partition_acceptance.json` (1377 unique test groups, negative guards),
`truth_successor_acceptance.json` + `truth_successor_verify.json` (1377-contrast
equivalence proof, no new R fits), `hvg_actual_build.json` /
`hvg_actual_verify.json`, `hvg_metrics_actual_acceptance.json`,
`tahoe_counts_source_acceptance.json`, `tahoe_historical_raw_diagnosis.json`,
`tahoe_training_v1_rejection.json`, `diagnostic_publication.json`,
`final_health/summary.json`.

Full ledger with per-claim digests: [`155-evidence-ledger.json`](155-evidence-ledger.json).
