# #155 gates: reviewable acceptance list for #156 and #157

Every gate names an observable outcome and the artifact that proves it. Gates are
grouped by phase: #156 publishes the stable layout, #157 removes superseded
payloads. Nothing here authorises deleting or rebuilding a product by itself.

## A. #156 publication gates

| id | gate | proof |
|---|---|---|
| G1 | Every product sits at exactly one stable path from [`155-migration-map.json`](155-migration-map.json); no versioned product directory remains under `output/`. | `find stage2/01_bases/output -maxdepth 1 -type d` matches the map's `stable_path` set; `git status` shows only intended moves. |
| G2 | No expression payload changed bytes during the move. Re-hash the large payloads once, serially, under the caps (≈20.4 GB read, <1 GiB RSS, ~1–2 min). | `sha256sum` of `tenv.h5ad`, `tahoe.h5ad`, `tenv.h5ad`/`tahoe_bulk.h5ad`, both `*.npy`, the truth statistic matrices equals the map's declared digests. |
| G3 | Every non-expression digest re-verified with the reproducible audit. | `python3 stage2/01_bases/log/155-reproduce.py --step-root stage2/01_bases --write /tmp/audit-after.json` reports `non_pass: []`, one expected gap only, and the same manifest digest set. |
| G4 | The truth successor envelope is re-issued, never edited: bound paths are the new stable paths, create and verify both PASS, truth bytes untouched. | new `output/tenv-truth-successor/manifest.json` digest recorded in `log/`; `truth_successor_verify.json` re-issued with `groups_compared: 1377`; the superseded envelope bytes are in `log/historical/seals/`. |
| G5 | `SuccessorTruth` and the standalone evaluator work through the new paths only. | `evaluate_predictions.py --successor output/tenv-truth-successor --training output/training-corpora …` and `verify_hvg_metrics.py` run against the stable paths; no fallback to an old path exists in code. |
| G6 | Exactly one published figure triple: `figure/diagnostic_tenv_pr.png`, `figure/data/diagnostic_tenv_pr.csv`, `figure/script/render_diagnostic_tenv_pr.py`; the renderer reads only the CSV. | renderer run in a directory where `output/` is unreadable produces a PNG byte-identical to `35d14266…`; the CSV contains every plotted value including the displayed macro-AP annotations (finding F6). |
| G7 | Step-root duplicate `diagnostic_tenv_pr.{png,csv}` removed; publication record keeps the old hashes. | hashes recorded before removal; no file of those names at the step root afterwards. |
| G8 | Audit evidence moved out of `output/`: HVG group tables, panel summary, replay artifacts, acceptance and visual-QA records under `log/diagnostic/`; execution/acceptance records under `log/records/`; repo-health records under `log/health/`; superseded seals and mappings under `log/historical/`. | `find stage2/01_bases/log -type d` matches the map; no `*.h5ad`/`*.npy` under `log/`. |
| G9 | No public runtime consumer references an old path family; historical path strings survive only in tagged evidence records. | `git grep -nE 'identity-v1|native-v1|shared-v[12]|training-v[12]|truth-tenv-v[12]|tahoe-counts-v1|hvg-tenv-v1|hvg-metrics-diagnostic-v1' -- stage2/01_bases ':!stage2/01_bases/log/records' ':!stage2/01_bases/log/health' ':!stage2/01_bases/log/historical' ':!stage2/01_bases/log/155-*'` returns only historical record files. |
| G10 | README and AGENTS describe the accepted products, limitations and the historical/current boundary in the new vocabulary; the public implementation location matches decision D1. | review of `stage2/01_bases/README.md`, `AGENTS.md`, `execution/hvg_metrics.md`, `execution/truth_successor.md`. |
| G11 | `docs/data-inventory.md` lists every stable product with refreshed sizes and the reconstruction note; it is no longer missing HVG, diagnostic and successor rows. | review of that section. |
| G12 | Repository health checks re-run for the new layout: `script/check_paths.sh` and `script/check_tracked.sh`. Neither hashes expression payloads. `check_paths` must run in the checkout that owns the external payloads (the main checkout) — an isolated worktree cannot resolve external literals — and must be budgeted in minutes, not seconds, because its whitelist loop is subprocess-bound (724 literals × 53 whitelist lines ≈ 77k spawns; #155 measured >330 s without completing). `script/verify_all.sh` is **not** run. | both records with exit code, scan counts, wall time and peak RSS under `log/`; `check_tracked` must show zero holes. Existing scope and cost evidence: [`155-health-check.md`](155-health-check.md). |
| G13 | Existing scientific acceptances still bind: truth labels unchanged, HVG rankings unchanged, training/raw digests unchanged, diagnostic values unchanged. | the re-run audit (G3) plus the re-issued successor envelope (G4); no re-fit, no retrain. |
| G14 | The rejected Tahoe payload is still present and still marked rejected; nothing reads it. | `tahoe_training_v1_rejection.json` and `tahoe_historical_raw_diagnosis.json` unchanged; no code path opens `historical-training-release/tahoe_bulk.h5ad`. |

## B. #157 removal gates

| id | gate | proof |
|---|---|---|
| D1 | Superseded Tahoe training payload (`historical-training-release/tahoe_bulk.h5ad`, 11.20 GB, digest `6e3ac2ec…`) may be removed only after the successor envelope is re-issued and a removal record states the old path, digest, rejection record and the compact counterexample that replaces it. | removal record in `log/records/`; `tenv.h5ad` and the sealed historical manifest remain byte-identical. |
| D2 | `truth-tenv-v1.building/` (1.21 GB, 567 files) may be removed; refit through the sealed R worker is the only reconstruction path and must be stated as such. | removal record lists `parts/` layout, `identity.json`, `planned_*` and the matrices; no accepted product reads it. |
| D3 | `.checkpoints/` plate state may be removed once moved out of `output/` and the Tahoe release build is documented complete. | removal record; no resumption claim. |
| D4 | `shared-expression/` payload (5.20 GB) may be retired only if `log/historical/seals/shared-expression-manifest.json` preserves the sealed manifest bytes so training-v2's declared input seal stays checkable, and one documented command regenerates it from sealed raw counts. | seal copy present; regeneration command recorded; no consumer outside the training build. |
| D5 | `shared-v1/` mapping evidence stays as historical audit evidence unless the identity product is re-sealed by explicit decision; if re-sealed, the 229 Tahoe native-symbol deltas are reviewed first. | decision record; identity manifest digest change explained. |
| D6 | After removals, no accepted product's manifest declares an absent file except where a removal record exists, and every historical manifest keeps its original bytes. | rerun of the reproducible audit plus a declared-vs-absent list in the removal record. |

## C. Exact tenv truth-fit substitution prerequisites

The current accepted truth (`tenv-truth`) was fitted from the **historical**
training release's tenv file. The successor envelope is the only artifact that
authorises using those labels with the current training release. Substitution is
proved only when **all** of the following hold; any failure means the affected
truth must be refitted, never relabelled.

| id | prerequisite | where it is enforced |
|---|---|---|
| T1 | The original truth has a complete `stage2-parent-dose-truth-verification-v1` report: `complete: true`, `result: PASS`, `dataset: tenv`, `group_count` 1377, `gene_count` 18164, ≥3 distinct independent R refits covering ordinary, multi-plate and tested-empty classes. | `bases/truth_successor.py:create_successor/check_successor`, `accepted_refits()` |
| T2 | The verification report's `input_seals` equal the truth manifest's `provenance` fields `native_counts_sha256`, `tenv_training_sha256`, `training_manifest_sha256`, `worker_sha256`. | `accepted_input_seals()` |
| T3 | Axis identity: the same 18164 integers in the same order (`axis_sha256` = `sha256(tenv-truth/genes.csv)`). | `check_successor()` |
| T4 | Native counts: the successor native file's producer manifest seals its digest; the successor native file equals the original (currently the same file, `efd07632…`). | `producer_binding()`, exact integer comparison |
| T5 | Sample identity and covariates: identical `sample_id` sets, no duplicates, and equal `source_profile_id`, `drug_parent`, `dose_um`, and where present `source`, `group_id`, `cell_identity`, `time_h`, `vehicle`, `source_inchikey`, `source_parent_key`. | `sample_identity()` |
| T6 | Matched-control membership: identical sorted control-id sets per sample; controls joined by stable profile id, never by row position. | `sample_identity()` |
| T7 | Contrast integrity: identical `counts_sha256`, `sample_ids_sha256`, `source_ids_sha256`, `control_ids_sha256`, `design`, `n_treated`, `n_controls` for each of the 1377 contrasts; identical BH values, test statuses and ranked lists. | `check_successor()` contrast loop, `verify_truth()` replay |
| T8 | Expression equivalence: `X`, `control`, `delta` finite and **bit-exact** (`numpy.array_equal`, no tolerance; one float32 ULP fails), and each layer independently recomputed from raw counts with full-axis CP10K/log1p and independently normalized matched-control means. | `compare_layers()` |
| T9 | Recipe and runtime: exact statistical recipe dict (two-sided DESeq2 Wald, within-contrast BH over the whole 18164 axis, `independent_filtering: false`, `prediction_significance: false`, plate design rule, replication rule), R 4.5.3, DESeq2 1.50.2, worker sha256, current Rscript binary digest. | `check_successor()` |
| T10 | Training manifest normalization equals the recipe's normalization string; the envelope's recorded binding paths resolve and hash to the recorded digests. | `check_successor()` |

Corollaries recorded for #156:

- Matching raw paths or approximately equal normalized expression never
  substitutes; proportional integer-count changes require a new fit.
- Re-issuing the envelope after a path move does **not** weaken T1–T10: the proof
  is recomputed, and the truth statistic bytes are referenced, not copied.
- The original truth manifest is never rewritten; the successor envelope lives at
  its own stable path.

## D. Rejected Tahoe exclusions (never eligible inputs)

| payload | verdict | evidence |
|---|---|---|
| `training-v1/tahoe_bulk.h5ad` (11.20 GB, `6e3ac2ec…`) | **rejected and superseded**; not a training input, not a label source | `tahoe_training_v1_rejection.json`: independent raw replay fails; the old source aggregator removed marker entries and then indexed compacted arrays with original offsets, misattributing and dropping counts. |
| the historical Tahoe native derived matrix (62710 genes) | deleted, not reconstructed, not an input | `tahoe_historical_raw_diagnosis.json`; `shared-v1/tahoe_mapping.csv` describes it and is therefore historical evidence only. |
| any Tahoe significance labels or PR | do not exist and are not published | ADR 0061: only tenv has ranked significance truth; `prediction_significance: false`; no Tahoe label build exists. |
| the 1380 legacy parent–dose groups | never used | `truth_launch.json`/`truth_partition_acceptance.json`: the accepted fit uses exactly 1377 groups. |
| the old Tahoe native column numbering (`native_column` = token − 3) | not a defect and not an input | `bases/raw.py:tahoe_axis_mapping()` documents `native_column` as the immutable `gene_metadata` row; the token id is a separate column. Recorded so the offset is not "repaired". |

The current eligible Tahoe raw counts (`tahoe-raw-counts/`) and the current
Tahoe training corpus (`training-corpora/tahoe_bulk.h5ad`) remain accepted for
**training only**; they still carry no significance labels.

## E. Bounded expression checks deferred to #156

#155 performed none of these. Each is bounded and must run serially under
`systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G --`, with
`OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=2` and `/usr/bin/time -v`.

| id | check | why | measured cost |
|---|---|---|---|
| E1 | Re-hash the large payloads once after the move (≤20.4 GB read). | proves the move did not alter bytes | ~1–2 min, <1 GiB |
| E2 | Re-issue + verify the truth successor envelope (create then verify, serial). | its recorded paths are runtime dependencies | 4:20 + 4:21, 867/870 MiB |
| E3 | `verify_hvg_metrics.py` against the stable paths (rebuilds rankings and re-checks the held-out isolation intervention). | confirms the HVG product still binds the moved training corpora | 10–11 s, 245 MiB |
| E4 | Render the diagnostic panel from the CSV alone and compare PNG bytes. | ADR 0062 CSV-only figure reproduction, currently impossible (F6) | seconds |
| E5 | `evaluate_predictions.py --successor … --diagnostic` end-to-end on the accepted successor. | proves the public evaluator works through stable paths | bounded, ≤1 GiB |
| E6 | Only if E1–E5 contradict an accepted record: a bounded per-contrast raw-count replay of the affected truth groups, then refit if a numerical defect is confirmed. | ADR 0062 repair rule | reserved slot |

Not required: refitting truth, retraining either corpus, re-hashing the 3388
Tahoe source shards, re-running `verify_all.sh`, or rescanning expression for a
second opinion.

## F. Evidence hygiene

- A successful build or a matching hash alone never establishes acceptance; the
  scientific criterion plus downstream read plus figure/CSV reproduction are the
  gate (ADR 0062).
- Reused evidence stays labelled as reused: `native_tenv_verify`,
  `shared_v2_verify`, `training_v2_{build,verify,tenv_replay,tahoe_replay}`,
  `truth_complete_verify`, `truth_partition_acceptance`,
  `truth_successor_{acceptance,verify}`, `hvg_actual_{build,verify}`,
  `hvg_metrics_actual_acceptance`, `tahoe_counts_source_acceptance`,
  `diagnostic_publication`. #155 re-checked their digest bindings only.
- Stale point-in-time fields (for example `v2_actual_acceptance.json`'s
  `pending` list) are read as history, not as current state.
- Tolerances stay separate from byte integrity: metric values are deterministic
  here (VQA records image hashes), so no stochastic tolerance is claimed.
