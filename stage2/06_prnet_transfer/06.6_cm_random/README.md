# 06.6_cm_random — 45 random-init CM candidates and the per-fold LR selection (ticket #6)

**Status: executed and verified** (2026-11-02 15:24–15:42 Asia/Shanghai). The full
random-initialisation candidate matrix — **5 folds × 3 seeds × 3 learning rates = 45
full-parameter CM finetunes** — ran on the frozen `06.3_input` view under the ticket-#4
trainer, 45/45 completed, every candidate is recorded, and the per-fold common learning
rate is locked. This substep trains **only** random-init candidates; the Tahoe-init
mirror is ticket #7, and the test split is untouched until ticket #9.

The decision is [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md); the study
design is [`stage2/06_prnet_transfer/README.md`](../README.md); the frozen input, folds,
consensus controls and locked config are [`06.3_input`](../06.3_input/README.md); the
trainer and every metric are [`06.4_pipeline`](../06.4_pipeline/README.md). This substep
adds **no training or evaluation code** — it defines the matrix, drives the shared
queue, and records/gates the outcome.

## What is delivered

| deliverable | where |
|---|---|
| the 45-unit matrix definition (folds × seeds × LRs from the locked config) | [`output/units.json`](output/units.json) (sha256 `503b0041…`) + [`execution/matrix.py`](execution/matrix.py) |
| every candidate's outcome, config, data identity, checkpoint path+sha256, validation curve | [`output/record/candidate_record.csv`](output/record/candidate_record.csv), [`validation_curves.csv`](output/record/validation_curves.csv) |
| per-fold LR selection (three-seed mean, lower LR on tie) | [`output/record/lr_scores.csv`](output/record/lr_scores.csv), [`lr_selection.csv`](output/record/lr_selection.csv) |
| the **15 selected checkpoints** ticket #9 scores | [`output/record/selected_checkpoints.csv`](output/record/selected_checkpoints.csv) |
| run-level evidence (scope, pre-flight, wall, host/GPU peaks, resume proof, gates) | [`output/record/matrix_run.json`](output/record/matrix_run.json), [`resume_check.json`](output/record/resume_check.json), [`record.json`](output/record/record.json), [`record_check.json`](output/record/record_check.json) |

## Matrix and recipe (locked, no re-tuning)

| factor | values |
|---|---|
| folds | 0–4 (frozen drug-key folds; same-plate matched DMSO controls for training, the fold's fixed consensus `c*` for validation) |
| seeds | 20261021, 20261022, 20261023 (all reported; none dropped or cherry-picked) |
| learning rates | 1e-4, 3e-4, 1e-3 (ascending candidate order) |
| candidates | 5 × 3 × 3 = **45** random-init full-parameter finetunes |
| structure | PRnet `PGM`: hidden [128], z 64, c 64, noise 10, drug 1024, adaptor [128], dropout 0.05 → **7,185,000 params** |
| optimisation | Gaussian NLL (unclamped), Adam, weight_decay 1e-8, batch 512, max_epochs 100, min_epochs 20, early-stop patience 20, `ReduceLROnPlateau(mode="max", factor 0.5, patience 10, min_lr 1e-8)` |
| selection signal | **validation direction AUPRC@200 (max)**, one fixed `c*` input; zero validation/inference noise |
| LR selection | per fold, the LR maximising the **mean over the three seeds** of each seed's best validation AUPRC@200; **lower LR on a tie** |

The matrix is derived from the frozen `locked_config.json` and fails closed if the folded
counts, the three seeds or the three LRs ever drift (`matrix.py`; unit-tested).

## Execution

The 45 units ran as **one scoped job**, through the ticket-#4 `run_queue.py` (resume-safe,
`--workers 2 --omp-threads 2`), each unit a separate `train.py` subprocess:

```
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G -- \
  python execution/run_matrix.py ... --workers 2 --omp-threads 2
```

| quantity | measured |
|---|---|
| wall time | **981.8 s (16 min 22 s)** for all 45 units, 2 workers |
| per-candidate wall | 23.0–81.7 s, mean **41.8 s** (sum 1,880 s ⇒ ≈1.9× parallel efficiency) |
| host peak RSS (process tree, 1 Hz) | **3,088 MiB** |
| scope cgroup peak memory | 3,459 MiB (cap 12 GiB — never approached) |
| GPU peak (whole GPU) | **4,003 MiB** of 16,376 (idle baseline 1,063 MiB ⇒ ≈2.94 GiB training) |
| per-candidate host peak | 1,322–1,340 MiB (matches ticket #4's 1.33 GiB) |
| pre-flight | 44.6 GiB `MemAvailable`, 0 other heavy GPU processes, CUDA present |
| failures | 0; `queue_summary = {units:45, completed:45, failed:0}` |

Only one heavy job ran in this scope; it is compatible with ticket #5's pretrain running
in parallel (≈3.1 GiB + ≈2.7 GiB ≪ 47 GiB), and the wrapper refuses to start if
`MemAvailable < 12 GiB` or if more than one other heavy GPU process is present.

## Results

### 45-candidate outcome summary

**40 `ok`** (early-stopped on validation AUPRC) and **5 `capped`** (ran the whole 100-epoch
budget without an AUPRC improvement trigger — the budget bound, not a failure);
**0 failed**. All 45 have a manifest, a checkpoint whose bytes match its recorded sha256,
and a complete validation curve. Epochs grow with the learning rate, as expected:

| LR | candidates | epochs min / mean / max | capped |
|---|---|---|---|
| 1e-4 | 15 | 22 / 29.7 / 100 | 1 |
| 3e-4 | 15 | 22 / 44.7 / 100 | 1 |
| 1e-3 | 15 | 27 / 67.7 / 100 | 3 |

Per fold × LR, the **three-seed mean validation AUPRC@200** (each mean is over all three
seeds; `n=3` for every cell):

| fold | 1e-4 | 3e-4 | 1e-3 | selected | margin over next |
|---|---|---|---|---|---|
| 0 | 0.095906 | 0.096418 | **0.099018** | 1e-3 | +0.002600 |
| 1 | 0.101828 | 0.101925 | **0.103012** | 1e-3 | +0.001087 |
| 2 | 0.069567 | 0.069959 | **0.075352** | 1e-3 | +0.005393 |
| 3 | 0.076336 | 0.077056 | **0.082602** | 1e-3 | +0.005545 |
| 4 | 0.089433 | 0.089516 | **0.089705** | 1e-3 | +0.000189 |

**Every fold selects LR = 1e-3**, with no tie and no fallback: the three-seed mean rises
monotonically with the learning rate in all five folds. Fold 4's margin over 3e-4 is
narrow (1.9e-4) and is reported as such; the rule was not applied ex post. Fold 1 is the
strongest valid fold (0.1030) and fold 2 the weakest (0.0754) — a fold spread to carry
into ticket #9's per-drug pairing, not a selection criterion.

### The 15 selected checkpoints (input to ticket #9)

One checkpoint per (fold, seed) at that fold's selected LR; all three seeds are kept.
Validation AUPRC is each seed's best epoch; full sha256 and worktree-relative paths are
in [`selected_checkpoints.csv`](output/record/selected_checkpoints.csv).

| fold | seed | valid AUPRC@200 | best epoch | epochs run | stop |
|---|---|---|---|---|---|
| 0 | 20261021 | 0.096834 | 10 | 31 | early stop |
| 0 | 20261022 | 0.097841 | 9 | 30 | early stop |
| 0 | 20261023 | 0.102379 | 62 | 83 | early stop |
| 1 | 20261021 | 0.102841 | 48 | 69 | early stop |
| 1 | 20261022 | 0.103244 | 49 | 70 | early stop |
| 1 | 20261023 | 0.102952 | 69 | 90 | early stop |
| 2 | 20261021 | 0.075247 | 60 | 81 | early stop |
| 2 | 20261022 | 0.075816 | 88 | 100 | **max_epochs (capped)** |
| 2 | 20261023 | 0.074992 | 66 | 87 | early stop |
| 3 | 20261021 | 0.082520 | 99 | 100 | **max_epochs (capped)** |
| 3 | 20261022 | 0.083244 | 97 | 100 | **max_epochs (capped)** |
| 3 | 20261023 | 0.082042 | 66 | 87 | early stop |
| 4 | 20261021 | 0.089499 | 6 | 27 | early stop |
| 4 | 20261022 | 0.089716 | 11 | 32 | early stop |
| 4 | 20261023 | 0.089899 | 8 | 29 | early stop |

15 checkpoints × 28,746,059 bytes = 431 MB. The training manifests keep
`best_valid_auprc`/`best_epoch`; each checkpoint file is the best-validation-epoch weights,
not the last epoch.

### Data identity

All 45 runs trained on the identical frozen bytes, and the record cross-checks it:
frozen `manifest.json` sha256 `3050d149…`, `locked_config.json` sha256 `8b48ba76…`,
data identity `a57b0e35…` (sha256 over manifest + locked config + the eight artifact
hashes). Each run's own manifest records the verified artifact set; the seeded hash
cache only memoised the 16.4 GB re-hash by `(size, mtime_ns)` — it was invalidated
and rebuilt on any mismatch, and the record's data identity is derived from the frozen
files, not from the cache.

## Evidence and gates

`make check` re-derives the whole record from the run directories and applies **20
gates**; all 20 pass (`output/record/record_check.json`, `passed: true`):

| gate group | what it proves |
|---|---|
| `matrix_complete`, `queue_status_present` | 45/45 candidates present in matrix order, no missing/extra unit, queue status rows for all 45 |
| `data_identity_uniform` | every complete candidate used the same frozen bytes; the identity equals the frozen input's |
| `config_one_to_one`, `checkpoint_one_to_one` | 45 distinct config fingerprints and 45 distinct checkpoint sha256 — no candidate reuses another's config or weights |
| `checkpoint_bytes_match` | every checkpoint file still hashes to its manifest value |
| `selection_recomputes`, `selection_complete`, `selected_checkpoint_count`, `selected_three_seeds_per_fold` | the recorded selection is exactly what the recorded three-seed means imply; 5 complete folds; 15 selected = 3 seeds per fold at the chosen LR |
| `validation_curves_consistent`, `epoch_budget_respected` | each curve has one row per epoch, its max equals the manifest best at the manifest best epoch, and 20 ≤ epochs ≤ 100 |
| `no_prediction_artifacts_exist`, `no_failed_candidates` | no run directory holds a `prediction/` or `evaluation/` tree, and no candidate failed |
| `record_content_sha256_matches`, `record_artifacts_match`, `record_files_present` | the tracked record on disk is byte-identical to a fresh derivation |
| `matrix_run_recorded`, `resume_check_passed` | the run-level record shows a completed 45-unit queue, and the resume proof passed |

**The test split was not read for training or selection.** Three independent arguments:
(1) the only process the queue launches is `train.py --stage finetune`, whose `CMDataset`
materialises only the fold's **train/valid** rows and treats test rows as forbidden —
ticket #4's `tests/test_data_leakage.py` (20 tests, re-run in this worktree: green)
includes `test_only_train_and_valid_rows_are_materialised` and
`test_training_loader_refuses_test_rows`; (2) no candidate run directory contains a
prediction/evaluation artifact (gate above); (3) selection uses only
`fit_history.json`'s validation AUPRC curve. The test role is scored once, after
selection, by ticket #9.

### Resume proof

`make resume` (recorded in `output/record/resume_check.json`, `passed: true`):

1. re-running the queue over the finished matrix reports **45 units, 45 skipped as
   complete, 0 planned** (the skip predicate is manifest + config fingerprint + checkpoint
   sha256, so a changed or tampered unit cannot be skipped);
2. a finished run copied with `manifest.json`/`fit_history.json`/`checkpoint/` removed —
   the state an interrupted unit leaves — is **re-run to completion** by the queue, and its
   re-derived `config_sha256` (`6241da7d…`) equals the original run's;
3. the recovered unit is then **skipped** on the next pass.

An interrupted candidate is therefore continued (re-run cleanly, not resumed mid-epoch);
completed candidates are never recomputed. The probe trees live in `output/_resume_probe/`
(local only) and are not part of the matrix.

## Reproduction

From `stage2/06_prnet_transfer/06.6_cm_random/execution/`:

```bash
make test          # 31 offline unit tests (arc env, synthetic runs; no GPU, no input)
make tables        # derived CSV cache of the frozen parquet tables (pinned env has no engine)
make queue         # output/units.json  (45 units from the frozen locked config)
make plan          # resume-aware dry run of the queue
make matrix        # THE heavy run: systemd-run 12 GiB scope, 2 workers, 1 Hz peaks
make record        # derive output/record/* from output/runs/*
make check         # re-derive + 20 gates
make resume        # resume proof (45/45 skip + interrupted-unit re-run), scoped
# full sequence: execution/run.sh all   (scoped -> record -> check)
```

`DATA_DIR` defaults to this worktree's `06.3_input/output` when the 16.4 GB CM view is
present, else to the canonical `/home/chang/arc/stage2/06_prnet_transfer/06.3_input/output`
(the frozen input is gitignored and exists once). `OUT`, `WORKERS`, `SCOPE` and both
interpreters are overridable; `GPU_PY` must be the pinned CUDA env.

## Environment, memory and cost

- GPU env `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python`
  (py3.9, torch 2.0.1+cu118); tests and the record writer use the `arc` env.
- The whole matrix is **3.0 GiB host / 4.0 GiB GPU** — one `systemd-run`
  `MemoryMax=12G` scope, measured, never near the cap. It does not duplicate ticket #5's
  full Tahoe-valid evaluation (8.5 GiB / 156 s), and it never materialises an
  observation × 18,164 dense array; per-fold row-block reads come through the ticket-#4
  loaders and the frozen tables are read from the derived CSV cache.
- The 1.3 GB of run directories (checkpoints, logs, curves) and the 4.9 GB Tahoe memmap are
  local/gitignored; only the 0.8 MB record is tracked.

## Boundaries and notes

- **No science claim here.** These are candidate *selection* numbers on the validation
  split. Whether random-init PRnet predicts unseen-drug CM responses, and how it compares
  with the Tahoe-init branch and the baselines, is ticket #9.
- `capped` (5 candidates) means "used the full epoch budget without early stopping"; they
  are usable candidates with their best-validation checkpoint, and they are labelled so
  they are never silently read as early-stopped runs.
- The early-stopped candidates all stop at exactly `best_epoch + 21` (patience 20 plus the
  stop epoch), confirmed for all 40 `ok` rows — the stopping rule is mechanical, not tuned
  per candidate. The lower LRs peak early (fold 0's 1e-4 seeds at epochs 9–10) and then
  spend the patience; the 1e-3 seeds often keep improving longest (up to epoch 99).
- Checkpoint files are gitignored local weights. The record carries their sha256 and
  worktree-relative paths, so ticket #9 must score a tree where they exist (copy
  `output/runs/` or keep this worktree) and re-verify the hashes.
- The run was launched from integration commit `e0e2d68` with this substep's driver scripts
  still uncommitted (they are committed immediately after); the training path was the
  committed, unmodified `06.4_pipeline`. `matrix_run.json` records `git.dirty: true` for
  that reason.
- The GitHub issue #6 stays open until the lead's integration; `06.5`/`06.7`/`06.9` are
  unaffected by this substep except through the 15 checkpoint files.

## Tracked vs local

Tracked: `AGENTS.md`, `README.md`, `execution/**` (matrix/queue/record scripts, Makefile,
`run.sh`, 31 tests), `output/units.json`, and the force-added
`output/record/` record of 2026-11-02 (`candidate_record.csv`, `validation_curves.csv`,
`lr_scores.csv`, `lr_selection.csv`, `selected_checkpoints.csv`, `record.json`,
`matrix_run.json`, `resume_check.json`, `record_check.json`).
Local/gitignored: `*.pt` weights under `output/runs/`, `output/input_cache/`,
`output/hash_cache.json`, `output/queue_status.csv`, `output/logs/`, `output/_resume_probe/`,
`output/record/*.tmp`.
