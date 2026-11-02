# 06.7_cm_tahoe — 45 Tahoe-init CM candidates and the per-fold LR selection (ticket #7)

**Status: executed and verified** (2026-11-02 19:15–19:31 Asia/Shanghai, after the
ticket-#5 formal sweep).  The full Tahoe-initialisation candidate matrix — **5 folds × 3
seeds × 3 learning rates = 45 full-parameter CM finetunes**, each warm-starting from its
seed's formal pretrain checkpoint — ran on the frozen `06.3_input` view under the
ticket-#4 trainer: **45/45 completed**, every candidate is recorded with its warm-start
provenance, and the per-fold common learning rate is locked.  The substep mirrors
[`06.6_cm_random`](../06.6_cm_random/README.md) exactly — same folds, seeds, learning-rate
candidates, recipe, budget, early stopping and per-fold LR rule — changing **only the
initialisation**.  This substep trains **only** Tahoe-init candidates; the test split is
untouched until ticket #9.

The decision is [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md); the study
design is [`stage2/06_prnet_transfer/README.md`](../README.md); the frozen input, folds,
consensus controls and locked config are [`06.3_input`](../06.3_input/README.md); the
trainer and every metric are [`06.4_pipeline`](../06.4_pipeline/README.md); the pretrain
checkpoints are [`06.5_tahoe_pretrain`](../06.5_tahoe_pretrain/README.md).  This substep
adds **no training or evaluation code** — it defines the matrix, validates the warm start,
drives the shared queue, and records/gates the outcome.

## What is delivered

| deliverable | where |
|---|---|
| the strict warm-start loader + provenance record | [`execution/warmstart.py`](execution/warmstart.py) (wraps `06.4`'s `ptlib.checkpoint`; the trainer itself records the same provenance in every run manifest) |
| the 45-unit matrix definition (folds × seeds × LRs from the frozen config, seed→checkpoint) | [`execution/matrix.py`](execution/matrix.py), `execution/build_queue.py` → `output/units.json` (written after the sweep) |
| every candidate's outcome, config, data identity, checkpoint, warm-start identity, validation curve | [`output/record/candidate_record.csv`](output/record/candidate_record.csv), `validation_curves.csv` |
| per-fold LR selection (three-seed mean, lower LR on tie) | `output/record/lr_scores.csv`, `lr_selection.csv` |
| the **15 selected checkpoints** ticket #9 scores | `output/record/selected_checkpoints.csv` |
| run-level evidence (scope, pre-flight, wall, peaks, resume proof, gates, warm-start block) | `output/record/matrix_run.json`, `resume_check.json`, `record.json`, `record_check.json` |
| bounded wrapper profile on the first available formal checkpoint | [`output/profile/profile_record.json`](output/profile/profile_record.json) (tracked) |

## Matrix and recipe (locked, no re-tuning)

| factor | values |
|---|---|
| folds | 0–4 (the same frozen drug-key folds; same-plate matched DMSO controls for training, the fold's fixed consensus `c*` for validation) |
| seeds | 20261021, 20261022, 20261023 — each seed's **own** `#5` pretrain checkpoint; all reported, none dropped |
| learning rates | 1e-4, 3e-4, 1e-3 (ascending candidate order) |
| candidates | 5 × 3 × 3 = **45 full-parameter Tahoe-init finetunes** |
| warm start | `tahoe_seed<seed>/checkpoint/model.pt` from `06.5_tahoe_pretrain/output/formal/runs/`, strict axis/structure-checked load (below) |
| structure | PRnet `PGM`: hidden [128], z 64, c 64, noise 10, drug 1024, adaptor [128], dropout 0.05 → **7,185,000 params** |
| optimisation | Gaussian NLL (unclamped), Adam, weight_decay 1e-8, batch 512, max_epochs 100, min_epochs 20, early-stop patience 20, `ReduceLROnPlateau(mode="max", factor 0.5, patience 10, min_lr 1e-8)` |
| selection signal | **validation direction AUPRC@200 (max)**, one fixed `c*` input; zero validation/inference noise |
| LR selection | per fold, the LR maximising the **mean over the three seeds** of each seed's best validation AUPRC@200; **lower LR on a tie** |

The matrix is derived from the same frozen `locked_config.json` and fails closed if the
fold count, the three seeds or the three LRs ever drift (`matrix.py`; unit-tested).  The
sibling `06.6` matrix is the same definition with `init=random`, so the branches are
paired by construction: **the same folds, seeds, LR candidates, recipe and budget**.

## The warm-start contract (ticket #7's core acceptance item)

`execution/warmstart.py` does not re-implement training or the weight load.  It wraps the
pipeline's own primitives and fails closed on:

1. **missing file / foreign payload** — `ptlib.checkpoint.read_checkpoint` (the
   `PAYLOAD_KEYS` contract: `model`, `state_dict`, `config`, `seed`, `n_genes`, …);
2. **gene-axis mismatch** — the checkpoint's `n_genes` and the locked structure's
   `n_genes` must both equal the run's gene axis (18,164);
3. **structure mismatch** — checkpoint structure vs the run structure through
   `ptlib.config.structure_compatible` (widths, latents, noise, drug dim, dropout);
4. **wrong provenance source** — the payload `config.stage` must be `pretrain`, the payload
   `seed` must equal the matrix seed, and the path must be
   `tahoe_seed<seed>/checkpoint/model.pt` (a fold can never warm-start another seed's weights);
5. **partial or silently random load** — `ptlib.checkpoint.load_provenance` compares the
   checkpoint state dict with the freshly built model **before** loading: every model
   parameter **and** buffer key must be present, shapely and dtype-identical, and finite;
   then it strict-loads and re-compares every tensor **bit-for-bit** after the load.
   Missing keys, unexpected keys, shape/dtype mismatches or non-finite tensors raise before
   a single optimisation step.  The manifest records the empty mismatch lists; a partially
   ignored load cannot be reported as a successful transfer;
6. **untrained checkpoint** — an anti-random-init canary compares the loaded weights with a
   fresh `random_init(model, seed)`: bit-identical weights mean the "pretrain" checkpoint
   was never trained, and the run is rejected instead of silently becoming a random-init
   candidate.

The provenance record (per seed, in `output/units.json` → `warmstart`, copied into
`record.json`; per candidate in the run's own `manifest.init` and in
`candidate_record.csv`) carries: file `sha256`, `size_bytes`, checkpoint payload seed and
stage, **config digest** (`sha256` of the exact checkpoint config block), checkpoint epoch /
best validation AUPRC, gene axis, parameter count, `state_dict_sha256` (format-independent
digest of keys+dtypes+shapes+bytes), the full **loaded layer/key list** (14 tensors: 11
parameter tensors + 3 BatchNorm buffers in the locked structure), the empty
missing/unexpected/shape/dtype/non-finite lists, `strict_load: true`,
`weights_identical: true`, and the anti-random-init comparison.  `write_record.py --check`
re-verifies all of it per candidate and re-hashes each seed checkpoint, so a changed seed
file or a swapped provenance block is a failed candidate, not a warning.

The `/implement-spec` first turn's evidence that this works on a **real** formal checkpoint
is [`output/profile/profile_record.json`](output/profile/profile_record.json): seed
20261021's file (sha256 `a2c83736…`, config digest `5ec00027…`, epoch 8) was loaded in the
training process with 14/14 keys, `missing=[] unexpected=[]`, `weights_identical=true`, and
the loaded weights are **not** a fresh random init.

## Execution

The executed sequence (reproduce it after the `#5` sweep's three `tahoe_seed*/checkpoint/model.pt`
are final):

```
cd stage2/06_prnet_transfer/06.7_cm_tahoe/execution
make init      # verify all three seed checkpoints + write output/init/tahoe_init.json
make queue     # output/units.json (45 units, each warm-starting its seed)
make matrix    # THE heavy run: systemd-run … -p MemoryMax=12G, 2 workers, 1 Hz peaks
make record    # derive output/record/* from output/runs/*
make check     # re-derive + ~27 gates
make resume    # resume proof (45/45 skip + interrupted-unit re-run), scoped
# full sequence: execution/run.sh all   (scoped -> record -> check)
```

The matrix command is
`systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G -- python execution/run_matrix.py
… --workers 2 --omp-threads 2`.  The wrapper refuses to start when `MemAvailable < 12 GiB`,
when more than one other heavy GPU process is running, or when the queue file is not the
complete 45-unit queue with a verified warm-start block.  `DATA_DIR`, `PRETRAIN_DIR`, `OUT`,
`WORKERS`, `SCOPE` are overridable; the frozen input and the `#5` checkpoints are gitignored
and live once (the Makefile falls back to the integration checkout).

### Bounded wrapper profile (executed 2026-11-02 16:00, while `#5` seed 20261021 was still training)

The first available formal checkpoint (seed 20261021, best epoch 8 at that moment) was
**snapshotted byte-identically** (`_stage_checkpoint`: source/copy sha256 equal, so the live
sweep's next `os.replace` cannot tear the read), verified by `check_init.py` as its own
process, and used for one unit (`cm_f0_s20261021_lr1e-3_tahoe`) with a 2-epoch bounded
budget through the same queue path as the matrix (`--workers 2`).  Profile runs live
under `output/profile/runs/`, so the formal `output/runs/` tree holds exactly the 45
matrix units:

| quantity | measured |
|---|---|
| command scope | `systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G` |
| warm-start proof | sha256 `a2c83736…`, config digest `5ec00027…`, epoch 8, gene axis 18,164 = 18,164, 14/14 keys (11 params + 3 buffers), missing `[]`, unexpected `[]`, shape `[]`, dtype `[]`, `weights_identical: true`, `matches_fresh_random_init: false` |
| probe of the load | the training process logged `warm start: loaded … (sha256 a2c83736…, config 5ec00027…, epoch 8, run seed 20261021, 14 keys: 11 params + 3 buffers, missing=[] unexpected=[], weights_identical=True)` |
| epochs | 2 (epoch 1: train 6.07 s cold, valid 0.41 s, valid AUPRC@200 0.0804; epoch 2: train 0.18 s, valid 0.40 s, 0.0826) |
| run wall | pipeline 9.47 s, wrapper 11.82 s |
| host peak RSS (process tree, 1 Hz) | **2,292.7 MiB** |
| scope cgroup peak memory | 1,582.8 MiB (cap 12 GiB — never approached) |
| run manifest peak RSS | 1,325.3 MiB (matches ticket #6's per-candidate 1.32–1.34 GiB) |
| GPU device peak | 2,913 MiB of 16,376 — **not attributable** to this unit: it overlaps the `#5` pretrain, and this WSL2 driver reports `[N/A]` for per-process GPU memory |
| pre-flight | `MemAvailable` 37.9 GiB, CUDA present, one other heavy GPU process (the `#5` sweep), allowed by the two-heavy-jobs rule |
| path comparison vs ticket #6 | same unit `(fold 0, seed 20261021, lr 1e-3)`: the resolved formal configs differ in **exactly one key — `init` (`random` → `pretrained`)**; all shared training flags are identical; random config sha `7a445671…` vs Tahoe expected config sha `c05d5ea1…` |
| result | `passed: true`, 0 problems |

This is **pipeline/throughput evidence only**: the 2-epoch validation score is not a model
claim and the run is not part of the matrix.

## Results

### 45-candidate outcome summary

**42 `ok`** (early-stopped on validation AUPRC), **3 `capped`** (consumed the whole
100-epoch budget without an AUPRC improvement trigger — the budget bound, not a failure),
**0 failed**.  All 45 have a manifest, a checkpoint whose bytes match its recorded sha256,
a complete validation curve, and a warm-start provenance block.  The Tahoe-init runs
generally peak early (the 1e-4 candidates at epochs 9–29) and the higher LRs run longer:

| LR | candidates | epochs min / mean / max | capped |
|---|---|---|---|
| 1e-4 | 15 | 35 / 39.9 / 50 | 0 |
| 3e-4 | 15 | 31 / 43.9 / 100 | 1 |
| 1e-3 | 15 | 28 / 55.3 / 100 | 2 |

The three `capped` candidates are all seed `20261021` (`fold2 lr3e-4`, `fold2 lr1e-3`,
`fold4 lr1e-3`); they are usable candidates with their best-validation checkpoint and are
labelled so they are never silently read as early-stopped runs.

### Per-fold LR selection (three-seed mean validation AUPRC@200)

| fold | 1e-4 | 3e-4 | 1e-3 | selected | margin over next |
|---|---|---|---|---|---|
| 0 | **0.090974** | 0.090711 | 0.089994 | 1e-4 | +0.000263 |
| 1 | **0.095643** | 0.094575 | 0.092251 | 1e-4 | +0.001068 |
| 2 | 0.065851 | 0.065605 | **0.069518** | 1e-3 | +0.003667 |
| 3 | 0.074176 | 0.074310 | **0.074471** | 1e-3 | +0.000162 |
| 4 | **0.084970** | 0.084633 | 0.083655 | 1e-4 | +0.000337 |

No tie-break was needed in any fold.  Unlike the random-init branch (which selected 1e-3
in every fold), Tahoe init selects the **lowest** LR in folds 0, 1 and 4 and the highest in
folds 2 and 3 — the optimum shifts with the already-trained weights, which is exactly what
the per-fold, three-seed rule is for.

**Paired reading against ticket #6 (not a selection input — the test split is still
unscored).** With the same folds/seeds/LRs/recipe/budget, the random-init branch's selected
three-seed means are higher in every fold:

| fold | random init (#6) | Tahoe init (#7) | difference |
|---|---|---|---|
| 0 | 0.099018 | 0.090974 | −0.008044 |
| 1 | 0.103012 | 0.095643 | −0.007369 |
| 2 | 0.075352 | 0.069518 | −0.005834 |
| 3 | 0.082602 | 0.074471 | −0.008131 |
| 4 | 0.089705 | 0.084970 | −0.004735 |

This is a **validation-selection** comparison only; whether the gap holds on the held-out
drugs is ticket #9's question, and a negative transfer result is a valid first-round
outcome under the study's negative-result clause.

### The 15 selected checkpoints (input to ticket #9)

One checkpoint per `(fold, seed)` at that fold's selected LR; all three seeds are kept.
Validation AUPRC is each seed's best epoch; full sha256 and checkout-relative paths are in
[`output/record/selected_checkpoints.csv`](output/record/selected_checkpoints.csv) (the
14-column schema ticket #9 consumes, identical to #6's).

| fold | seed | valid AUPRC@200 | best epoch | epochs run | stop |
|---|---|---|---|---|---|
| 0 | 20261021 | 0.093502 | 29 | 50 | early stop |
| 0 | 20261022 | 0.088336 | 15 | 36 | early stop |
| 0 | 20261023 | 0.091084 | 22 | 43 | early stop |
| 1 | 20261021 | 0.098518 | 17 | 38 | early stop |
| 1 | 20261022 | 0.093596 | 14 | 35 | early stop |
| 1 | 20261023 | 0.094815 | 17 | 38 | early stop |
| 2 | 20261021 | 0.069965 | 83 | 100 | **max_epochs (capped)** |
| 2 | 20261022 | 0.063297 | 9 | 30 | early stop |
| 2 | 20261023 | 0.075292 | 34 | 55 | early stop |
| 3 | 20261021 | 0.076267 | 29 | 50 | early stop |
| 3 | 20261022 | 0.072772 | 17 | 38 | early stop |
| 3 | 20261023 | 0.074375 | 15 | 36 | early stop |
| 4 | 20261021 | 0.087221 | 25 | 46 | early stop |
| 4 | 20261022 | 0.082700 | 16 | 37 | early stop |
| 4 | 20261023 | 0.084988 | 20 | 41 | early stop |

15 checkpoints × 28,746,059 bytes = **431,190,885 bytes**; selected_checkpoints.csv sha256
`8a9a6d6c5e7e1063b5177ea119780b6ebc1a13d8dbaf89fc6bb950ad77fe3c17`.

### Warm-start provenance (executed)

All 45 candidates loaded the matching seed's pretrain checkpoint exactly; the record's
`warmstart` block (also in `units.json`) carries the proof, and every candidate manifest
carries the same proof in `init`:

| seed | checkpoint sha256 | config digest | pretrain epoch | pretrain val AUPRC | keys loaded | missing / unexpected | weights identical | fresh-random canary |
|---|---|---|---|---|---|---|---|---|
| 20261021 | `a2c83736cf4faed9…` | `5ec0002761b5b21a…` | 8 | 0.200751 | 14 (11 params + 3 buffers) | `[]` / `[]` | true | false |
| 20261022 | `e911e0545586a2bf…` | `9c3c9514c1af8c3d…` | 74 | 0.222222 | 14 (11 params + 3 buffers) | `[]` / `[]` | true | false |
| 20261023 | `5cf773515b8e8f0f…` | `55b19d1b7b81a210…` | 32 | 0.231263 | 14 (11 params + 3 buffers) | `[]` / `[]` | true | false |

Three-seed provenance sha256 `3b52ea219b97db1d5a9f891a25e9adc724d181d6778ee8e3c4b39c648e26a58d`;
gene axis 18,164 = 18,164 for all; every checkpoint is a `stage: pretrain` payload whose
payload seed equals the matrix seed.  The three seed files are distinct and unchanged since
verification (`init_checkpoint_bytes_match`), and each fold×seed used exactly its own seed's
checkpoint (`seed_checkpoint_one_to_one`).

### Data identity and branch symmetry

Data identity `a57b0e35b1126d0b88da361a43ff27479f42e005a5291cd7641d3fa52cb7b376` — the
same frozen bytes as ticket #6.  The `matrix_symmetry_with_ticket6` gate confirms the
folds, seeds, LRs, recipe, structure, metric, loss, optimiser, noise policy, budget and
selection rule are identical to #6's record; only `init` (`random` → `pretrained`) differs.

### Execution (one heavy run, alone, 12 GiB scope)

| quantity | measured |
|---|---|
| command | `systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G -- … run_matrix.py … --workers 2 --omp-threads 2` |
| wall time | **923.6 s (15 min 24 s)** for all 45 units, 2 workers |
| queue summary | `{units:45, completed:45, skipped_complete:0, failed:0}` |
| host peak RSS (process tree, 1 Hz) | **3,514 MiB** |
| scope cgroup peak memory | 3,748 MiB (cap 12 GiB — never approached) |
| GPU device peak | **3,412 MiB** of 16,376 (idle baseline ≈1,300 MiB) |
| pre-flight | 44.3 GiB `MemAvailable`, no other heavy GPU process, CUDA present |
| resume proof (`make resume`) | 45/45 skipped when complete; an interrupted unit re-ran with the identical config fingerprint (`a69e12cb…`) and was then skipped |

**First attempt (kept for the audit trail, not the recorded run).** The first matrix
launch completed 44/45 units and lost `cm_f0_s20261021_lr1e-4_tahoe` to a genuine
concurrency defect in the shared frozen-artifact hash memo: two workers wrote the same
`hash_cache.json.tmp`, the first `replace` moved it away and the second raised
`FileNotFoundError` inside `open_frozen_input` (907 s, peaks 3,527 MiB host / 3,265 MiB
GPU).  Fixed in `06.4_pipeline` (`d691fbf`: unique `tempfile.mkstemp` temp per writer,
`os.replace`, regression-tested with 160 concurrent saves); the recorded run is the clean
45/45 re-run on the fixed code.  First-pass matrix/status/log kept locally in
`output/first_attempt/`.

## Evidence and gates

`make check` re-derives the record and applies the 06.6 gate set plus the ticket-#7
additions; **all 27 gates pass** (`output/record/record_check.json`, `passed: true`, 0
failed):

| gate group | what it proves |
|---|---|
| `matrix_complete`, `queue_status_present` | 45/45 candidates present in matrix order; queue status rows for all 45 |
| `data_identity_uniform` | every complete candidate used the same frozen bytes |
| `config_one_to_one`, `checkpoint_one_to_one`, `checkpoint_bytes_match` | 45 distinct configs and 45 distinct candidate checkpoints, bytes still match |
| `warmstart_verified`, `warmstart_provenance_complete`, `per_run_init_provenance` | three seed checkpoints verified; each candidate's manifest proves the exact matching-seed load (axis, config digest, sha256, keys, empty mismatch lists, bit-identical weights) |
| `init_checkpoint_bytes_match`, `init_checkpoints_distinct`, `seed_checkpoint_one_to_one` | the three seed files are unchanged since verification, distinct, and each fold×seed uses exactly its own seed's checkpoint |
| `matrix_symmetry_with_ticket6` | folds/seeds/LRs/recipe/budget/selection rule identical to ticket #6's record — only the initialisation differs |
| `selection_recomputes`, `selection_complete`, `selected_checkpoint_count`, `selected_three_seeds_per_fold` | the recorded selection follows from the three-seed means; 15 selected = 3 seeds per fold |
| `validation_curves_consistent`, `epoch_budget_respected` | one curve row per epoch, curve max equals the manifest best at the best epoch, 20 ≤ epochs ≤ 100 |
| `no_prediction_artifacts_exist`, `no_failed_candidates` | no run directory holds a prediction/evaluation tree; no candidate failed (a failed candidate fails the matrix) |
| `record_content_sha256_matches`, `record_artifacts_match`, `record_files_present` | the tracked record is byte-identical to a fresh derivation |
| `matrix_run_recorded`, `resume_check_passed` | a completed 45-unit run is recorded; the resume proof passed |

**The test split is not read for training or selection.**  As in `06.6`: the only process
the queue launches is `train.py --stage finetune` (its `CMDataset` materialises only the
fold's train/valid rows and forbids test rows; ticket #4's leakage suite covers it), no
candidate run directory holds a prediction/evaluation artifact, and selection uses only
`fit_history.json`'s validation AUPRC curve.  The test role is scored once, after
selection, by ticket #9.

## Reproduction

From `stage2/06_prnet_transfer/06.7_cm_tahoe/execution/`:

```bash
make test          # 67 offline stdlib tests (arc env, CPU, synthetic checkpoints/runs)
make tables        # derived CSV cache of the frozen parquet tables (pinned env has no engine)
make init          # verify the three #5 checkpoints (fail closed)
make queue         # output/units.json (45 units, every seed's own checkpoint)
make plan          # resume-aware dry run of the queue
make matrix        # THE heavy run: systemd-run 12 GiB scope, 2 workers, 1 Hz peaks
make profile       # bounded one-unit wrapper profile (2 epochs, staged snapshot)
make record        # derive output/record/* from output/runs/*
make check         # re-derive + gate the record
make resume        # resume proof, scoped
# full sequence: execution/run.sh all   (scoped -> record -> check)
```

## Environment, memory and cost

- GPU env `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python` (py3.9, torch
  2.0.1+cu118) for the warm-start verification and every training unit; tests and the
  record writer use the `arc` env.
- The matrix is expected to be the same footprint as `06.6`: ≈3.0 GiB host / ≈4.0 GiB GPU
  for two workers, in one `systemd-run MemoryMax=12G` scope.  The bounded profile measured
  a single unit at 1.33 GiB (manifest) / 2.29 GiB process tree.  The matrix must run alone
  (after the sweep): the repo-wide rule is two heavy jobs at most and combined ≤ ~28 GiB.
- Warm-start verification loads a 28.7 MB CPU checkpoint and builds one 7.185 M-parameter
  model per seed; it does not touch the expression matrices.  Run manifests are read through
  the ticket-#4 loaders; nothing here materialises an observation × 18,164 dense array.
- The run directories (45 × ~29 MB checkpoints + logs/curves) are local/gitignored; the
  tracked record is small (≤ 1 MB).

## Boundaries and notes

- **No science claim here.**  These are candidate *selection* numbers on the validation
  split.  Whether Tahoe init improves prediction of unseen-drug CM responses, and how the
  two branches compare with the baselines, is ticket #9.
- The same negative-result clause as the study: a Tahoe branch that does not beat random
  init is a valid first-round outcome, reported as such.
- All three seeds and all three LRs are reported; no seed or LR is dropped, and the LR
  selection never uses test scores.
- The warm-start file is a **stage `pretrain`** checkpoint only; a CM checkpoint or a
  checkpoint from another seed is rejected by the loaders, not repaired.
- The GitHub issue #7 stays open until the Lead's integration and the 45-unit matrix; the
  `#5` sweep and `#9` evaluation are unaffected by this substep except through the 15
  selected checkpoints.

## Tracked vs local

Tracked: `AGENTS.md`, `README.md`, `execution/**` (matrix/warm-start/queue/record scripts,
Makefile, `run.sh`, 67 tests), `output/README.md`, and the force-added evidence:
`output/profile/profile_record.json` (bounded wrapper profile), `output/units.json` (the
45-unit matrix + verified warm-start block), `output/init/tahoe_init.json` (three-seed
provenance), and all nine files of `output/record/`.  Local/gitignored: `*.pt` weights
under `output/runs/`, `output/profile/staged_checkpoints/`, `output/profile/runs/`,
`output/input_cache/`, `output/hash_cache.json`, `output/queue_status.csv`, `output/logs/`,
`output/_resume_probe/`, `output/profile/logs/`, `output/first_attempt/`.

## Review deviations

- **Deliberate per-substep duplication.** The `common.py` helpers and the matrix/record
  definitions are repeated here instead of shared: ADR 0030 §3 gives every substep an
  independently re-derivable execution tree, and the single shared trainer/evaluator lives
  in [`06.4_pipeline`](../06.4_pipeline/README.md). Accepted as a maintainability cost.
- **Results under `output/`, not the substep root.** Record CSVs (and, in 06.9, figures)
  live in `output/record/` / `output/figures/` with a tracked record force-added, rather
  than directly in the substep directory as ADR 0030 §3 asks. Accepted deviation: moving
  them after the records were locked would break the recorded artifact paths.
- **Checkpoint identity travels as columns.** `config_sha256`, `checkpoint_sha256`,
  `checkpoint_bytes` and `checkpoint_worktree_relpath` are explicit CSV/JSON fields rather
  than a dedicated identity type. Accepted for CSV/record compatibility.
