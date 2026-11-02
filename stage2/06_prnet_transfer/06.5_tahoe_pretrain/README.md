# 06.5_tahoe_pretrain — three seeded formal Tahoe pretrains (ticket #5)

This substep owns the **formal Tahoe pretraining** evidence for the PRnet
formal-transfer study (ADR 0032, decision 3/5): three independent pretrains of
the author PRnet `PGM` on the frozen 06.3 Tahoe view, one per frozen seed,
selected and early-stopped on **validation direction AUPRC@200 (max)**, with the
matched `(plate, cell_line)` DMSO controls and the frozen `plate13`+`plate14`
validation split. Each seed's checkpoint is reused across the five CM folds by
ticket #7.

It contains **no training code**. The one trainer/evaluator lives in
[`06.4_pipeline/execution/`](../06.4_pipeline/README.md); this substep is a thin
wrapper plus the evidence record around three `train.py --stage pretrain`
invocations (ticket #4's CLI). The wrapper imports `06.4/execution/ptlib` for the
config resolution and digest — it never copies or re-derives the recipe.

## Status

| | |
|---|---|
| substep | **built, wrapper verified on the real data** (2026-11-02). The three-seed formal sweep is **ready to launch** and is Lead-run; no formal checkpoint exists yet in this branch. |
| profile evidence | `output/profile/profile_record.json` — two short real-data units (the CLI's bounded 2-epoch smoke, and one **full-width** epoch). Peaks: host **11,273 MiB**, GPU device **4,003 MiB** (concurrent with ticket #6); the 12 GiB `systemd-run` cap was never tripped |
| wrapper | `execution/run_pretrain.py` (sweep) + `execution/plan.py` (06.4 config resolution & digests) + `execution/record.py`/`check_record.py` (record writer/checker) |
| tests | `make test` — 23 offline stdlib `unittest` cases on synthetic fixtures (record, budget, digests, curve, tamper cases); no GPU, no frozen data |
| what remains | the Lead launches `execution/run.sh formal` (three frozen seeds, sequential, locked 100-epoch budget) and force-adds the resulting record + curves + model cards |

A **profile** record is pipeline/throughput evidence only: it is produced with
2-epoch-small-subset and 1-epoch-full-width budgets, never with the formal locked
budget. The checker refuses to treat it as the formal record.

## Layout

```
06.5_tahoe_pretrain/
├── AGENTS.md                       # routing (ADR 0030, ≤ 30 lines)
├── README.md                       # this file
├── execution/
│   ├── run_pretrain.py             # the sweep wrapper (preflight → per-seed 06.4 run → record)
│   ├── plan.py                     # resolves each unit through 06.4 ptlib (digests, budgets)
│   ├── record.py                   # stdlib record writer/reader/checker (schema 06.5…/record/1)
│   ├── check_record.py             # strict record validator
│   ├── run.sh Makefile             # one-command entry points, systemd-run memory caps
│   └── tests/                      # offline stdlib unittest (`make test`)
└── output/
    ├── README.md                   # tracked-vs-local table (force-added)
    ├── profile/                    # tracked short-profile record + curves + model cards
    └── formal/                     # the three-seed record + curves + model cards (after the sweep)
```

## The formal contract (frozen; never overridden by this wrapper)

| | |
|---|---|
| stage | `train.py --stage pretrain` from `06.4_pipeline` (one invocation per seed) |
| input | frozen `06.3_input/output/tahoe_transfer.h5ad` (sha256 `00f8d906…`), 67,018 × 18,164 |
| model input | matched `(plate, cell_line)` QC-pass DMSO control mean, log1p (649 groups, 60,861 pairs) |
| target | `layers['log1p_cp10k']` treated expression |
| validation split | frozen `tahoe_split` valid = `plate13` + `plate14` (8,524 usable rows / 8,110 conditions) |
| loss / optimiser | author Gaussian NLL (`gaussian_nll`, may be negative, never clamped) / Adam |
| learning rate | `lrs.pretrain` = **1e-3** from the locked config (no LR search in pretraining) |
| batch | 512 |
| budget | max **100** epochs, min **20**, early stop after **20** stale epochs |
| scheduler | `ReduceLROnPlateau(mode="max", factor=0.5, patience=10, min_lr=1e-8)` on validation AUPRC@200 |
| selection | validation direction **AUPRC@200 (max)**; checkpoint at the best epoch; the test/valid split never informs LR, epoch, seed or capacity |
| seeds | `20261021, 20261022, 20261023` (all three reported, no seed cherry-picking) |
| derived cache | `06.4_pipeline/output/cache/tahoe_log1p_cp10k_v1.npy` (4.87 GB float32 memmap; its sidecar must bind the frozen `tahoe_transfer.h5ad` sha256) |

The wrapper re-asserts this table (`record.TICKET_EXPECTATIONS`) against the
frozen `locked_config.json` before training and refuses to run if the frozen
input, the cache or the recipe disagree with it. It passes **no** recipe, LR,
epoch or row-limit flag in formal mode; `train.py` resolves the locked config
itself. The record stores the resulting config digest per seed and the checker
re-verifies it.

## How to run

From the repository root (or use the equivalent Makefile targets):

```bash
SUB=stage2/06_prnet_transfer/06.5_tahoe_pretrain/execution

# what would run? (no data read, no GPU)
$SUB/run.sh plan

# formal sweep: the three frozen seeds, sequential, locked budget, 18G systemd cap
$SUB/run.sh formal

# one seed only (writes tahoe_pretrain_record_partial.json — not formal evidence)
$SUB/run.sh formal --seeds 20261021

# short real-data profile: bounded 2-epoch smoke + one full-width epoch, 12G cap
$SUB/run.sh profile

# validate a written record
$SUB/run.sh check-formal
$SUB/run.sh check-profile
```

`run.sh` wraps the wrapper in
`systemd-run --user --scope -p MemoryMax=18G -p MemorySwapMax=18G` (12G for the
profile) by default; `--cap=24G` raises it and `--no-cap` assumes an outer
scope. The exact per-seed trainer command is printed before each run and stored
in the record (for seed `20261021` in the integration worktree it is
`…/06.4_pipeline/execution/train.py --stage pretrain --data-dir
…/06.3_input/output --out …/06.5_tahoe_pretrain/output/formal/runs/tahoe_seed20261021
--seed 20261021 --device cuda --cache-mode require --tahoe-cache
…/06.4_pipeline/output/cache/tahoe_log1p_cp10k_v1.npy --table-cache
…/06.4_pipeline/output/input_cache --verify-cache deep`).

### Prerequisites the wrapper refuses to run without

1. **Frozen input** `06.3_input/output/` with its `manifest.json` (13 artifacts,
   15.3 GB). Every `train.py` run re-verifies the artifact sha256s before
   training; the wrapper records that verification.
2. **Derived log1p memmap cache** bound to the frozen view. The memmap in this
   worktree's sibling `06.4_pipeline/output/cache/` is not checked in; a cache
   built from the 04 pretrain view must be rebound once (it is byte-identical,
   proven by a full-array checksum + sample equality):
   ```bash
   /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python \
     stage2/06_prnet_transfer/06.4_pipeline/execution/tools/build_cache.py \
     --data-dir stage2/06_prnet_transfer/06.3_input/output \
     --cache    stage2/06_prnet_transfer/06.4_pipeline/output/cache/tahoe_log1p_cp10k_v1.npy \
     --rebind-frozen
   ```
   The wrapper prints this exact fix (with absolute paths) when the sidecar's
   `source_sha256` is not the frozen `tahoe_transfer.h5ad` hash. Rebuilding from
   the frozen input instead is correct but is a ~9.7 GiB peak job — run it alone.
3. **Derived frozen-table CSV cache** `06.4_pipeline/output/input_cache/` (the
   pinned GPU env has no parquet engine): `make frozen-tables` in
   `06.4_pipeline/execution` with the `arc` env.
4. **RAM floor**: the wrapper aborts when `/proc/meminfo` `MemAvailable` <
   12 GiB (and warns below 16 GiB for a formal run), and aborts when another
   `train.py --stage pretrain` process is running (`--allow-concurrent` overrides
   only after re-reading `docs/agents/execution-environment.md`).

## Profile run on the real data (2026-11-02, commit at the time of writing)

Both units ran through `run.sh profile` inside
`systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G`; the record is
[`output/profile/profile_record.json`](output/profile/profile_record.json). The
machine-readable record is authoritative; the table is its summary.

| unit | epochs | train / valid rows | wall | train step | train epoch | valid epoch | epoch total | host peak | stop |
|---|---|---|---|---|---|---|---|---|---|
| bounded `--smoke` (2 ep) | 2 | 3,000 / 1,000 | 27.2 s | 0.04 s | 0.25 s | 5.9 s | 6.1 s | **2,703 MiB** | max_epochs |
| one full-width epoch | 1 | **50,520 / 8,524** | **75.9 s** | 0.03 s | **2.73 s** | **48.2 s** | **50.9 s** | **11,273 MiB** | max_epochs |

Both runs used the derived memmap (`cache-mode require`, `verify=metadata`,
`sample_max_abs_diff: 0`), the frozen split and the 649 matched-control groups,
and re-verified the 8 required frozen artifacts (16.4 GB) inside the wall time
(the frozen manifest hashes 13 artifacts; the trainer verifies its 8 inputs). The
one-epoch unit is **full width**: all 50,520 train and all 8,524 valid rows,
9,900+ conditions — it is the first full-width residency measurement of this
pretrain (ticket #4 deferred it to #5).

GPU peak: **4,003 MiB device-level** (3,889 for the one-epoch unit), per-process
attribution unavailable — this worktree's WSL2 driver reports `[N/A]` for
`--query-compute-apps used_memory`, and the window overlapped ticket #6's CM
finetunes, so the device peak includes them (the record does **not** claim the
`device − baseline` delta, 1,716 / 877 MiB, as the pretrain footprint). The
pretrain's isolated GPU footprint measured by ticket #4 is 2,945 MiB; a formal
run must be alone on the GPU, so its device peak is then attributable.
Host peak: **11,273 MiB** (`ru_maxrss` from the run's own manifest, the
authoritative number) — about 2.7 GiB base + the resident portion of the 4.87 GB
memmap + the validation metric's `(8,110 × 18,164)` float64 accumulators
(~1.2 GB each for `sum_pred`/`sum_true`). 18 GiB (`run.sh formal`'s default cap)
is the right budget.

Determinism cross-check: the profile's 2-epoch smoke checkpoint sha256
`df7f9635…` is **identical** to ticket #4's Phase-B smoke pretrain checkpoint
(`df7f9635…`, same seed / subset / config / frozen input) — the pretrain is
reproducible across worktrees.

### Formal budget implication

Measured full-width epoch = **50.9 s** (2.73 s train + 48.2 s validation
AUPRC@200/aux metrics). A 100-epoch run without early stop is therefore
≈ **85 min per seed**, ≈ 4.2 h for the three seeds; early stop cuts it. This
corrects the 06.4 README's 2.8 min/epoch extrapolation (which assumed the
Phase-B cold-cache 1.23 s/step; the warm measured step here is 0.03 s, and
48 of the 51 s are the validation metric, not training). Treat 06.4's 4.5 h/seed
as an upper bound and re-measure the first formal epoch — if the 4.87 GB memmap is
cold, the first epoch pays the random-read cost once.

The single-epoch validation AUPRC@200 in the profile (0.0988) is **not
comparable** to any number in ticket #4's smoke (subset of 1,000 valid rows) and
is not a model claim: one epoch, no selection.

## The record (per seed and per sweep)

`record.py` writes `output/<mode>/…json` + `…_curve.csv`; `check_record.py`
re-validates both against the frozen input. Per **unit** (one seed):

- **identity**: `unit_id`, seed, budget, the exact command, run dir, exit code,
  status, failures (with a log tail on failure);
- **config**: the resolved config block, `config_sha256`, the expected digest
  from `plan.py`, and `config_matches_expected`;
- **input**: the frozen artifact sha256s, gene axis (18,164 genes), matched
  control groups (649), the run's own data counts (`n_train` 50,520,
  `n_valid` 8,524, `control_groups` 649), the derived-cache identity and the
  artifact hash-verification record;
- **checkpoint**: path + sha256 + `reload_identical` + best epoch/AUPRC
  (the same seed's checkpoint is what ticket #7 warm-starts from);
- **validation curve**: the full per-epoch history (train/valid NLL, direction
  AUPRC@200 overall/up/down, Δ RMSE/Pearson/Spearman, LR, train/valid seconds),
  also flattened into the CSV;
- **timings**: pipeline wall, train-loop seconds, mean/median/min/max epoch
  seconds, train/valid split and steps per epoch;
- **resources**: host peak RSS (pipeline manifest = authoritative, plus the
  sampler cross-check), GPU peak attributed/device/baseline/after;
- **stop reason** (`max_epochs`, `early_stop_no_auprc_improvement_20`,
  `no_finite_auprc`) and `epochs_run`.

Per **sweep**: mode (`profile` / `formal`), `complete`, the frozen and requested
seeds, `seeds_completed`/`seeds_failed`/`seeds_missing`, the locked reference
budget, the locked-config sha256, the input identity and `no_seed_cherry_picking`.
A formal record must contain exactly the three frozen seeds; a partial sweep is
written as `tahoe_pretrain_record_partial.json` and the checker rejects it as
formal evidence.

## Execution rules

- One pretrain at a time (the wrapper refuses a concurrent `--stage pretrain`);
  `run.sh` caps the sweep with `systemd-run … MemoryMax=18G`. Every seed's peak is
  recorded; the 05 bring-up measured 15.1 GiB on the raw CSR, this pipeline
  measured **11.3 GiB** on the memmap.
- Never `anndata.read_h5ad(backed="r")` on the frozen view; the wrapper only
  reads the 06.3 manifest/locked config/contract and the derived caches, and the
  06.4 loader is header/`h5py`-based.
- Two heavy jobs repo-wide at most; when one other job runs, keep
  `OMP/MKL/OPENBLAS/NUMEXPR=2` (`--threads 2`, the wrapper default).
- A checkpoint stays local (`*.pt` is gitignored); the record carries its path
  and sha256. The Lead-run sweep's tracked evidence is the record, curves and
  model cards.

## Boundaries

- This substep does not select a seed or a checkpoint for the study: all three
  are reported and all three are reused by ticket #7, as ADR 0032 requires.
- No CM finetune, no cold-drug claim, no transfer conclusion here: that is
  tickets #6–#9.
- The frozen 06.3 input and both derived caches are read-only; the 06.4 pipeline
  is imported, never edited. A negative result remains a valid outcome.
