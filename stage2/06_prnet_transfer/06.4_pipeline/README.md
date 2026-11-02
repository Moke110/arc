# 06.4_pipeline — reusable PRnet formal-transfer pipeline (ticket #4)

**Status:** built and verified (2026-11-02) on the synthetic frozen-input fixture
and on a pinned-env CPU run of both stages; the **real-data GPU smoke is Phase B
and not yet executed**. No science result exists here: a smoke run is a pipeline
check with 2 epochs and small subsets, not a model result. Ticket #3's frozen
input (`06.3_input/output/`) is the only data this pipeline trains on.

This substep is the **single place** the study's training and evaluation code
lives (ADR 0032; step `AGENTS.md`). Tickets #5–#11 invoke it through the CLIs
below and record their own results; #10/#11 only pass `--structure` /
`--train-drugs`, so no pipeline edits are needed.

## What the pipeline does

Two stages, one trainer:

| stage | data | model input | target | selection |
|---|---|---|---|---|
| `pretrain` | frozen Tahoe view (`tahoe_transfer.h5ad`), split = frozen `tahoe_split` (plate13/14 valid) | matched `(plate, cell_line)` DMSO_TF control mean (log1p) | treated `log1p_cp10k` | validation direction AUPRC@200 (**max**) |
| `finetune` | one CM fold (`cm_transfer.h5ad`, `folds.parquet`) | **train**: same-plate matched DMSO control mean from the fold's consensus membership; **validation/test/inference**: the fold's fixed consensus `c*` | treated `log1p_cp10k` | validation direction AUPRC@200 (**max**) |

Both stages optimise the author Gaussian NLL (`GUSS`; may be negative, never
clamped), use Adam, `ReduceLROnPlateau(mode="max")` on validation AUPRC@k, and
save the identical checkpoint payload
`{model, state_dict, config, seed, n_genes, best_valid_auprc, epoch}`. The
experimental truth Δ is computed only by the evaluator
(CM: treated − per-plate matched DMSO mean; Tahoe: treated − matched control).
The model never receives a held-out control pool, and `predict.py` cannot read
treated expression at all.

Metric (frozen): within a condition (CM `drug×dose`; Tahoe
`drug×dose×cell_line`) aggregate replicate predictions/Δ **first**, then per
direction (`top-200` / `bottom-200` true-Δ genes) AP, macro over directions, then
macro over conditions. Δ RMSE / Pearson / Spearman are stored as auxiliary.

## Layout

```
06.4_pipeline/
├── AGENTS.md                       # routing (≤30 lines)
├── README.md                       # this file
├── execution/
│   ├── train.py predict.py evaluate.py run_queue.py make_model_card.py
│   ├── smoke_real.py               # `make gpu-smoke` (Phase-B real-data smoke)
│   ├── run.sh Makefile
│   ├── ptlib/                      # the library (#5–#11 import this, not each other)
│   │   ├── config.py data.py model.py train.py predict.py metrics.py
│   │   ├── evaluate.py manifest.py                       # required modules
│   │   ├── cache.py tables.py checkpoint.py queue.py
│   │   ├── modelcard.py runtime.py h5io.py paths.py
│   ├── tools/{fixture.py,make_fixture.py,build_cache.py,convert_frozen_tables.py}
│   └── tests/                      # stdlib unittest (107 tests)
└── output/
    ├── smoke/                      # Phase-B smoke record (small files force-added)
    ├── input_cache/                # derived CSV copy of the frozen parquet tables
    └── cache/                      # derived Tahoe log1p-CP10k memmap (4.9 GB, local)
```

## CLI surface

All commands run from the repository root (or via `execution/run.sh <cmd> ...`).
`--out` is always an explicit run directory; each run writes `manifest.json`,
`config.json`, `fit_history.json`, `summary_table.csv`, `run.log`,
`checkpoint/model.pt`, `model_card.md`.

```bash
SUB=stage2/06_prnet_transfer/06.4_pipeline/execution
GPY=/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
PY=/home/chang/miniconda3/envs/arc/bin/python

# contract check (hash-verified input, resolved config, architecture; no training)
$GPY $SUB/train.py --stage finetune --fold 0 --out /tmp/x --check

# pretrain (Tahoe) / finetune (CM fold)
$GPY $SUB/train.py --stage pretrain --out <run> --seed 20261021 --device cuda
$GPY $SUB/train.py --stage finetune --fold 0 --out <run> --seed 20261021 \
     --lr 3e-4 --init random --device cuda
$GPY $SUB/train.py --stage finetune --fold 0 --out <run> \
     --init pretrained --init-checkpoint <pretrain-run>/checkpoint/model.pt

# prediction (target-free, byte-identical across repeats) and evaluation
$GPY $SUB/predict.py --stage finetune --fold 0 --role test \
     --checkpoint <run>/checkpoint/model.pt --out <pred-dir>
$GPY $SUB/evaluate.py --stage finetune --fold 0 --role test \
     --pred-dir <pred-dir> --out <eval-dir>

# resume-safe parallel queue
$GPY $SUB/run_queue.py --queue <units.json> --out-root <runs> --workers 2

# regenerate a model card (optionally with an evaluation attached)
$GPY $SUB/make_model_card.py --run-dir <run> --eval-dir <eval-dir>
```

Key `train.py` flags: `--stage`, `--fold`, `--seed`, `--device`
(`cuda` is a hard failure when unavailable — no silent CPU fallback),
`--structure hidden128|latent128|hidden256|<json>`, `--init random|pretrained`,
`--init-checkpoint`, `--train-drugs main|expanded`, `--epochs`, `--lr`,
`--batch-size`, `--min-epochs`, `--early-stop-patience`, `--sched-patience`,
`--max-train-rows`, `--max-valid-rows` (Tahoe subset), `--tahoe-cache`,
`--cache-mode auto|require|build|off`, `--verify-cache metadata|deep`,
`--hash-cache`, `--table-cache`, `--smoke`, `--check`, `--no-model-card`.

### Config keys

The frozen `06.3_input/output/locked_config.json` is normalised and validated
(aliases accepted and recorded in the run manifest; missing recipe/structure
values fall back to the documented frozen defaults, also recorded). Canonical
blocks:

| block | keys |
|---|---|
| `structure` | `n_genes`, `hidden` (`[128]`), `z_dim` (64), `c_dim` (64), `n_dim` (10), `drug_dim` (1024), `adaptor` (`[128]`), `dropout` (0.05) |
| `recipe` | `batch_size` (512), `max_epochs` (100), `min_epochs` (20), `early_stop_patience` (20), `sched_factor` (0.5), `sched_patience` (10), `min_lr` (1e-8), `weight_decay` (1e-8), `lr` (resolved per stage) |
| `metric` | `name='direction_auprc'`, `k=200`, `plateau_mode='max'` |
| top level | `loss='gaussian_nll'`, `optimizer='adam'`, `prediction_noise='zero'`, `dose_encoding` (FCFP4 r2, 1024 bits, ×log10(dose_µM+1)), `n_folds=5`, `seeds`, `lrs.{pretrain,cm}` |

`--structure latent128` → `z_dim=128` (7,201,448 params); `hidden256` →
`hidden=[256]` (14,194,280 params). `--train-drugs expanded` adds the
eligibility-excluded QC-pass wells to **training only**; they never get a fold or
a score.

### Run manifest (every run)

`schema=ptlib.run_manifest/v1`, `kind=train|predict|evaluate`, git SHA/branch/dirty,
`data.dir` + every frozen artifact sha256, locked-config sha256 + alias/default
mapping, resolved config + `config_sha256`, seed, structure, device (torch/CUDA,
cudnn deterministic), initialisation (mode, source path, checkpoint sha256),
checkpoint sha256 + strict-reload check, noise policy, wall time, **per-epoch
validation cost split (gather/forward/AUPRC)**, peak RSS, epochs run and stop
reason. `predict` additionally records row count, prediction sha256 and
`repeated_inference_identical`.

## How the other tickets invoke it

Queue JSON: `{"units":[{"unit_id","stage","args":[...]}]}` (one unit per training
run; `--out`/`--data-dir` are owned by the queue). `run_queue.py` skips a unit
whose manifest is complete **and** whose `config_sha256` matches the unit's
resolved args, and re-runs changed/tampered units.

```bash
# #5 — 3 seeded Tahoe pretrains (once: make tahoe-cache; sequential/two-wide)
#   units: {"unit_id":"tahoe_seed20261021","stage":"pretrain",
#           "args":["--seed","20261021","--device","cuda","--cache-mode","require",
#                   "--tahoe-cache",CACHE,"--hash-cache",HASH]}
$GPY $SUB/run_queue.py --queue 06.5_units.json \
    --out-root <repo>/stage2/06_prnet_transfer/06.5_tahoe_pretrain/output/runs --workers 2

# #6 — 45 random-init CM candidates: --fold K --seed S --lr {1e-4,3e-4,1e-3}
#      --init random --device cuda  (queue --workers <= 4, --omp-threads 2)
# #7 — 45 Tahoe-init CM candidates: same, but --init pretrained
#      --init-checkpoint <#5 run>/checkpoint/model.pt   (strict, axis+structure checked)
# #8 — baselines: write prediction dirs in the predict.py format (pred_expr.npy,
#      delta.npy, meta.csv with a `row` column, manifest.json) and score them with
#      evaluate.py; the evaluator is the shared truth/metric path.
# #9 — for each of the 30 selected checkpoints: predict.py --role test, then
#      evaluate.py; aggregate across folds/seeds/units and regenerate cards.
# #10 — own pretrain with --structure latent128|hidden256, then finetune as #7.
# #11 — finetune with --train-drugs expanded (test drugs stay excluded).
```

Evaluation artifacts (`<eval-dir>/`): `metrics.csv` (macro rows per k/direction/
metric), `per_condition.csv`, `aggregate.json`, `deltas.npz`
(`conditions`, `true_delta`, `pred_delta`, `rows`), `manifest.json`. #9 consumes
these directly (per-drug pairing/bootstrapping is #9's job, not this pipeline's).

## Execution rules (reproducibility)

- Seed python/numpy/torch from `--seed`; cuDNN deterministic, benchmark off.
- Training noise: per-batch `torch.randn` from a CPU generator seeded with the run
  seed. Validation/inference: **zero noise**, `model.eval()`, fixed row order →
  deterministic; `predict.py` repeats inference (`--repeats 2`) and **fails**
  unless the arrays are byte-identical.
- Selection/early stop on validation direction AUPRC@200 (max semantics); the
  test split is only read by `evaluate.py`, after a checkpoint exists.
- `pretrain` never sees CM data; `finetune` never sees Tahoe test/valid rows.

## Derived caches (not frozen evidence)

1. **Tahoe log1p memmap** — `make tahoe-cache` builds a dense float32 `.npy`
   (67,018 × 18,164 = 4.87 GB) from the frozen CSR layer plus a sidecar with the
   source parquet-h5ad sha256, the full-array checksum and a 16-row equality
   check against the frozen layer. Training reads it through a shared page-cache
   memmap, so concurrent pretrain seeds no longer each hold a private 14.5 GiB
   CSR copy. Verified on load (`metadata`; `deep` re-hashes the whole file).
2. **Frozen-table CSV cache** — the pinned GPU env has **no parquet engine**, so
   `make frozen-tables` (run with the `arc` env) converts `folds`,
   `consensus_membership` and `tahoe_control_groups` parquet to CSV under
   `output/input_cache/`, recording each source parquet's frozen sha256. The
   loader refuses a stale entry and fails closed with the fix command when the
   cache is absent. The frozen input and the pinned environment are not modified.

## Timings / memory

Measured so far (CPU, `arc` env, synthetic fixture; from `make test` and the
Phase-A checks). **Phase-B real-GPU smoke numbers are pending** and will be added
to [`output/smoke/`](output/smoke/README.md) with the model cards.

| run | rows | time | notes |
|---|---|---|---|
| `make test` (107 tests) | — | ≈15 s | one 18,164-gene fixture + one 512-gene fixture, CPU |
| CPU e2e (2-epoch pretrain + 2 × 2-epoch CM finetune + predict + evaluate + poisoning) | 120/60 Tahoe, 107/36 CM | ≈9 s | real 18,164 axis, real AUPRC@200 |
| pinned-env CPU chain (pretrain 1 ep → warm-start finetune 1 ep → predict → evaluate) | 16/8 + CM fold | 0.1–0.4 s/step | `arc-stage2-bench-prnet`, CPU only |

Phase-B smoke will record, per run: wall time, peak RSS, per-epoch validation cost
split, cold vs warm cache epoch times (2 epochs = one cold + one warm), and the
measured concurrency rule for pretrain seeds (2-concurrent probe). Until those
numbers exist this README makes **no** concurrency or throughput claim.

## Deviations (D0–D3)

| id | deviation | why |
|---|---|---|
| D0 | Derived, verified caches instead of per-process frozen reads: the Tahoe log1p memmap and the CSV copy of the frozen parquet tables. The frozen input is never written. | the pinned env has no parquet engine, and the 05 bring-up's 14.5 GiB CSR RSS blocks concurrent seeds; both caches record source sha256 + equality checks |
| D1 | The two stages are one trainer with a stage switch and a target-free predictor; `predict.py` writes predicted *expression* and Δ = prediction − input baseline (c* for CM, matched control for Tahoe). | keeps a single training/evaluation code path, and makes "test is touch-once" structural: only `evaluate.py` reads truth Δ |
| D2 | Every epoch also stores per-condition Δ RMSE / Pearson / Spearman and the validation cost split, not just the selection AUPRC. | the ticket asks for the auxiliary metrics and for AUPRC's per-epoch cost to be visible; Spearman is the expensive part |
| D3 | The 05 pretrain checkpoint is **not** reused as an initialisation; a formal pretrain must come from `--stage pretrain` in this pipeline. Tahoe validation uses the frozen plate13/14 split. | 05's checkpoint used raw-count controls against log targets (step README boundary); plate13/14 is not a cold-drug claim |

## Environment

- training/GPU: `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python`
  (py3.9, torch 2.0.1+cu118; **no parquet engine** — use `make frozen-tables`).
- tests/fixture/table conversion: `/home/chang/miniconda3/envs/arc/bin/python`.
- `make test` (offline), `make gpu-smoke` (Phase B), `make tahoe-cache`,
  `make frozen-tables`, `make fixture`, `make check`.

## Tracked vs local

Tracked: `AGENTS.md`, `README.md`, `execution/**` (code, tools, tests,
Makefile, run.sh), `output/smoke/` record files and model cards (force-added).
Local/gitignored: `*.pt` weights, `output/cache/*.npy` (4.9 GB),
`output/input_cache/`, `output/smoke/` raw predictions/checkpoints/logs,
`run.log`, `__pycache__`.

## Out of scope

No formal training, no 90-candidate sweep, no test-set model selection, no
capacity/env-compound conclusions: those are #5–#11. A negative result is a valid
outcome. This substep publishes the pipeline and its smoke evidence only.
