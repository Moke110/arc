# 06.4_pipeline — reusable PRnet formal-transfer pipeline (ticket #4)

**Status:** Phase A (the pipeline) is built and verified on the synthetic fixture;
**Phase B — the real-data GPU smoke of both initialisations — is executed and green**
(2026-11-02 15:04–15:09 Asia/Shanghai, commit `99650d8`, frozen input freshly
hash-verified, 120 tests).  Pretrain → CM finetune (random **and** Tahoe-initialised) →
predict Δ → evaluate direction AUPRC@200 → model card all ran on the real frozen input;
the checkpoints save+reload, inference is byte-identical across repeats, and the
evaluator's c\* baseline guard passes.  No science result exists here: the smoke is
2 epochs on a small Tahoe subset — pipeline evidence, not a model result.  Ticket #3's
frozen input (`06.3_input/output/`) is the only data this pipeline trains on.

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
│   └── tests/                      # stdlib unittest (120 tests)
└── output/
    ├── smoke/                      # Phase-B smoke record + model cards (force-added)
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

## Phase-B real-data smoke (ticket #4 acceptance)

`execution/smoke_real.py` (via `make gpu-smoke`) ran on the frozen input on
2026-11-02 15:04–15:09 (Asia/Shanghai) at commit `99650d8`, as the only heavy job
on the VM, inside `systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G`,
with a 1 Hz `nvidia-smi` sampler and a host-RSS sampler.  The frozen input passed
fresh hash verification (8 artifacts, 15.30 GB, hash-cache hits 0).  All 17 steps
passed and wrote [`output/smoke/smoke_record.json`](output/smoke/smoke_record.json)
+ `smoke_summary.csv` + the three model cards.  The command was
`python execution/smoke_real.py --data-dir <06.3_input/output> --tahoe-cache
<output/cache/tahoe_log1p_cp10k_v1.npy> --table-cache <output/input_cache>
--warm-rerun --concurrency-probe`; the machine-readable record is authoritative for
every number below.

| step | rows | wall | per-epoch train / valid | peak RSS |
|---|---|---|---|---|
| pretrain (2 ep, subset) | 3,000 train / 1,000 valid | 16.0 s | 7.48 / 7.38 s; valid 5.72 / 5.69 s | 2.70 GiB |
| pretrain warm-cache (1 ep) | 3,000 / 1,000 | 7.4 s | 7.39 s; valid 5.68 s | 2.70 GiB |
| 2-concurrent probe | 2 × (1 ep, 3,000 rows) | 10.2 s | both exit 0 | 2.70 + 2.71 GiB |
| CM fold 0, random init (2 ep) | 571 train / 188 valid / 192 test | 9.3 s | 6.65 s (one-off) then 0.64 s; valid 0.41 s | 1.33 GiB |
| CM fold 0, Tahoe warm start (2 ep) | same | 9.5 s | 6.59 s then 0.78 s; valid 0.40 s | 1.33 GiB |
| predict CM valid / test | 188 / 192 | 1.6 s each | — | ~0.5 GiB |
| predict Tahoe valid | 8,524 | 4.9 s | — | (parent) |
| evaluate CM valid / test | 188 / 192 | 4.0–5.2 s each | — | 0.76 GiB |
| evaluate Tahoe valid | 8,524 (8,110 conditions) | 156 s | — | **8.54 GiB** |
| model cards | — | < 0.1 s each | — | — |

Peaks: **GPU 2,945 MiB** of 16,376, **host 8.54 GiB** — the Tahoe-valid evaluation,
which is the run's peak; training itself peaks at 2.70 GiB (pretrain) / 1.33 GiB (CM).
The 12 GiB scope cap was never tripped, no swap was used, and total wall was **4:29**
(including the fresh 15.30 GB artifact hash; 246 s of steps).

Checks recorded by the same run:

- **Checkpoint save + reload:** both training manifests record `reload_identical:
  true`; the Tahoe-initialised CM run loads the pretrain checkpoint `df7f9635…` in a
  **separate process** through the strict axis/structure-checked loader.
- **Deterministic inference:** all five `predict.py` runs record
  `repeated_inference_identical: true` (zero noise, `model.eval()`, fixed batch order).
- **Leak boundaries:** the commit's 120-test suite (20 leakage tests) passes — the CM
  training loader refuses held-out rows, poisoning test labels does not move a training
  run, c\* is rebuilt only from its recorded membership, and the prediction path has no
  target access.  Every CM evaluation additionally re-derives the input baseline and
  **fails** unless the prediction was computed against this fold's fixed c\* (max abs
  diff ≤ 1e-4); all four CM evaluations passed that guard.
- **Pipeline numbers, not science:** headline direction AUPRC@200 was 0.0466 on the
  Tahoe valid split (8,110 conditions) and 0.0894 / 0.0894 for the random / Tahoe CM
  test split; CM test Δ regression RMSE was 0.44 / 0.42.  These are 2-epoch smoke
  outputs and make **no model or transfer claim**.

### Real-data throughput (for #5–#7 budgets)

Measured, not extrapolated, on the real frozen input + derived memmap cache:

| quantity | measured |
|---|---|
| pretrain train step (batch 512 × 18,164) | ≈ 1.23 s (6 steps → 7.4 s/epoch at 3,000 rows) |
| pretrain validation at 1,000 rows | 5.7 s, of which 5.6 s is the AUPRC@k + auxiliary Δ metric |
| pretrain epoch, subset scale | 7.4 s cold ≈ 7.4 s warm (218 MB of memmap touched, so no warm-cache gain is visible) |
| pretrain host peak | 2.70 GiB / process |
| pretrain-valid evaluation (full 8,524 rows, 8,110 conditions) | 155 s and **8.54 GiB** — one-off after training; do not run two concurrently |
| CM finetune step (fold 0: 571 train rows) | ≈ 0.32 s from epoch 2 (epoch 1 carries a ~6 s one-off target read) |
| CM finetune host peak | 1.33 GiB |
| CM predict / evaluate (fold 0) | 1.6 s / ~5 s; evaluate peak 0.76 GiB |
| GPU peak, all steps | 2,945 MiB (batch 512 makes training GPU-light) |

Scaling estimate for #5 (clearly an extrapolation, to be replaced by the first real
run): a full pretrain epoch is 99 steps ≈ 2 min of training plus ≈ 45 s of full-valid
metric ≈ **2.8 min/epoch**, so 100 epochs ≈ 4.5 h per seed; a full epoch also touches
all 50,520 rows (≈ 3.7 GB of the memmap), which adds to the 2.7 GiB base RSS.  A formal
CM candidate is 20–100 epochs × ≈ 1 s ≈ **0.5–2 min**.  #5/#7 must re-measure.

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

Real-data GPU numbers are in [Phase-B real-data smoke](#phase-b-real-data-smoke-ticket-4-acceptance)
above.  CPU/offline numbers (from `make test` and the Phase-A checks):

| run | rows | time | notes |
|---|---|---|---|
| `make test` (120 tests) | — | ≈14 s | one 18,164-gene fixture + one 512-gene fixture, CPU; peak RSS ≈1.2 GiB |
| CPU e2e (2-epoch pretrain + 2 × 2-epoch CM finetune + predict + evaluate + poisoning) | 120/60 Tahoe, 107/36 CM | ≈9 s | real 18,164 axis, real AUPRC@200 |
| pinned-env CPU chain (pretrain 1 ep → warm-start finetune 1 ep → predict → evaluate) | 16/8 + CM fold | 0.1–0.4 s/step | `arc-stage2-bench-prnet`, CPU only |

Measured derived-cache numbers (pinned env, pre-built 2026-11-02 from the 04
pretrain view; the frozen copy is re-verified before use — see below):

| cache action | measured |
|---|---|
| build 4.87 GB float32 memmap from the 16.3 GB CSR layer | 41 s wall (17.1 s write pass + 16.3 GB source hash), **peak RSS 9.7 GiB** (dirtied memmap pages) |
| `--verify-cache deep` (full-array sha256) | 3.7 s |
| read 512 random rows (37 MB) | 0.27 s first touch, 0.006 s repeat; process RSS 0.9 GiB |
| sidecar identity | rebuilt 2026-11-02 from the 04 pretrain view, then **rebound to the frozen input**: `source_sha256=00f8d906…` (frozen `tahoe_transfer.h5ad`), `layer_sha256=f39c4cb3…` (full-array checksum), 16-row sample max abs diff = 0; the earlier source stays under `rebound_from`.  Rebind took 18 s wall (the full-array checksum of the frozen layer); `--verify-cache metadata` re-checks the sample and the source hash on every load. |
| pre-fix note | the 9.7 GiB build peak and the 04-view source hash above predate the `h5io` metadata fix; opening the h5ad for metadata no longer materialises layers, and the trainer now records `sample_max_abs_diff: 0.0` with the frozen source hash (see the smoke record). |

The Phase-B record carries, per run: wall time, peak RSS, the per-epoch validation
cost split, the cold/warm epoch pair, and the 2-concurrent probe.

**Operator rule (updated by the Phase-B probe):** two concurrent *smoke-subset*
pretrain seeds fit comfortably — measured 2 × 2.70 GiB host, both exit 0, 10.2 s
wall, GPU 2.9 GiB for the whole run — so the shared memmap works.  The probe did
**not** measure full-width residency: a full epoch touches all 50,520 rows (≈3.7 GB
of the 4.87 GB memmap) per process, so budget ≈6 GiB RSS per seed and rely on the
shared page cache.  #5 should start with **one** seed, record its peak, and only then
open a second slot; the cache build (~9.7 GiB peak, before the `h5io` fix) must run
alone.

## Deviations (D0–D4)

| id | deviation | why |
|---|---|---|
| D0 | Derived, verified caches instead of per-process frozen reads: the Tahoe log1p memmap and the CSV copy of the frozen parquet tables. The frozen input is never written. | the pinned env has no parquet engine, and the 05 bring-up's 14.5 GiB CSR RSS blocks concurrent seeds; both caches record source sha256 + equality checks |
| D1 | The two stages are one trainer with a stage switch and a target-free predictor; `predict.py` writes predicted *expression* and Δ = prediction − input baseline (c* for CM, matched control for Tahoe). | keeps a single training/evaluation code path, and makes "test is touch-once" structural: only `evaluate.py` reads truth Δ |
| D2 | Every epoch also stores per-condition Δ RMSE / Pearson / Spearman and the validation cost split, not just the selection AUPRC. | the ticket asks for the auxiliary metrics and for AUPRC's per-epoch cost to be visible; Spearman is the expensive part |
| D3 | The 05 pretrain checkpoint is **not** reused as an initialisation; a formal pretrain must come from `--stage pretrain` in this pipeline. Tahoe validation uses the frozen plate13/14 split. | 05's checkpoint used raw-count controls against log targets (step README boundary); plate13/14 is not a cold-drug claim |
| D4 | The pretrain stage's prediction/evaluation always cover the **full** frozen valid split (8,524 rows / 8,110 conditions); `--max-train-rows`/`--max-valid-rows` only limit training/selection. | prediction is target-free and selection is validation-only, so scoring the whole valid split is honest; it is also the #5/#9 budget number (155 s, 8.54 GiB per evaluation, one at a time). |

## Real-data fixes found by the Phase-B smoke

The synthetic fixture hid four real-data defects; all are fixed inside this substep
(commits `1a3f41d`, `6c49569`, `99650d8`):

1. **Metadata reads materialised layers.** `h5io.open_obs_var` opened h5ad with
   `anndata.read_h5ad(backed="r")`, which eagerly loads every layer — ~8.1 GB RSS on
   the frozen Tahoe view before training.  It now decodes only the `obs`/`var` element
   groups with `read_elem` and takes the shape from the X header (measured: 0.90 s /
   484 MiB on the 16.2 GB file).  This is the shared metadata loader #8 also uses.
2. **Usability must come from the frozen mask.** Re-deriving Tahoe usability
   re-included the 4,757 treated rows without a QC-pass matched control, and their
   control-group lookup crashed (`KeyError: 'plate3|CVCL_1715'`).  The loader now uses
   the contract column `tahoe_usable` and fails closed when the side files disagree;
   the real counts match ticket #3 exactly (50,520 train / 8,524 valid / 649 groups,
   60,861 pairs cross-checked).  The fixture keeps the derived fallback.
3. **O(n_rows × n_obs) condition keys.** A per-row `astype(str)` inside the condition
   comprehension made `TahoeDataset.prepare` take 133.7 s on the real view; hoisting it
   makes it 0.72 s.
4. **Smoke harness and evaluator.** `--device` was passed to `evaluate.py` (which has
   no such flag), the summary writer mishandled the pretrain checkpoint shape, and the
   evaluator built a dead ~1.2 GB float64 baseline for the pretrain stage (peak 9.92 →
   8.54 GiB).

## Environment

- training/GPU: `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python`
  (py3.9, torch 2.0.1+cu118; **no parquet engine** — use `make frozen-tables`).
- tests/fixture/table conversion: `/home/chang/miniconda3/envs/arc/bin/python`.
- `make test` (offline), `make gpu-smoke` (Phase B), `make tahoe-cache`,
  `make frozen-tables`, `make fixture`, `make check`.

## Tracked vs local

Tracked: `AGENTS.md`, `README.md`, `execution/**` (code, tools, tests,
Makefile, run.sh), and the force-added Phase-B record under `output/smoke/`:
`smoke_record.json`, `smoke_summary.csv` and the three model cards
(`pretrain_seed*/`, `cm_fold0_random/`, `cm_fold0_tahoe/`).
Local/gitignored: `*.pt` weights, `output/cache/*.npy` (4.9 GB),
`output/input_cache/`, the `output/smoke/` run trees (predictions, evaluations,
checkpoints, logs), `hash_cache.json`, `run.log`, `__pycache__`.

## Out of scope

No formal training, no 90-candidate sweep, no test-set model selection, no
capacity/env-compound conclusions: those are #5–#11. A negative result is a valid
outcome. This substep publishes the pipeline and its smoke evidence only.
