# 05_prnet_pretrain — PRnet pretrain bring-up on the Tahoe-100M view (GPU)

Stage 2 has to answer whether a drug-perturbation model can be **pretrained on
Tahoe-100M** before it is finetuned on the small CM corpus. This step is the
pipeline bring-up for that question: train PRnet's author model (`PGM`) on the
model-facing Tahoe view until a loss-epoch curve exists on real GPU hardware.

It is deliberately narrow — **no CM finetune, no warm start, no cold-drug
evaluation, no comparison against `03_models_cmbench`**. The decision is
[ADR 0031](../../docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md) §2
(contract as planned in
[`.scratch/stage2-05-prnet-pretrain/spec.md`](../../.scratch/stage2-05-prnet-pretrain/spec.md)
§2). Since the 02+04 merge
([ADR 0033](../../docs/adr/0033-stage2-dataset-corpus.md)) this step also **owns**
the CM-schema view (moved byte-identically) and its builder (moved, then
repointed at the merged corpus and given a `--verify` replay), both inherited
from the retired `stage2/04_tahoe_data/`; the view's deviations (D0–D3) are
documented below and mirrored into `output/model_card.md`.

## Status

| | |
|---|---|
| execution | **EXECUTED 2026-11-01** (fresh user go-ahead) — after the inherited view passed QA 9/9 |
| pretrain run | **68 epochs**, early stop fired at epoch 68 (patience 21 > 20 once the 20-epoch floor was met); best valid **−2.628080** at epoch 47; `NVIDIA GeForce RTX 4070 Ti SUPER` (torch 2.0.1+cu118) |
| data | `stage2/05_prnet_pretrain/output/pretrain_view/tahoe_cm_schema.h5ad` (moved byte-identically from `04_tahoe_data`; `content_sha256` `476812c0…`) — 67,018 × 18,164; 65,318 samples (**55,922 train / 9,396 valid**, split by plate: `plate13` + `plate14` valid); **300** treated rows dropped (no resolvable structure); 700 control groups; 375 unique SMILES; 1,400 controls / 65,618 treated |
| cost | **956 s (15.9 min)** wall between the run's own start/finish stamps (`17:01:45`→`17:17:41`), over 68 epochs; the in-process monotonic timer reads 883.2 s for the same span — it runs low against the wall clock on this host (see `model_card.md`); peak host RSS **15,150 MiB** on a 31 GiB host |
| loss | Gaussian NLL is **negative throughout** (train −0.714 → −2.814; best valid −2.628). It is a log-likelihood, so its sign carries no meaning; it is never clamped. It says nothing about transfer, and it is **not comparable** to any `03_models_cmbench` figure |
| pipeline | **smoke-verified on CUDA** (2 epochs on a synthetic fixture, ticket 04) — retained as the pipeline proof |
| fixture encoding | carries the production view's `dataframe/0.2.0` `uns['control_mean']` (generated with the `arc` env) |

The plan of record is
[`.scratch/stage2-05-prnet-pretrain/`](../../.scratch/stage2-05-prnet-pretrain/)
(`spec.md` contract, `map.md` decisions, `issues/` tickets — all six resolved).
The accepted decision record is
[`docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md`](../../docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md).

Real run (2026-11-01). `execution/run.sh` on the inherited view: 68 epochs in 956 s (15.9 min
between the start and finish stamps), with early
stopping firing at epoch 68 after 21 stale epochs; the valid Gaussian NLL fell from
−1.359030 at epoch 1 to a best **−2.628080** at epoch 47 and was still the best when the
run stopped, and the plateau scheduler halved the learning rate to 5e-4 at epoch 58. The
curve (`output/summary.png`) is drawn by the run itself: integer epoch axis, train +
valid series, best-epoch star at 47, early-stop dotted line at 68 (drawn inside the
panel). `output/model_card.md`
records the recipe, the four deviations, the row counts, the device, the wall time and
the peak memory. This is a pretraining **pipeline** result: no finetune, no cold-drug
evaluation, no claim that pretraining helps, and no comparison with `03_models_cmbench`.

Smoke result (`--smoke`, synthetic 300 × 18,164 fixture, RTX 4070 Ti SUPER):
epoch 1 `train 0.038362 / valid 0.048439`, epoch 2 `train 0.019514 / valid
0.047446`; `fit_history.json`, `summary.png` and a 7,185,000-parameter checkpoint
written. These are **pipeline** numbers from synthetic data — not a result, and not
comparable to any `03_models_cmbench` figure.

## Reproduce

Run from the repository root. Use the pinned GPU environment — **never** the
shared `arc` env (its torch 2.12 dropped `ReduceLROnPlateau(verbose=...)` and its
scanpy/numba cache is broken):

```bash
PY=/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
ARC=/home/chang/miniconda3/envs/arc/bin/python

# 0. contract check: resolved root/view, shape, split, dropped rows, device
$PY stage2/05_prnet_pretrain/execution/train_prnet.py --check

# 1. the real bring-up (100-epoch cap, early stop after 20 stale epochs, min 20)
#    2026-11-01 run: early stop at epoch 68, best valid -2.628080 (see Status)
stage2/05_prnet_pretrain/execution/run.sh

# 2. pipeline check without the real view: build the synthetic fixture, then smoke
#    (the fixture is written with the `arc` env so it carries the production
#     view's exact dataframe/0.2.0 encoding for uns['control_mean']; see Environment)
$ARC stage2/05_prnet_pretrain/execution/_smoke/make_fixture.py
$PY stage2/05_prnet_pretrain/execution/train_prnet.py --smoke \
    --view stage2/05_prnet_pretrain/execution/_smoke/tahoe_cm_schema_smoke.h5ad

# 3. redraw the curve from the recorded history (no training, no view read)
$PY stage2/05_prnet_pretrain/execution/train_prnet.py --replot
```

`--view` / `--output` / `--device` / `--seed` are overridable; defaults are this
step's own view, `stage2/05_prnet_pretrain/output/`, `cuda`, `20261021`.

## View ownership (since the 02+04 merge)

`execution/build_pretrain_view.py` was moved here from the retired
`stage2/04_tahoe_data/`; its path constants were then repointed at the merged
corpus + tables and `--verify` was added, so the builder is **not** a
byte-identical move. The view file itself was **moved, not rebuilt**:
`output/pretrain_view/tahoe_cm_schema.h5ad` is 16,319,543,984 B with
`content_sha256 = 476812c0…`, unchanged (ADR 0033 decision 5). The builder reads
the merged corpus and its identity tables, and writes this step's own
`output/pretrain_view/`:

| input | path |
|---|---|
| Tahoe corpus (obs + row order) | `stage2/02_dataset_corpus/output/corpus/v2/tahoe.h5ad` |
| per-plate CSR corpus parts | `stage2/02_dataset_corpus/execution/_parts/` |
| 18,164-gene CM axis | `stage2/02_dataset_corpus/output/corpus/v2/tables/genes.csv` |
| CM drug table (overlap routes) | `stage2/02_dataset_corpus/output/corpus/v2/tables/compounds.csv` |

A rebuild is not part of this round (it is hours of I/O and must reproduce the
same hash). **`--verify`** is the check that the moved file still equals what the
corpus implies: it replays the projection over the corpus parts in the builder's
own order (two streaming passes, one part in memory, no second 16 GB write) and
compares the recomputed `content_sha256` to the recorded `476812c0…`. It is
**default off**, reads only, and needs the merged corpus + tables; missing inputs
exit with a message (no traceback), and a hash mismatch exits non-zero.

```bash
$ARC stage2/05_prnet_pretrain/execution/build_pretrain_view.py --check   # paths only
$ARC stage2/05_prnet_pretrain/execution/build_pretrain_view.py --verify  # replay + compare
```

## Input contract

One file only: `stage2/05_prnet_pretrain/output/pretrain_view/tahoe_cm_schema.h5ad`
(67,018 × 18,164; schema per spec §1.2–1.4). The trainer never touches the corpus,
the per-plate parts, or CM h5ads, and it never densifies the matrix — it reads obs
and `uns` via anndata **backed** mode and the CSR layer through h5py.

Per training sample (a non-control, structure-bearing row):

| tensor | source | shape |
|---|---|---|
| `control` | `uns['control_mean']` row of the row's `control_group` | (18,164) float32 |
| `fingerprint` | RDKit rFCFP4 of that row's SMILES, × `log10(dose_µM + 1)` | (1,024) float32 |
| `target` | that row's `layers['log1p_cp10k']` | (18,164) float32 |

Rows with `is_control` are **not** samples (they are control inputs). Rows with
`has_structure == False` and not control are dropped, with the count logged.
SMILES fingerprints are cached once per unique SMILES (375 in the real view: 379 drugs,
two of them structure-less) and dose-scaled per row. Split is **by plate**:
`plate13` + `plate14` are validation; validation loss drives the scheduler but never
early-stops before epoch 20.

## Model and training recipe — same-source copy, do not diverge

The upstream tree is read-only, so the classes and the loop are **re-implemented
in-file** as exact equivalents (no import from `pert_model/` or from
`03_models_cmbench/`):

| component | mirrored from | setting |
|---|---|---|
| `PAdaptor` | `pert_model/PRnet/source/models/PRnet.py:238-271` | 1024 rFCFP4 → 128 → 64 |
| `PEncoder` | `.../PRnet.py:144-181` | [18164 + 64] → 128 → 64 (mean only) |
| `PDecoder` | `.../PRnet.py:183-233` | [64 + 64 + 10] → 128 → 2 × 18164; `relu(head[:, :dim])` mean, `softplus(head[:, dim:])` variance |
| drug encoder | `pert_model/PRnet/source/data/_utils.py:44-68` | `GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True, nBits=1024)` × `log10(dose+1)` |
| loss | `.../trainer/PRnetTrainer.py:193-194` (author `GUSS`) | `F.gaussian_nll_loss(mean, target, var, reduction="mean")` |
| optimiser | `.../PRnetTrainer.py:161-165` | Adam(lr 1e-3, wd 1e-8); `ReduceLROnPlateau(mode='min', factor 0.5, patience 10, min_lr 1e-8)` on **valid** loss |
| loop | `stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/prnet/run.py:188-258` | batch 512, max 100 epochs, early stop 20, min epochs 20, seed 20261021, dropout 0.05, `noise_dim=10` |

The hyperparameters, the early-stop rule and the seed are copied from the frozen
03.1 adapter so the **recipe** is identical; the reported losses are still not
commensurable with the benchmark's (different corpus, rows and split). `--check`
asserts the six audited layer shapes before any GPU time is spent.

Per-epoch losses are **batch-size-weighted** (`Σ NLL·n / Σn`), so the final
partial batch cannot bias an epoch; a trailing single-row batch is folded into
the previous batch, so no row is dropped. `run.log` prints the split counts and,
per epoch, the losses, the scheduler lr and the early-stop state (`patience
k/20`, whether `min_epochs` is met) — a finished curve can be explained without
re-reading the code.

**The loss may be negative.** A Gaussian NLL is a log-likelihood, so it is
unbounded below and its sign carries no meaning; the program never clamps,
`abs()`-es or otherwise "fixes" it, and `summary.png` labels the axis
`Gaussian NLL (may be negative)`.

## Outputs (all under `output/`, overwritten on every rerun)

| file | content |
|---|---|
| `run.log` | full log, overwrite mode (never appends) |
| `fit_history.json` | per-epoch train/valid loss, lr, seconds + best/stop metadata, plus the run's start/finish, wall time and peak RSS |
| `summary.png` | loss-epoch curve: train blue, valid orange, best-epoch star, early-stop marker |
| `summary_table.csv` | one row per epoch |
| `config.json` | every hyperparameter, the view, the split, row counts, device |
| `checkpoint/model.pt` | `{model, state_dict, config, seed, n_genes, best_valid_loss, epoch}` |
| `model_card.md` | what was trained, on what, with which deviations, plus wall time and peak memory |

The five human-facing artifacts — `config.json`, `fit_history.json`,
`summary_table.csv`, `summary.png`, `model_card.md` — are force-added to version
control (the repo-wide `output/` rule needs `git add -f`). `run.log` and the 28.7 MB
`checkpoint/model.pt` stay local, per ticket 05: the curve, history, table, config and
model card are the evidence a later reader needs, and the weights are regenerable.

## Deviations

Inherited with the moved view and still recorded here (the builder now lives in
this step):

- **D0 layout** — CSR instead of CM's dense `int32`; the builder never materialises a
  dense 18,164 × 67,018 buffer, so what differs is storage layout, not schema.
- **D1 control keying** — `control_group` = `tahoe::<plate>::<cell_line_id>` (700 groups)
  so the control mean never mixes cell lines.
- **D2 gene axis** — 18,164 slots, `FAM153B` all-zero (no Tahoe measurement).
- **D3 provenance gaps** — `time_h`/`cell_type`/`donor` NaN, `atc`/`moa`/`mesh_classes`
  empty, `pharm_*` False; recorded, not invented.

Two further notes for this step:

- **The control input and the target are not in the same space** — the control is
  the group **mean of raw counts** while the target is `log1p_cp10k`. The author's
  dataset draws both from one matrix, so this is a recorded deviation from parity,
  not a match with it. It does not affect the pipeline claim, but it must be
  resolved before any finetune or transfer claim.
- `uns['control_mean']` is the production view's DataFrame encoding (`dataframe/0.2.0`,
  written by anndata 0.12.16 in the `arc` env). The pinned 0.10.8 reader here
  handles it and materialises a `(700, 18164)` `DataFrame`; `config.json` records
  which encoding was read. The trainer also accepts a transposed DataFrame or a
  bare `(700, 18164)` array plus a label array, as robustness — it never assumes
  the DataFrame is absent.

## Environment

Two pinned environments, split by role:

- **05 trains/reads** with `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python` —
  Python 3.9.25, torch 2.0.1+cu118 (CUDA available), anndata 0.10.8, RDKit
  2025.09.2, matplotlib 3.9.4 (Agg).
- **The view builder** writes the view with `/home/chang/miniconda3/envs/arc/bin/python`
  — anndata 0.12.16 (`execution/build_pretrain_view.py`; it moved here from 04).
  That is why the smoke fixture is generated with `arc`: it must
  carry the identical `dataframe/0.2.0` encoding for `uns['control_mean']` that
  production carries. (anndata 0.10.8 cannot *write* that encoding — HDF5
  overflows — but it *reads* it, which is all this step needs.)

GPU: NVIDIA GeForce RTX 4070 Ti SUPER, 16 GiB, driver 617.42. Host RAM 31 GiB: the
view is loaded as CSR and never densified, but it is **not** small — reading the real
`layers['log1p_cp10k']` (1,010,892,356 nonzeros) leaves **14,562 MiB** host RSS (14.2 GiB),
and the measured run peak is **15,150 MiB** (`output/model_card.md` is the authoritative
figure). Targets are gathered per batch.

## Out of scope this round

CM finetune / warm start, cold-drug re-benchmark, the 5k HVG fork, any
pretraining *protocol* beyond ADR 0031's bring-up (transfer, warm-start or
finetune experiments), any edit inside
`pert_model/`, any regeneration of `03_models_cmbench`, and rebuilding the
`02_dataset_corpus` corpus or this view.
