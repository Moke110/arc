# 05_prnet_pretrain — PRnet pretrain bring-up on the Tahoe-100M view (GPU)

Stage 2 has to answer whether a drug-perturbation model can be **pretrained on
Tahoe-100M** before it is finetuned on the small CM corpus. This step is the
pipeline bring-up for that question: train PRnet's author model (`PGM`) on the
model-facing Tahoe view until a loss-epoch curve exists on real GPU hardware.

It is deliberately narrow — **no CM finetune, no warm start, no cold-drug
evaluation, no comparison against `03_models_cmbench`**. The contract is
[`.scratch/stage2-05-prnet-pretrain/spec.md`](../../.scratch/stage2-05-prnet-pretrain/spec.md)
§2; the view it consumes is built by
[`stage2/04_tahoe_data/`](../04_tahoe_data/README.md) and its deviations are
recorded there (D0–D3, mirrored into `output/model_card.md`).

## Status

| | |
|---|---|
| execution | **HALTED by user decision 2026-11-01** — plan first, do not run |
| pipeline | **smoke-verified on CUDA**, 2 epochs on a synthetic fixture |
| fixture encoding | carries 04's production `dataframe/0.2.0` `uns['control_mean']` (generated with the `arc` env) |
| full 100-epoch run | **NOT EXECUTED** — needs a fresh go-ahead, and 04's view must pass QA first |

The plan of record is
[`.scratch/stage2-05-prnet-pretrain/`](../../.scratch/stage2-05-prnet-pretrain/)
(`spec.md` contract, `map.md` decisions, `issues/` tickets, `adr-0031-draft.md`).
Tickets: **01** build/QA 04's view, **03** the full run, **04** the doc/ADR flip.

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

# 1. the real bring-up (100 epochs, early stop after 20 stale epochs, min 20)
stage2/05_prnet_pretrain/execution/run.sh

# 2. pipeline check without the real view: build the synthetic fixture, then smoke
#    (the fixture is written with the `arc` env so it carries 04's exact
#     dataframe/0.2.0 encoding for uns['control_mean']; see Environment)
$ARC stage2/05_prnet_pretrain/execution/_smoke/make_fixture.py
$PY stage2/05_prnet_pretrain/execution/train_prnet.py --smoke \
    --view stage2/05_prnet_pretrain/execution/_smoke/tahoe_cm_schema_smoke.h5ad
```

`--view` / `--output` / `--device` / `--seed` are overridable; defaults are the 04
view, `stage2/05_prnet_pretrain/output/`, `cuda`, `20261021`.

## Input contract

One file only: `stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad`
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
SMILES fingerprints are cached once per unique SMILES (≈377) and dose-scaled per
row. Split is **by plate**: `plate13` + `plate14` are validation; validation loss
drives the scheduler but never early-stops before epoch 20.

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
03.1 adapter so the pretrain run and the benchmark run stay comparable. `--check`
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
| `fit_history.json` | per-epoch train/valid loss, lr, seconds + best/stop metadata |
| `summary.png` | loss-epoch curve: train blue, valid orange, best-epoch star, early-stop marker |
| `summary_table.csv` | one row per epoch |
| `config.json` | every hyperparameter, the view, the split, row counts, device |
| `checkpoint/model.pt` | `{model, state_dict, config, seed, n_genes, best_valid_loss, epoch}` |
| `model_card.md` | what was trained, on what, with which deviations |

The whole directory is gitignored by the repo-wide `output/` rule; like
`output/pseudobulk/` in 04, these are local artifacts, not tracked evidence.

## Deviations

Inherited from the 04 view and recorded in `model_card.md`:

- **D0 layout** — CSR instead of CM's dense `int32` (dense 18,164 × 67,018 = 4.87 TB);
  schema-identical.
- **D1 control keying** — `control_group` = `tahoe::<plate>::<cell_line_id>` (700 groups)
  so the control mean never mixes cell lines.
- **D2 gene axis** — 18,164 slots, `FAM153B` all-zero (no Tahoe measurement).
- **D3 provenance gaps** — `time_h`/`cell_type`/`donor` NaN, `atc`/`moa`/`mesh_classes`
  empty, `pharm_*` False; recorded, not invented.

Two further notes for this step:

- The control input is the group **mean of raw counts**, not `log1p_cp10k`; this
  mirrors the author's dataset, where control and target rows share a space.
- `uns['control_mean']` is 04's production DataFrame encoding (`dataframe/0.2.0`,
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
- **04 writes** the view with `/home/chang/miniconda3/envs/arc/bin/python` —
  anndata 0.12.16. That is why the smoke fixture is generated with `arc`: it must
  carry the identical `dataframe/0.2.0` encoding for `uns['control_mean']` that
  production carries. (anndata 0.10.8 cannot *write* that encoding — HDF5
  overflows — but it *reads* it, which is all this step needs.)

GPU: NVIDIA GeForce RTX 4070 Ti SUPER, 16 GiB, driver 617.42. Host RAM 31 GiB: the
view is loaded as CSR (~3–4 GB for the real layer) and targets are gathered per
batch — never densified.

## Out of scope this round

CM finetune / warm start, cold-drug re-benchmark, the 5k HVG fork, a full
pretraining-protocol ADR (write it once the real curve exists), any edit inside
`pert_model/`, any regeneration of `03_models_cmbench`, and rebuilding the 04
corpus or its view.
