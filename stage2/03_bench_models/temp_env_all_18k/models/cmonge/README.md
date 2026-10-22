# CMonge adapter — run card (ticket 15)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.0732** (95% CI 0.0698–0.0767), rank 9/12, `ideal_fit: full`. It is `full` in form; the 50-d OT bottleneck and the 10k-iteration substitute are recorded in `FIDELITY.md`. The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/cmonge/` (vendored, commit `195f5cd1`, MIT). No file
  under `pert_model/` is modified or generated (`PYTHONDONTWRITEBYTECODE=1` and
  `sys.dont_write_bytecode = True`).
- **Environment:** `stage2/_envs/cmonge` — a `cp -a` copy of the
  setup bringup venv `pert_model/cmonge/.venv` (Python 3.11.15, jax/jaxlib
  `0.4.38`, flax 0.10.4, optax 0.2.5, ott-jax 0.5.2, rdkit 2023.09.6,
  scanpy 1.11.5). JAX backend: **CPU** (`jax.default_backend() == "cpu"`); the
  CUDA 12 jaxlib install did not complete inside the 1 h environment budget, so
  the row is recorded as a `reproducibility_break` in
  [`FIDELITY.md`](FIDELITY.md). A full 50-step AE + 10,000-iteration fit measured
  32 min on fold 0 (CPU, shared machine), far inside the 6 h/fold cap. Lock:
  `envs/cmonge.lock`.
- **Config:** [`config.json`](config.json) — AE `data_dim=18164`, latent 50,
  2x512; conditional Monge map 4x64, 50-d OT latent, **10,000 iterations** (the
  `conditional-monge-sciplex` setting, not the 500,000 all-drug setting);
  Sinkhorn `epsilon_fitting=1` + Monge gap `epsilon_regularizer=1e-2`; AdamW
  1e-4 / wd 1e-5; batch 512; 200 RDKit props + log-dose.
- **Code:** `stage2/03_bench_models/temp_env_all_18k/scripts/adapters/cmonge/run.py`; shared plumbing in
  `adapters/common.py`; model code is imported from `pert_model/cmonge`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (fold 0, shared GPU/queue lock)
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock /tmp/opencode/stage2_bench_gpu.lock \
  bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model cmonge --fold 0 --timeout 1800

# full recipe (ticket 16 runs folds 0..4)
bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model cmonge --fold 0 --timeout 21600
```

Artifacts land in `stage2/03_bench_models/temp_env_all_18k/output/cmonge/fold_0/`; the fold
summary is appended to `results/leaderboard.csv`.

## Smoke result (fold 0, 4 AE steps + 20 iterations)

| quantity | value |
|---|---|
| fit wall clock | 80.7 s |
| predict wall clock (valid / test) | 4.2 s / 4.7 s |
| total unit wall clock | 105.0 s |
| peak RSS | 2.2 GB |
| best valid author loss (Sinkhorn + Monge gap) | 0.5070 |
| test AUPRC k=200 (headline) | 0.0760 |
| test AUPRC k=100 / k=500 | 0.0561 / 0.1100 |
| test PCC / Spearman | 0.0199 / 0.1633 |
| test R² / RMSE | −16.64 / 0.4603 |
| checkpoint | 77 MB (AE 19,195,686 params + map 31,745 params) |
| test groups | 276 group means over 92 test drugs × 3 doses, 837 wells |

All 837 test wells produced finite 18,164-vectors and were averaged to the 276
test groups. Four AE steps and 20 Monge iterations are far from the author
schedule (50 AE steps + 10,000 iterations); the score is a wiring proof, not
evidence about CMonge.

## Full-recipe probe (fold 0, not a benchmark unit)

A direct adapter run with the real budget (50 AE steps + 10,000 Monge
iterations + ten 1,000-iteration valid evaluations) finished in **32 min** on
CPU on this shared machine (fitting only; no test data was read). The cold-valid
author loss improved from 0.358 at 1,000 iterations to a best of **0.3287 at
iteration 6,000** and ended at 0.3299 at 10,000 — the checkpoint grid keeps the
6,000-iteration state, which is exactly why ticket 09 asks for grid selection
instead of the endpoint. Five folds should take ~2.5–3 h, inside the 16 h/model
cap.

## How the adapter works

- **Stage 1 (AE).** The author's `AutoEncoder` (2x512, latent 50, GELU, MSE,
  AdamW 1e-4/wd 1e-5, warmup-cosine 0→1e-4→1e-5, batch 512, 50 author
  "epochs" = 50 optimizer steps) is fitted on the fold's 2,523 training treated
  wells + all 192 shared DMSO cells, then frozen.
- **Condition.** For each (drug, dose) condition: 200 RDKit 2D descriptors
  (`cmonge.models.rdkit.rdkit_feats`), first column dropped, 15 zero-variance
  columns dropped on the fold's 276 training drugs, z-scored with train-only
  mean/std, then `log(dose_nM)` appended → 186-d context. The network splits it
  into the author's drug/dose modality chunks (`context_entity_bonds`).
- **Stage 2 (CMonge).** One global map; each iteration samples a training
  condition uniformly and 512 source/target cells with replacement in the
  frozen 50-d latent space, minimising Sinkhorn divergence + Monge gap exactly
  as `cmonge.metrics` and `ConditionalMongeTrainer` do.
- **Selection.** Every 1,000 iterations the same author loss is evaluated on
  the 276 cold valid drug-dose conditions (valid expression only, never test);
  the best state is kept. Checkpoints are written every grid point, so a
  SIGTERM'd fit resumes from its last state.
- **Prediction.** One matched DMSO cell per test well
  (`Controls.cell_for(control_group, replicate)`) is encoded, transported under
  the test drug+dose condition, decoded to 18,164 genes; Δ = decoded absolute
  `log1p_cp10k` − the matched plate DMSO mean. No test expression enters
  training or selection.

## Caveats

- The AE is intentionally the author's 50-step budget, so it reconstructs
  poorly and the decoded Δ carries a mean offset (smoke Δ mean ≈ −0.23); this
  lowers R²/RMSE but not the direction metrics. It is a property of the
  two-stage recipe at this budget, not an adapter bug.
- The 50-d OT bottleneck is the model: 18,164 genes are summarised through 50
  dimensions. This is the largest fidelity caveat, not an adapter bug.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 1936 | 1915.6 | 4.1 | 3.9 | 2738 | 10000 | 5000 |
| 1 | completed | 1970 | 1947.8 | 4.1 | 3.9 | 2744 | 10000 | 2000 |
| 2 | completed | 1935 | 1914.0 | 4.0 | 5.1 | 2717 | 10000 | 2000 |
| 3 | completed | 2134 | 2114.6 | 3.8 | 3.7 | 2582 | 10000 | 2000 |
| 4 | completed | 2090 | 2069.6 | 3.9 | 3.7 | 2534 | 10000 | 3000 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
