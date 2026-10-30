# biolord adapter — run card (ticket 14)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.0974** (95% CI 0.0914–0.1040), rank 5/12, `ideal_fit: partial`. It is `partial` because prediction is transductive (registered DMSO control rows). The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/biolord/source/` (vendored, commit `b7688790`, BSD-3-Clause;
  no published weights). No file under `pert_model/` is modified or generated
  (`sys.dont_write_bytecode` is set before the model import).
- **Environment:** `stage2/_envs/biolord` (venv on the `arc`
  interpreter, Python 3.13.13, `--system-site-packages`), overlay pins
  `scvi-tools==1.1.6 anndata==0.10.9 zarr==2.18.7 biolord==0.0.3` from
  `stage2/01_setup/biolord/requirements.txt`, plus `descriptastorus==2.8.0`,
  `pandas_flavor==0.7.0` and `xarray==2024.3.0` for the author's RDKit2D
  features. torch is the `arc` wheel `2.12.0+cu130`; `torch.cuda.is_available()`
  is `True`. The model source is imported from
  `pert_model/biolord/source/src` (the PyPI `biolord==0.0.3` install supplies only the
  distribution metadata the version lookup expects).
- **Config:** [`config.json`](config.json) records the SI Table 4/5 author
  recipe (`rdkit2d_dose` ordered + `cell_type` categorical, τ=1e4, σ=20,
  decoder 4096×4, `n_latent` 256, batch 512, 500 epochs, Adam latent/decoder
  1e-4, attribute NN 2048×2 @ 1e-2, cosine to 1e-5) and the literal
  18,164-gene output.
- **Code:** `stage2/03_bench_models/temp_env_all_18k/scripts/adapters/biolord/run.py`; shared plumbing in
  `adapters/common.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, two batches, fold 0) — holds the shared GPU lock
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock /tmp/opencode/stage2_bench_gpu.lock \
  bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model biolord --fold 0 --timeout 1800
```

Without `BENCH_SMOKE=1` the adapter runs the frozen author budget (500 epochs,
checkpoint grid every 50 epochs). Artifacts land in
`stage2/03_bench_models/temp_env_all_18k/output/biolord/fold_0/`; the fold summary is appended to
`results/leaderboard.csv`.

## Smoke result (fold 0, 1 epoch / 2 batches)

| quantity | value |
|---|---|
| fit wall clock | 17.2 s |
| predict wall clock (valid / test) | 8.3 s / 7.5 s |
| whole fold wall clock | 49 s |
| peak RSS | 5.9 GB |
| registered rows | 2,715 (2,523 treated + 192 DMSO) |
| valid author loss (grid epoch 1) | 2,329.97 |
| test AUPRC k=200 (headline) | 0.0764 |
| test AUPRC k=100 / k=500 | 0.0564 / 0.1108 |
| test PCC / Spearman | 0.0209 / 0.1757 |
| test R² / RMSE | −14.52 / 0.4405 |
| valid AUPRC k=200 / PCC | 0.0776 / 0.0420 |
| checkpoint | 202.9 M params, 811.6 MB state (`best.pt`), plus a 811.6 MB `last.pt` and a 200.3 MB `train_adata.h5ad` |

All 837 test wells produced finite 18,164-vectors and were averaged to the 276
test groups (92 drug keys × 3 doses). Two independent smoke runs produced
byte-identical `pred_test.npy` (seed `20261021`).

A one-epoch model is far from the author schedule (500 epochs); the score is
**not** evidence about biolord.

**Full-recipe check (informal, not a benchmark artifact).** A direct 500-epoch
fit on fold 0 (same seed, run directory under `/tmp/opencode`, not written to
`output/`) took 4 min 35 s and selected epoch 250 by the cold-valid author loss
(1,407.4; the valid loss rises sharply afterwards — 2,338 at 350 and 24,829 at
500 — which is why the 50-epoch grid matters). Its fold-0 test AUPRC k=200 was
0.0886 (PCC 0.0756, Spearman 0.2213, R² −7.70, RMSE 0.3247). One fold, one
seed; not a ranking and not the ticket-16 run.

## How the adapter works

biolord is transductive: a cell's baseline is a per-row latent code learned at
registration (`_module.py:322-329`), not an encoder that can embed a new cell.
The adapter therefore:

1. **Registers one AnnData at fit**: all 2,523 training-role treated wells plus
   all 192 shared DMSO cells (`controls.cells`, 16 plates × 12 replicates).
   `X` is the frozen `log1p_cp10k`; `obsm["rdkit2d_dose"]` is the retained
   standardised RDKit2D descriptors (186 on fold 0; the column count is
   re-fitted per fold) + raw dose in µM; `obs["cell_type"]` is the
   constant `ipsc_cm`. The 192 control rows carry the DMSO vehicle attribute
   (`CS(C)=O` descriptor, dose 0) while training — at prediction it is
   overwritten. The row order and the `(control_group, replicate) → row` map are
   saved in `checkpoint/meta.json`.
2. **Fits** the single joint objective from random init. All registered rows
   are training data (`train_size=1.0`); Lightning validation is disabled
   (`limit_val_batches=0`). A callback evaluates the **cold valid drugs**
   every `checkpoint_grid` epochs: it predicts from each valid well's matched
   DMSO control row with the valid drug's attribute and computes the authors'
   loss (Gaussian NLL + τ·MSE + λ·‖z_unknown‖²). The best grid point is kept;
   `last.pt` is refreshed every grid so a timed-out fit can still predict. No
   AUPRC/PCC/test signal is used for selection.
3. **Predicts** by loading `best.pt` and `train_adata.h5ad`, then for each case
   row: index that well's matched control cell (`Controls.cell_for`), overwrite
   `obsm["rdkit2d_dose"]` for those registered rows with the test drug's vector,
   call `model.get_dataset(adata, indices=rows)` +
   `model.module.get_expression(...)` and take the decoder mean. This is the
   same construction `_compute_prediction` uses internally
   (`_model.py:848-858`); the categorical-only counterfactual helper is not
   used. Δ is the absolute prediction minus the **matched plate DMSO mean**
   from `controls.npz`, the benchmark Δ contract.

The RDKit2D descriptor set is descriptastorus `MakeGenerator(("RDKit2D",))`
with element 0 dropped (200 columns; 186 kept after dropping near-constant
columns) — see the attribute-encoding row in `FIDELITY.md`. Descriptor
standardisation is fitted on the unique **training** drugs only; the dose
column is appended raw, as the authors do.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 270 | 239.3 | 9.3 | 7.1 | 7446 | 500 | 250 |
| 1 | completed | 265 | 236.0 | 9.1 | 6.8 | 7160 | 500 | 500 |
| 2 | completed | 289 | 261.7 | 7.5 | 6.8 | 7091 | 500 | 500 |
| 3 | completed | 265 | 234.7 | 9.6 | 6.9 | 7113 | 500 | 500 |
| 4 | completed | 289 | 258.8 | 10.1 | 7.1 | 7094 | 500 | 500 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
