# biolord adapter — 5k sensitivity run card

> Derived from `stage2/03_models_cmbench/03.1_cold_drug_18k/models/biolord/`; the gene space is re-declared for the 5,000-HVG sensitivity axis (ADR 0028). Architecture, optimizer, schedule, seed and budget are unchanged.

**Status:** 5k sensitivity run complete (2026-10-22/23): AUPRC k=200 **0.1467** (95% CI 0.1398–0.1539), rank 6/12, `ideal_fit: partial`. Paired gap vs `delta_mean`: **-0.1025** (paired CI -0.1214..-0.0838); essentially unchanged from 18k (-0.1014). This is a sensitivity result, not a replacement for the 18k leaderboard.

- **Source:** `pert_model/biolord/source/` (vendored, commit `b7688790`, BSD-3-Clause;
  no published weights). No file under `pert_model/` is modified or generated
  (`sys.dont_write_bytecode` is set before the model import).
- **Environment:** `stage2/_envs/biolord` (venv on the `arc`
  interpreter, Python 3.13.13, `--system-site-packages`), overlay pins
  `scvi-tools==1.1.6 anndata==0.10.9 zarr==2.18.7 biolord==0.0.3` from
  `stage2/01_setup/execution/biolord/requirements.txt`, plus `descriptastorus==2.8.0`,
  `pandas_flavor==0.7.0` and `xarray==2024.3.0` for the author's RDKit2D
  features. torch is the `arc` wheel `2.12.0+cu130`; `torch.cuda.is_available()`
  is `True`. The model source is imported from
  `pert_model/biolord/source/src` (the PyPI `biolord==0.0.3` install supplies only the
  distribution metadata the version lookup expects).
- **Config:** [`config.json`](config.json) records the SI Table 4/5 author
  recipe (`rdkit2d_dose` ordered + `cell_type` categorical, τ=1e4, σ=20,
  decoder 4096×4, `n_latent` 256, batch 512, 500 epochs, Adam latent/decoder
  1e-4, attribute NN 2048×2 @ 1e-2, cosine to 1e-5) and the literal
  5,000-gene output.
- **Code:** `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/adapters/biolord/run.py`; shared plumbing in
  `adapters/common.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, two batches, fold 0) — holds the shared GPU lock
cd /home/chang/arc
BENCH_SMOKE=1 PYTHONDONTWRITEBYTECODE=1 flock "${TMPDIR:-/tmp}/stage2_bench_gpu.lock" \
  bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model biolord --fold 0 --timeout 1800
```

Without `BENCH_SMOKE=1` the adapter runs the frozen author budget (500 epochs,
checkpoint grid every 50 epochs). Artifacts land in
`stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/units/biolord/fold_0/`; the fold summary is appended to
`results/leaderboard.csv`.

## Smoke result (fold 0, shortened budget)

| quantity | value |
|---|---|
| fit wall | 16.7 s |
| total wall | 38 s |
| peak RSS | 3.5 GB |
| test AUPRC k=200 | 0.1240 |
| output | `(837, 5000)` finite |

The smoke proves wiring and output shape, not ranking; it is not part of the five-fold result.

<!-- ticket16-run -->
## Full five-fold run (5k)

| fold | status | wall s | fit s | pred valid s | pred test s | peak RSS GB | final/best step |
|---|---|---|---|---|---|---|---|
| 0 | completed | 136 | 114 | 8.5 | 9.9 | 3.5 | 500 / 300 |
| 1 | completed | 131 | 111 | 8.1 | 7.5 | 3.6 | 500 / 500 |
| 2 | completed | 133 | 110 | 7.4 | 11.5 | 3.6 | 500 / 500 |
| 3 | completed | 134 | 111 | 8.3 | 7.3 | 3.6 | 500 / 500 |
| 4 | completed | 131 | 111 | 8.4 | 7.5 | 3.6 | 500 / 500 |

Machine-readable copy: `results/run_cards.csv`.
