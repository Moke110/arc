# CMonge setup smoke

**Status:** official demo train + evaluate on bundled sciplex data

- source: `pert_model/CMonge/source/` (commit `195f5cd1`), MIT
- data: `pert_model/CMonge/source/data/dummy_data.h5ad` (4,000 × 100, three conditions)
- env: shared bench venv `stage2/_envs/cmonge` (Python 3.11, jax 0.4.38, ott-jax 0.5.2); the
  model-local `.venv` was a duplicate and was removed in ticket 04
- command:

```bash
bash stage2/01_setup/execution/cmonge/run.sh
```

Expected output: `output/meta.json`, `output/run.log`. `demo_logs.yml`, `ae/` and
`embed/` from the 2026-10-19 run were deleted in ticket 06; a rerun recreates them.

Caveat: the demo is 100 genes; CMonge's gene space is dataset-defined, but its
dose is a categorical `drug-dose` condition, not a continuous covariate.

**Result (2026-10-19):** passed. Official demo AE + 10 CMonge iterations; validation R² 0.785 / 0.871 / 0.808 for the three conditions.
