# CPA setup

**Status:** native checkpoint smoke (Kang IFN-β, not a chemical perturbation)

- source: `pert_model/CPA/source/` (commit `fbd7c025`), BSD-3-Clause
- weights: `pert_model/CPA/source/weights/Pretrained_Models/Kang/` (`model.pt`,
  `CPA_info.json`, `history.csv`)
- data: `pert_model/CPA/source/data/kang_normalized_hvg.h5ad` (13,576 × 5,000; 8 cell
  types, ctrl/stimulated); `Norman2019_normalized_hvg.h5ad` is stored alongside
- env: `arc-stage2-bringup-cpa` (Python 3.9, torch 2.0.1, scvi-tools 0.20.3,
  pandas 1.5.3)
- command:

```bash
bash stage2/01_setup/execution/cpa/run.sh
```

Expected output: `output/run.log`. The 2026-10-19 run wrote `cpa_pred.npy`
(2,000 × 5,000) and `cpa_pred_obs.csv`; deleted in ticket 06 — a rerun recreates
them.

**Result (2026-10-19 re-run):** passed. The Kang checkpoint loaded and predicted
**2,000 cells × 5,000 genes** (250 cells from each of 8 cell types; CPU; 5
samples). Required pins: scvi-tools 0.20.3, anndata 0.9.2, pandas 1.5.3,
numpy 1.23.5, scipy 1.10.1, jax/jaxlib 0.4.23, chex 0.1.86, pytorch-lightning
1.9.5, lightning 2.2.5, ray 2.9.3, pyarrow 14.0.2, setuptools 69.5.1. The
source checkout is used via `PYTHONPATH`; `smoke.py` supplies a runtime
`cpa-tools` `0.0.0` version fallback so the package need not be installed. (The
tracked `pydist/` stub dist-info that used to do this was deleted in ticket 06.)

Caveat: Kang is cytokine stimulation, not a drug/dose screen. This smoke proves
the CPA load/predict path, not cold-drug capability.
