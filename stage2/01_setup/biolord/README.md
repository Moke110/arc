# biolord setup

**Status:** from-scratch smoke (no public checkpoint exists)

- source: `pert_model/biolord/` (commit `b7688790`), BSD-3-Clause
- weights: none published
- data: `pert_model/biolord/data/GSM_new.h5ad` (14,811 × 4,999; 5 drugs × 7
  doses, A549)
- env: local `.venv` (Python 3.13 from `arc`, torch 2.14.1+cu130, scvi-tools
  1.1.6, anndata 0.10.9, zarr 2.18.7, biolord 0.0.3 metadata + repo source via
  `PYTHONPATH`)
- command:

```bash
bash stage2/01_setup/biolord/run.sh
```

Environment recreation (only needed if `.venv` is missing; the overlay reuses
the `arc` base for torch/scanpy and pins the 2024-era biolord dependencies):

```bash
cd /home/chang/arc
python -m venv --system-site-packages stage2/01_setup/biolord/.venv
stage2/01_setup/biolord/.venv/bin/pip install -r stage2/01_setup/biolord/requirements.txt
```

The smoke trains a small biolord model for a few epochs with `drug` as a
categorical attribute and `dose_val` as an ordered attribute, then predicts
counterfactuals from control cells. Expected output: `output/pred_mean.npy`,
`output/pred_obs.csv`, `output/model/`, and `output/run.log`.

**Result (2026-10-19 re-run):** passed. Trained on GPU for 3 epochs; produced
**640 counterfactual rows × 4,999 genes** (128 control cells × 5 drug
categories). biolord's counterfactual API only accepts categorical target
attributes, so `dose_val` was held at each source cell's value; the run card
records this. `anndata`/`zarr`/`scvi-tools` had to be pinned for the 2024-era
biolord code (see `requirements.txt`).

Caveat: biolord learns per-drug categories; it is not a chemical-structure
cold-drug model. The smoke proves the training/prediction path only.
