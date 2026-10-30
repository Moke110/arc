# Stage 2 setup

Nine deployed non-preprint models, one run card each. The goal is to prove each
model runs in its own pinned environment — not to compare them. Results below
are the 2026-10-19 bring-up record; every model was later exercised for real in
the benchmark step (`03_models_cmbench/`, re-laid-out in ticket 07).

- **HTML report for non-computational readers:**
  [`model_architecture_weights_benchmarks.html`](model_architecture_weights_benchmarks.html) —
  architecture (drug encoder vs. main predictor), weights, and published
  benchmark claims per model.

## Layout

```
01_setup/
├── AGENTS.md                                     routing index for agents
├── README.md                                     this file (the one human summary)
├── model_architecture_weights_benchmarks.html    human report, stays at step level
└── execution/<model>/
    ├── README.md       run card: source commit, weights, data, env, command, result
    ├── run.sh          entry point (repo root found by walking up to .git)
    ├── smoke.py        the actual smoke; prnet uses run_prnet.py instead
    ├── env.lock        recorded environment (conda `@EXPLICIT`/pip freeze, 2026-10)
    └── output/         run.log + meta.json (the textual evidence, tracked)
```

`output/` keeps only the textual evidence: `run.log` for all nine models and
`meta.json` for the six that emit one. The prediction arrays of the 2026-10-19
run (`*.npy`, `*.npz`, the four `demo_split_*_array.csv`, `pred_obs.csv`, the
`model/` checkpoint, CMonge's `ae/` store and `demo_logs.yml`) were deleted in
the ticket-06 restructure; a rerun regenerates them and `.gitignore` keeps them
out of git. The pre-trim output tree was 116.6 MB, the surviving evidence 78 KB.

## Run cards

| model | env | status | what was actually exercised | surviving evidence |
|---|---|---|---|---|
| PRnet | `arc-stage2-bringup-prnet` | ✅ 661 × 978 | 978-gene LINCS demo forward; **setup-only** | `execution/prnet/output/run.log` |
| CPA | `arc-stage2-bringup-cpa` | ✅ 2,000 × 5,000 | Kang checkpoint load + predict (IFN-β) | `execution/cpa/output/run.log` |
| chemCPA | `arc-stage2-bringup-chemcpa` | ⚠️ load-only | public L1000 977/2,000-gene checkpoint forward | `execution/chemcpa/output/{meta.json,run.log}` |
| biolord | local `.venv` | ✅ 640 × 4,999 | 3-epoch from-scratch train + counterfactuals | `execution/biolord/output/run.log` |
| PerturbNet | `arc-stage2-bringup-perturbnet` | ✅ checkpoints + latent | released LINCS/sci-Plex cINN + ChemicalVAE load | `execution/perturbnet/output/{meta.json,run.log}` |
| scVIDR | `arc-stage2-bringup-scvidr` | ✅ 2,000-gene train + 3 doses | synthetic train + continuous-dose predict | `execution/scvidr/output/{meta.json,run.log}` |
| Squidiff | `arc-stage2-bringup-squidiff` | ✅ 2,000-gene forward | real SMILES/rFCFP drug-structure forward | `execution/squidiff/output/{meta.json,run.log}` |
| CMonge | `stage2/_envs/cmonge` | ✅ demo R² 0.785 / 0.871 / 0.808 | official sciplex demo train + evaluate | `execution/cmonge/output/{meta.json,run.log}` |
| DECANT | `arc-stage2-bringup-decant` | ✅ 2,000-gene, two doses | synthetic forward with dose context | `execution/decant/output/{meta.json,run.log}` |

## Reproduce

Run from anywhere; each `run.sh` resolves the repository root by walking up to
the directory containing `.git` (ADR 0030 §4), so no `cd` and no
depth-dependent `../..` is involved. From the repository root:

```bash
# cheapest full smoke: CPU torch, no downloaded weights (~seconds)
bash stage2/01_setup/execution/decant/run.sh

# the other eight, all CPU
bash stage2/01_setup/execution/chemcpa/run.sh
bash stage2/01_setup/execution/cmonge/run.sh
bash stage2/01_setup/execution/cpa/run.sh
bash stage2/01_setup/execution/perturbnet/run.sh
bash stage2/01_setup/execution/prnet/run.sh
bash stage2/01_setup/execution/scvidr/run.sh
bash stage2/01_setup/execution/squidiff/run.sh
bash stage2/01_setup/execution/biolord/run.sh   # needs the local .venv, see below
```

Each card writes `<model>/output/run.log` and, where the model emits one,
`<model>/output/meta.json`, then asserts the evidence file exists. Environments:
six `arc-stage2-bringup-*` conda envs plus the shared `stage2/_envs/cmonge` venv
and the biolord `.venv` overlay. The biolord venv is a machine-local overlay
(Python 3.13 on the `arc` base plus 2024-era pins); it is not in git — recreate
it with

```bash
python -m venv --system-site-packages stage2/01_setup/execution/biolord/.venv
stage2/01_setup/execution/biolord/.venv/bin/pip install \
  -r stage2/01_setup/execution/biolord/requirements.txt
```

`env.lock` in each card is a 2026-10 record of the environment the run used; the
envs are already built, so the locks are documentation, not re-resolved.

## What "9/9" means: completion, not usability

The count measures whether each model's code path runs in its own environment,
not whether it can produce a whole-transcriptome CM perturbome today.

| model | data / weights | ready for selection? |
|---|---|---|
| PRnet | real LINCS 978 checkpoint | no: 5,000-gene input data missing, 978 quarantined |
| CPA | real public checkpoint (IFN-β) | no: needs chemical/dose retraining |
| chemCPA | public L1000 977/2,000 weights | no: load-only, needs retraining |
| biolord | real GSM data, no public weights | no: per-dataset retraining, no cold-drug encoder |
| PerturbNet | released 978/5,087 weights | partial: full inference needs TF 1.15 |
| scVIDR | synthetic data, no weights | no: architecture-only |
| Squidiff | random weights (none released) | no: architecture-only |
| CMonge | bundled demo | no: demo-only |
| DECANT | random weights, no data | no: architecture-only |

**No model produced a whole-transcriptome prediction on CM data at setup time.**
The stage therefore had to choose a gene space and a retraining strategy before
any model could be compared.

## What is ready

- **PRnet**: code and the 978 LINCS checkpoint run end-to-end. The 978 demo is
  **setup-only** (ADR 0024) and quarantined from evaluation. The 5,000-gene
  Sci-Plex checkpoint ships with the clone, but the `dataset/Sci_Plex.h5ad` that
  `test_sciplex.py` reads is not distributed, so only the 978 path was exercised.
- **CPA**: the public Kang checkpoint loads and predicts. Kang is IFN-β
  stimulation, not a chemical screen; a chemical/dose training run is still
  needed. `smoke.py` carries a runtime `cpa-tools` version fallback so the
  source checkout need not be installed.
- **biolord**: training and counterfactual prediction work on `GSM_new`; the
  counterfactual API is categorical-only, and it learns per-drug categories
  rather than a structure encoder.
- **chemCPA**: architecture/checkpoint compatibility is proven, but the public
  weights are **L1000 977/2,000-gene** `ComPert` models (188 drugs, 3
  covariates). Real-data inference on `combo_sciplex` is blocked without
  retraining.
- **PerturbNet**: the released LINCS/sci-Plex cINN + ChemicalVAE weights load
  into the in-repo architecture and run a latent forward; full inference still
  needs the Python 3.7 / TF 1.15 cell-VAE.
- **scVIDR**: 2,000-gene synthetic train and continuous-dose prediction work; the
  dose curve is a log-linear latent interpolation, not a learned response.
- **Squidiff**: a 2,000-gene drug-structure forward works; dose is an rFCFP
  fingerprint scaled by `log10(dose+1)`.
- **CMonge**: the official demo trains and evaluates on the bundled sciplex data;
  its dose is a categorical `drug-dose` condition and the demo is 100 genes.
- **DECANT**: the model runs at 2,000 genes and responds to `dose_log1p`; there
  is no processed data or checkpoint, so it remains a representation model.

## Selection-stage implications

1. chemCPA's public weights cannot be a whole-transcriptome entrant without a
   retraining run; only its chemical encoder is reusable.
2. PRnet and CPA have usable public checkpoints only in 978/5,000-gene spaces;
   the 20,570-gene evaluation requires retraining or adaptation.
3. biolord has no public weights and no structural cold-drug encoder.
4. PerturbNet was the only new deployed model with released >5,000-gene weights
   already staged locally; scVIDR and Squidiff release weights on external
   archives, but their full inference environments/data were not staged here.
5. The next decision was therefore a **gene-space/retraining strategy**, not a
   leaderboard.

**Key finding:** the public chemCPA weights are L1000 **977/2,000-gene**
`ComPert` state dicts (188 drugs, 3 covariates), not 5,000/20,570-gene models;
the run card therefore uses the load-only fallback.

**Update (2026-10-22):** superseded by the benchmark step, where the
cold-drug-capable models produce 18,164-gene CM predictions and are scored over
five cold-drug folds (60/60 units; 5 CPU baselines included). The setup
conclusion above remains the correct statement for the setup stage itself.

## Excluded from setup

- **Preprints**: X-Pert, DORA, BaiZe, PerturbLDM — no Stage 2 vendor clone or
  run card; `pert_model/XPert/` is the Stage 1 XPert library, not an entrant.
- **Gate-fail references**: WAVE (no dose), ExPO (978-only).

Common conventions: one environment per model (`arc` itself stays clean); no
upstream model code is edited — adapters live in `execution/<model>/`; outputs
stay in `execution/<model>/output/` and out of git except the textual evidence.
