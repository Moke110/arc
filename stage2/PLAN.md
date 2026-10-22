# Stage 2 — model setup and selection

## Goal

Select the open-source whole-transcriptome drug-perturbation model that should
produce the perturbomes feeding future MIRI drug scoring. Setup is deliberately
separate from benchmark design: it answers “does the model run in its own
environment?”, not “which model is best?”.

## Setup stage (2026-10-19)

Nine non-preprint, gate-pass models are exercised. Preprint models (X-Pert, DORA,
BaiZe, PerturbLDM) are excluded; WAVE and ExPO are non-preprint but fail the
revised gates (no dose / 978-only) and remain reference-only.

| model | native asset | setup target | caveat |
|---|---|---|---|
| PRnet | `checkpoint/lincs_best_epoch_all.pt` + `dataset/demo.h5ad` | `test_demo.py` on the 978-gene LINCS demo | **setup-only**: 978 is quarantined from Stage 2 evaluation; the 5,000-gene `sciplex_best_epoch_all.pt` ships too, but its processed `Sci_Plex.h5ad` is absent |
| CPA | `Pretrained_Models/Kang/model.pt` + `kang_normalized_hvg.h5ad` | load checkpoint, predict held-out Kang cells | Kang is IFN-β stimulation, not a chemical perturbation |
| chemCPA | `chemCPA_models/*.pt` | load and forward on an L1000-shaped batch | public weights are **L1000 977 / 2,000 genes, 188 drugs** |
| biolord | none published | train a small model on `GSM_new.h5ad` | per-drug categories, not a chemical-structure model |
| PerturbNet | LINCS 978 + sci-Plex 5,087 cINN checkpoints | load cINN + ChemicalVAE and run a latent forward | full inference needs Python 3.7 / TF 1.15 |
| scVIDR | none published; synthetic AnnData | 2,000-gene train + continuous-dose prediction | dose is log-linear latent interpolation |
| Squidiff | none published | 2,000-gene drug-structure forward | dose is `rFCFP × log10(dose+1)` |
| CMonge | bundled sciplex demo | official demo train + evaluate | 100-gene demo; categorical dose conditions |
| DECANT | none published | 2,000-gene synthetic forward with two doses | no processed cell bags or checkpoint |

**Definition of done (per model):**

1. a dedicated environment exists (conda, or the approved venv for biolord /
   CMonge);
2. the model loads or trains without editing the upstream clone;
3. one smoke run produces a prediction/forward artifact and a run log;
4. a run card records the exact command, weights, data, environment, and caveat.

Setup does **not** compare models, compute AUPRC, define folds, or adapt
anything to CM data. Those are selection-stage decisions.

## Environments

One environment per model. Seven models use the conda prefix
`arc-stage2-bringup-<model>`; biolord and CMonge use dedicated venv overlays
(ADR 0021 update). The shared `arc` environment is not modified.

| model | environment |
|---|---|
| PRnet | `arc-stage2-bringup-prnet` |
| CPA | `arc-stage2-bringup-cpa` |
| chemCPA | `arc-stage2-bringup-chemcpa` |
| biolord | `stage2/01_setup/biolord/.venv` |
| PerturbNet | `arc-stage2-bringup-perturbnet` |
| scVIDR | `arc-stage2-bringup-scvidr` |
| Squidiff | `arc-stage2-bringup-squidiff` |
| CMonge | `pert_model/cmonge/.venv` |
| DECANT | `arc-stage2-bringup-decant` |

## 978 quarantine

PRnet's native demo is the 978-landmark LINCS matrix. It is allowed **only** in
`stage2/01_setup/prnet/` and its outputs are labelled `setup-only`. No later
Stage 2 evaluation reads that artifact; the whole-transcriptome guard still
applies to Stage 2 selection.

## Layout

```
stage2/
├── PLAN.md                 # this file
├── 01_setup/               # one run card + run.sh per model
│   ├── prnet/
│   ├── cpa/
│   ├── chemcpa/
│   ├── biolord/
│   ├── perturbnet/
│   ├── scvidr/
│   ├── squidiff/
│   ├── cmonge/
│   └── decant/
├── 02_cm_data/             # frozen CM corpus + datasets-1 (ADR 0026)
├── 03_bench_models/        # executed cold-drug benchmark (ADR 0027)
├── audit/                  # the 14-model selection audit
└── _archive/cm_bench_v1/   # archived old benchmark; not a dependency
```

Stage 1 lives at the repository root in `stage1/` and is not part of this plan.

## Exit criteria for setup

- 9/9 run cards with a successful command and a saved artifact, or a documented
  blocked status with the exact error.
- No edits to `pert_model/<model>/` beyond documented data/weights and
  `VENDORED.txt`; adapters live in the run-card directory.
- `01_setup/SUMMARY.md` states which models are usable and which are blocked.

## Selection stage — frozen benchmark contract (2026-10-21)

The benchmark contract is fixed in
[`docs/adr/0027-stage2-benchmark-contract.md`](../docs/adr/0027-stage2-benchmark-contract.md)
and built in `stage2/03_bench_models/`. In brief:

- **Data:** `02_cm_data/datasets-1/temp_env_all-smiles.h5ad` — 4,478 wells,
  460 drug keys, 18,164 genes, 16 plates, three doses, 24 h, `ipsc_cm`;
  Δ = `log1p_cp10k(treated) − same-plate DMSO mean`; group = `(drug_key, dose)`.
- **Split:** five cold-drug folds, 6:2:2 (276 train / 92 valid / 92 test drugs),
  connectivity-block clusters, seed `20261021`; every drug is test exactly once.
- **Entrants:** PRnet, CPA, chemCPA, biolord, PerturbNet, Squidiff, CMonge.
  scVIDR and DECANT exit Stage 2 (no cold-drug structure path). One leaderboard;
  `ideal_fit` is a `full`/`partial` column.
- **Training:** random init, author defaults, no HPO; validation selects the
  epoch/checkpoint by the author loss; 18,164-gene output for every model;
  15.5 h soft / 16 h hard per model, 6 h/fold, one GPU, seed `20261021`.
- **Metric:** direction-aware macro AUPRC at k=200 (top-200 up + bottom-200
  down), group macro then fold mean; PCC/Spearman/R²/RMSE secondary; 92-drug
  stratified bootstrap (2,000 reps, seed `20261022`).
- **Execution:** `bash stage2/03_bench_models/run_queue.sh` runs the five CPU
  baselines and the seven entrants over five folds; per-unit statuses, run cards,
  fidelity ledgers, and the generated leaderboard live under
  `stage2/03_bench_models/{output,results,models}/`.

The benchmark does not fix the eventual MIRI adaptation; it produces the
selection evidence for it.

**Outcome (2026-10-22).** All 60 units completed. The mean baseline wins the
pre-registered primary metric (`delta_mean` 0.1989); the best entrant is
PerturbNet (0.1861, `partial`), followed by PRnet (0.1643). No entrant beats the
mean baseline, and Squidiff did not converge at 18,164 genes. See
`03_bench_models/REPORT.md` for the leaderboard, per-drug/robustness view, and
the recommendation not to name a model winner from this run.

