# Stage 2 setup summary — 2026-10-19

**9/9 run cards completed** (chemCPA is load-only). Nine deployed non-preprint
models were exercised in their own environments. No upstream model code was
edited; adapters live under `stage2/01_setup/<model>/`.

| model | env | result | artifact |
|---|---|---|---|
| PRnet | `arc-stage2-bringup-prnet` | ✅ 661 × 978 | `prnet/output/demo_split_y_pre_array.csv` |
| CPA | `arc-stage2-bringup-cpa` | ✅ 2,000 × 5,000 | `cpa/output/cpa_pred.npy` |
| chemCPA | `arc-stage2-bringup-chemcpa` | ⚠️ load-only 8 × 977 | `chemcpa/output/pred.npy` |
| biolord | `biolord/.venv` | ✅ 640 × 4,999 | `biolord/output/pred_mean.npy` |
| PerturbNet | `arc-stage2-bringup-perturbnet` | ✅ checkpoints + latent forward | `perturbnet/output/meta.json` |
| scVIDR | `arc-stage2-bringup-scvidr` | ✅ 2,000-gene train + 3 doses | `scvidr/output/pred_doses.npz` |
| Squidiff | `arc-stage2-bringup-squidiff` | ✅ 2,000-gene forward | `squidiff/output/pred.npy` |
| CMonge | `pert_model/cmonge/.venv` | ✅ demo R² 0.785 / 0.871 / 0.808 | `cmonge/output/meta.json` |
| DECANT | `arc-stage2-bringup-decant` | ✅ 2,000-gene forward | `decant/output/pred_gene_dose{1,100}.npy` |

## What "9/9" means: completion, not usability

The count measures whether each model's code path runs in its own environment,
not whether it can produce a whole-transcriptome CM perturbome today.

| model | what was actually exercised | data / weights | ready for selection? |
|---|---|---|---|
| PRnet | 978-gene demo forward | real LINCS 978 checkpoint | no: 5,000-gene input data missing, 978 quarantined |
| CPA | Kang checkpoint load + predict | real public checkpoint (IFN-β) | no: needs chemical/dose retraining |
| chemCPA | checkpoint load + synthetic forward | public L1000 977/2,000 weights | no: load-only, needs retraining |
| biolord | 3-epoch train + counterfactual | real GSM data, no public weights | no: per-dataset retraining, no cold-drug encoder |
| PerturbNet | cINN + ChemicalVAE latent forward | released 978/5,087 weights | partial: full inference needs TF 1.15 |
| scVIDR | synthetic train + dose prediction | synthetic data, no weights | no: architecture-only |
| Squidiff | real SMILES/rFCFP forward | random weights (none released) | no: architecture-only |
| CMonge | official 100-gene demo train + eval | bundled demo | no: demo-only |
| DECANT | 2,000-gene forward at two doses | random weights, no data | no: architecture-only |

**No model currently produces a whole-transcriptome prediction on CM data.**
The next stage must choose a gene space and a retraining strategy before any
model can be compared.

## What is ready

- **PRnet**: code and the 978 LINCS checkpoint run end-to-end. The 978 demo is
  **setup-only** (ADR 0024). The 5,000-gene Sci-Plex checkpoint
  (`sciplex_best_epoch_all.pt`) ships with the clone, but the
  `dataset/Sci_Plex.h5ad` that `test_sciplex.py` reads is not distributed, so
  only the 978 path was exercised. A Stage 2 evaluation needs that processed
  data or a retraining run.
- **CPA**: the public Kang checkpoint loads and predicts. Kang is IFN-β
  stimulation, not a chemical screen, so CPA still needs a chemical/dose training
  run before it can compete on the MIRI task.
- **biolord**: training and counterfactual prediction work on `GSM_new`. Its
  counterfactual API supports categorical targets only, and it learns per-drug
  categories rather than a chemical-structure encoder.
- **chemCPA**: architecture/checkpoint compatibility is proven, but the public
  weights are **L1000 977/2,000-gene** `ComPert` models (188 drugs, 3 covariates).
  There is no 5,000/20,570-gene public checkpoint; real-data inference on
  `combo_sciplex` is blocked without retraining.
- **PerturbNet**: the released LINCS/sci-Plex cINN + ChemicalVAE weights load into
  the in-repo architecture and run a latent forward. Full inference still needs
  the Python 3.7 / TF 1.15 cell-VAE.
- **scVIDR**: a 2,000-gene synthetic train and continuous-dose prediction work.
  The dose curve is a log-linear latent interpolation, not a learned response.
- **Squidiff**: a 2,000-gene drug-structure forward works. Dose is encoded as an
  rFCFP fingerprint scaled by `log10(dose+1)`.
- **CMonge**: the official demo trains and evaluates on the bundled sciplex data.
  Its dose is a categorical `drug-dose` condition; the demo is 100 genes.
- **DECANT**: the model runs at 2,000 genes and responds to `dose_log1p`. There is
  no processed data or checkpoint; it remains a representation model.

## Selection-stage implications

1. chemCPA's public weights cannot be used as a whole-transcriptome entrant
   without a Tier S retraining run; only its chemical encoder is reusable.
2. PRnet and CPA have usable public checkpoints only in 978/5,000-gene spaces
   (and PRnet's 5,000-gene path lacks its processed input data); the 20,570-gene
   evaluation still requires retraining or adaptation.
3. biolord has no public weights and no structural cold-drug encoder.
4. PerturbNet is the only new deployed model with released >5,000-gene weights
   already staged locally; scVIDR and Squidiff release weights on external
   archives, but their full inference environments/data are not staged here.
5. The next decision is therefore a **gene-space/retraining strategy**, not a
   leaderboard: decide which model(s) can be trained at the Stage 2 gene space
   within the 16 h/model budget.

See `../PLAN.md` for the deferred selection-stage decision list.

## Excluded from setup

- **Preprints**: X-Pert, DORA, BaiZe, PerturbLDM — excluded from Stage 2 by
  policy; none has a Stage 2 vendor clone or run card. The Stage 1 XPert library
  remains at `pert_model/xpert/` and is not a Stage 2 entrant.
- **Gate-fail references**: WAVE (no dose), ExPO (978-only).
