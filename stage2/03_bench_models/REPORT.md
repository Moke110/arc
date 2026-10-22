# Stage 2 whole-transcriptome benchmark — report

**Status:** executed 2026-10-21/22. 60/60 units completed (5 CPU baselines +
7 real models, five cold-drug folds each). Contract: ADR 0027. Artifacts:
`stage2/03_bench_models/{output,results}/`.

## Executive summary

- **No entrant beats the constant training-mean baseline on the pre-registered
  primary metric.** `delta_mean` scores 0.1989 macro AUPRC at k=200;
  `delta_mean_dose` 0.1948. The best entrant is **PerturbNet at 0.1861**, then
  PRnet at 0.1643. PerturbNet's paired difference vs `delta_mean` is
  **−0.0128 (95% CI −0.0248 to −0.0007)**; the benchmark therefore does not
  support selecting any entrant over the mean baseline.
- **The entrant ranking is still informative.** PerturbNet is the best entrant
  on the primary metric and on PCC/Spearman; PRnet is second and has the best
  Spearman. CPA is the best entrant on R²/RMSE but only 7th on AUPRC — global
  correlation and top/bottom direction ranking disagree.
- **Squidiff failed to converge at 18,164 genes**: cold-valid ε-MSE stayed at
  0.976 (random = 1.0) for the whole 10^5-step budget, so its DDIM samples are
  out of distribution (test RMSE 80.45). Its row is retained as a documented
  failure, not a meaningful ranking.
- **The result is robust to the environmental-chemical subset.** On the
  descriptive `atc-moa` subset the top four are unchanged: `delta_mean` 0.1662,
  `delta_mean_dose` 0.1632, PerturbNet 0.1432, PRnet 0.1337.

## Dataset and split

- View: `stage2/02_cm_data/datasets-1/temp_env_all-smiles.h5ad` — 4,478 wells
  (4,206 perturbed / 272 control), 460 drug keys, 18,164 genes, 16 physical
  plates, 3 doses (0.2/1/10 µM), 24 h, `ipsc_cm`.
- Δ = `log1p_cp10k(treated) − mean(log1p_cp10k(same-plate DMSO))`; group =
  `(drug_key, dose)`; 1,380 groups (3–6 wells each).
- Five cold-drug folds, 6:2:2: 276 train / 92 valid / 92 test drugs; every drug
  is test exactly once and valid exactly once; `valid(k) = test(k−1)`.
- Test bundles are metadata-only; no test expression enters fit or checkpoint
  selection. Split seed `20261021`; training seed `20261021`.

## Entrants and fidelity labels

| model | cold-drug path | ideal_fit | key declared deviations |
|---|---|---|---|
| PRnet | rFCFP4 × log-dose | full | 18k vs 978/5k genes; 2.5k wells; random init; cp10k-log1p vs median-normalized counts; matched-plate DMSO mean; Gaussian mean |
| CPA | frozen Morgan-2048 + unseen-drug dose scaler | full | 18k vs 5k genes; cold-drug encoder extension; full-fold AnnData setup (training rows only); 2,000-epoch budget |
| chemCPA | RDKit2D molecule encoder | partial | 978 L1000 pretrain omitted; 18k from scratch; 201 epochs |
| biolord | ordered RDKit2D + dose, transductive | partial | 18k vs 7.5k genes; registered control-row prediction; 50-epoch grid; descriptastorus descriptors |
| PerturbNet | ChemicalVAE → cINN | partial | 18k vs 5k/978 genes; ChemicalVAE ZINC → train-fold SMILES; no control transcriptome; 3-bin dose; deterministic ZINB rate; faithful py3.7/torch 1.13/scVI 0.7.1 |
| Squidiff | rFCFP4 × log-dose | full (form) | 200 → 18,164 genes (90×); no HVG; 10^5 fixed steps; matched DMSO cell per well; **did not converge** |
| CMonge | RDKit2D + dose → Monge map | full | 1,000 HVG → 18,164; 500k → 10k iterations; 50-d OT bottleneck; CPU JAX |

The complete per-model ledgers are `models/<model>/FIDELITY.md`; the
machine-readable summary is `results/fidelity.csv` (104 declared deviations:
48 `config_only`, 25 `recipe_substitution`, 18 `adapter_extension`,
12 `recipe_omission`, 1 `reproducibility_break`).

## Leaderboard (five cold-drug folds, k=200)

| rank | model | AUPRC@200 | 95% CI | k=100 | k=500 | PCC | Spearman | R² | RMSE | ideal_fit |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---|
| 1 | `delta_mean` | **0.1989** | 0.1826–0.2162 | 0.1847 | 0.2202 | 0.2291 | 0.1899 | 0.0198 | 0.1406 | none |
| 2 | `delta_mean_dose` | 0.1948 | 0.1789–0.2116 | 0.1806 | 0.2160 | 0.2272 | 0.1860 | 0.0171 | 0.1407 | none |
| 3 | `perturbnet` | 0.1861 | 0.1699–0.2031 | 0.1576 | 0.2290 | 0.2816 | 0.3288 | −1.3878 | 0.1709 | partial |
| 4 | `prnet` | 0.1643 | 0.1545–0.1750 | 0.1333 | 0.2171 | 0.2681 | **0.3840** | −0.7099 | 0.1714 | full |
| 5 | `biolord` | 0.0974 | 0.0914–0.1040 | 0.0750 | 0.1368 | 0.1077 | 0.2371 | −5.4742 | 0.2762 | partial |
| 6 | `ridge_fp` | 0.0921 | 0.0828–0.1024 | 0.0807 | 0.1154 | 0.0761 | 0.0584 | −0.7552 | 0.1702 | none |
| 7 | `cpa` | 0.0910 | 0.0863–0.0960 | 0.0683 | 0.1418 | **0.3151** | 0.3457 | **0.0376** | **0.1395** | full |
| 8 | `nearest_drug` | 0.0807 | 0.0717–0.0902 | 0.0694 | 0.1039 | 0.0635 | 0.0502 | −1.5777 | 0.1988 | none |
| 9 | `cmonge` | 0.0732 | 0.0698–0.0767 | 0.0547 | 0.1034 | 0.0366 | 0.1121 | −12.6827 | 0.4134 | full |
| 10 | `chemcpa` | 0.0608 | 0.0578–0.0640 | 0.0453 | 0.0994 | 0.2704 | 0.2866 | −0.0863 | 0.1462 | partial |
| 11 | `squidiff` | 0.0212 | 0.0207–0.0218 | 0.0143 | 0.0394 | 0.0003 | 0.0023 | −491,660.9 | 80.4537 | full |
| 12 | `no_change` | 0.0110 | 0.0110–0.0110 | 0.0055 | 0.0275 | — | — | −0.0195 | 0.1445 | none |

AUPRC is direction-macro, then group-macro, then fold-mean. `no_change`
reproduces the 1.1% prevalence exactly; its NaN PCC/Spearman is the contract's
constant-prediction policy. `squidiff`'s row is a convergence failure (below).

### Fold-level AUPRC@200

| model | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 |
|---|---:|---:|---:|---:|---:|
| `delta_mean` | 0.2145 | 0.1789 | 0.1947 | 0.2089 | 0.1974 |
| `delta_mean_dose` | 0.2091 | 0.1757 | 0.1909 | 0.2057 | 0.1925 |
| `perturbnet` | 0.1946 | 0.1766 | 0.1902 | 0.1977 | 0.1712 |
| `prnet` | 0.1552 | 0.1571 | 0.1753 | 0.1671 | 0.1669 |
| `biolord` | 0.0886 | 0.0964 | 0.1067 | 0.0951 | 0.1004 |
| `cpa` | 0.0914 | 0.0857 | 0.0867 | 0.1066 | 0.0845 |
| `cmonge` | 0.0725 | 0.0740 | 0.0809 | 0.0681 | 0.0705 |
| `chemcpa` | 0.0597 | 0.0632 | 0.0550 | 0.0654 | 0.0608 |
| `squidiff` | 0.0214 | 0.0211 | 0.0224 | 0.0204 | 0.0210 |

PerturbNet is the only entrant whose per-fold values approach the mean baseline
on every fold; PRnet is consistently ~0.03–0.06 below. No entrant has a fold
above the corresponding `delta_mean` fold.

### Paired bootstrap vs `delta_mean`

2,000 stratified 92-drug resamples, seed `20261022`; negative = worse than the
mean baseline.

| model | mean Δ AUPRC | 95% CI | P(entrant > mean) |
|---|---:|---|---:|
| `perturbnet` | −0.0128 | −0.0248 to −0.0007 | 0.021 |
| `prnet` | −0.0346 | −0.0474 to −0.0222 | 0.000 |
| `cpa` | −0.1081 | −0.1266 to −0.0902 | 0.000 |
| `biolord` | −0.1017 | −0.1183 to −0.0852 | 0.000 |
| `cmonge` | −0.1259 | −0.1434 to −0.1090 | 0.000 |
| `chemcpa` | −0.1382 | −0.1560 to −0.1206 | 0.000 |
| `squidiff` | −0.1778 | −0.1951 to −0.1613 | 0.000 |

PerturbNet is the only entrant whose CI excludes zero by a small margin; it is
slightly but consistently below the mean baseline.

### Per-drug view (460 test drug evaluations, five folds)

| model | drugs better than `delta_mean` | drugs worse | median Δ AUPRC | mean Δ AUPRC |
|---|---:|---:|---:|---:|
| `perturbnet` | **289** | 171 | **+0.0200** | −0.0128 |
| `prnet` | 234 | 226 | +0.0040 | −0.0345 |
| `cpa` | 183 | 277 | −0.0375 | −0.1079 |
| `biolord` | 152 | 308 | −0.0693 | −0.1014 |
| `cmonge` | 141 | 319 | −0.0866 | −0.1257 |
| `chemcpa` | 144 | 316 | −0.0741 | −0.1380 |
| `squidiff` | 87 | 373 | −0.1249 | −0.1776 |

PerturbNet wins on the **median** drug (289/460) but loses on the mean because a
minority of drugs are catastrophic. The pre-registered metric is the mean, so
the leaderboard keeps the mean; the median is a descriptive warning that a
robust/trimmed metric could tell a different story. This is not a second
leaderboard.

## Squidiff convergence failure

Squidiff is a genuine cold-drug model, but the author recipe did not transfer to
18,164 genes:

- The author's own SciPlex run used **200 HVGs**; this benchmark forbids HVG
  selection, so the diffusion MLP sees 18,164 outputs.
- The 10^5-step / lr 1e-4 / batch-64 budget is the author default. Cold-valid
  ε-MSE was 0.986 at step 10k and 0.977 at step 100k (best grid 0.9766;
  random = 1.0); training loss moved from 1.12 to 0.97.
- A checkpoint diagnostic on a training batch gives an ε-output std of 0.165
  against a target noise std of 1.0 and a correlation of 0.17. The model
  learned only a small fraction of the noise scale.
- DDIM then turns the under-predicted ε into `x_t/√ᾱ` at high t, producing
  samples with mean 61.7 and std 90.0 (the data are ~0–10), hence test RMSE
  80.45 and R² ≈ −4.9×10^5.

This is the 90× gene-space extrapolation flagged in ticket 09, not a scoring
bug. The run is kept in the leaderboard as a documented failure; it should not
be read as "Squidiff is a weak model" without that caveat.

## Execution record

- **60/60 units completed**; no `failed`, no `timed_out`, no `skipped_budget` in
  the final run. Every test group received one finite 18,164-vector; every
  model has five folds and is `winner_eligible` in the sense of fold coverage.
- **Wall-clock (five folds):** PRnet 169 s; chemCPA 383 s; PerturbNet 2,690 s
  (44.8 min); biolord 1,378 s (23.0 min); CMonge 10,065 s (2.80 h); CPA 11,591 s
  (3.22 h); Squidiff 16,965 s (4.71 h). All well inside the 16 h/model cap.
- **Peak total device memory** (including the ~2.5 GB desktop baseline):
  PerturbNet 9.46 GB; Squidiff 7.63 GB; biolord 7.45 GB; CPA 3.96 GB;
  PRnet 3.88 GB; chemCPA 3.85 GB; CMonge 2.74 GB. All below the 12.5 GB guard.
- **Reproducibility:** PRnet fold 0 was rerun from scratch into a temp
  directory; `pred_test.npy` was **bit-identical** (max abs diff 0.0) and the
  metrics JSON was identical (AUPRC 0.1552056269, PCC 0.2577438544).
- **Operational incidents (resolved, no effect on completed units):**
  1. The first queue exported
     `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`; torch 2.0.1
     (PRnet/CPA) and 1.13.1 (PerturbNet) reject it at CUDA init. The flag was
     removed and the affected units reran; failed dirs are archived at
     `output/_failed_alloc/`.
  2. The second queue was launched through an outer `tee` pipeline and was
     killed by SIGHUP during CMonge fold 0. It was relaunched with
     `setsid nohup` and file redirection; the partial CMonge dir is archived at
     `output/_failed_pipe/`. The final run is one uninterrupted queue.
- The `peak_rss_mb` column in `run_status.csv` is the cumulative child max from
  `RUSAGE_CHILDREN`, not a per-unit RSS; use `peak_gpu_mb` for the memory record.

## Baseline comparison

The two mean baselines are strong: predicting the training group-mean Δ (or the
per-dose mean) beats every entrant on AUPRC. This is not a bug in the metric:
the true top/bottom Δ sets share a large common component across drugs, and the
mean baseline captures that component. The entrant models add per-drug
structure, which helps correlation (PCC/Spearman) but not the direction-aware
top/bottom ranking under this protocol.

- `ridge_fp` (0.0921) and `nearest_drug` (0.0807) are weaker than the mean
  baselines, as expected; they are included to show that the leaderboard is not
  trivially solved by a simple structure baseline.
- `no_change` reproduces the 1.1% random prevalence; it is the sanity anchor,
  not a competitor.

## Exclusions (documented, not rows)

- **scVIDR**: no drug input; prediction is latent arithmetic between observed
  control and treated cells of the *same* drug, so an unseen drug has no rows.
- **DECANT**: no drug input in the forward pass, and the prediction target is
  the mean of the same treated bag it is handed; with 3 wells per group it sees
  its own answer. A cold-drug version is a new architecture, out of scope.

Both are excluded from Stage 2 by ticket 08. There is no warm-drug table.

## What this shows — and does not show

**Shows**
- Under one fixed cold-drug protocol at 18,164 genes, the seven capable models
  can be compared on the same groups, splits, and metric.
- PerturbNet is the best entrant; PRnet is second; CPA is the best regression
  model but not the best ranker; chemCPA is weak without its 978 pretrain;
  CMonge and biolord are near the simple baselines; Squidiff fails to converge.
- The mean baseline is the benchmark to beat; no entrant clears it.

**Does not show**
- That any model is ready for MIRI drug scoring or has clinical utility.
- That these numbers reproduce any author's reported benchmark; the gene space,
  corpus, splits, and budgets differ by design (ADR 0027).
- That 18,164 genes is the right output space for every model; the per-model
  ledgers record the extrapolations.
- That a different training budget, HVG policy, or HPO would not change the
  ranking. The contract deliberately forbids per-model tuning; a second effort
  can test those axes.
- A warm-drug result, a cold-dose interpolation result, or a deployment
  decision.

## Recommendation

1. **Do not name a model winner from this run.** The pre-registered primary
   metric puts `delta_mean` first. PerturbNet is the best entrant and the only
   one whose per-drug median is positive, but its mean is significantly below
   the mean baseline and its advantage is not robust.
2. **Carry the mean baseline forward as the benchmark to beat.** Any future MIRI
   adaptation must show an improvement over `delta_mean` on this protocol, not
   just a better rank than the other entrants.
3. **If a model must be chosen for the next adaptation,** PerturbNet is the
   best-supported entrant: highest entrant AUPRC (0.1861), strong PCC (0.282)
   and Spearman (0.329), a faithful py3.7/torch 1.13/scVI 0.7.1 recipe, and the
   lowest five-fold wall-clock among the large models (44.8 min). It is
   `partial` (no control transcriptome) and should be carried with that caveat.
   Its outlier drugs need a robustness check before it is used.
4. **Fix or bound the metric question before the next run.** The mean-vs-median
   disagreement and the CPA AUPRC-vs-regression disagreement are real signal:
   decide whether the MIRI-facing objective is direction-aware top/bottom
   ranking (current AUPRC), global correlation (PCC/Spearman), or reconstruction
   (R²/RMSE). Do not choose the metric after seeing the winner.

## Reproduce

```bash
# data contract + harness tests
cd stage2/03_bench_models && make test

# full queue (baselines + seven entrants, five folds)
bash stage2/03_bench_models/run_queue.sh

# independent acceptance checks
/home/chang/miniconda3/envs/arc/bin/python stage2/03_bench_models/validate_run.py

# regenerate the run cards and leaderboard from the raw fold predictions
/home/chang/miniconda3/envs/arc/bin/python stage2/03_bench_models/summarize_run.py
/home/chang/miniconda3/envs/arc/bin/python stage2/03_bench_models/aggregate_fidelity.py
/home/chang/miniconda3/envs/arc/bin/python stage2/03_bench_models/evaluate.py \
  --output-dir stage2/03_bench_models/output \
  --results-dir stage2/03_bench_models/results
```

Raw artifacts: `output/<model>/fold_<k>/` (predictions, checkpoints, status);
tables: `results/leaderboard.csv`, `results/metrics.csv`, `results/per_drug.csv`,
`results/bootstrap.npz`, `results/run_cards.csv`.
