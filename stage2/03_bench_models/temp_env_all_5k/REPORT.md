# Stage 2 5,000-gene HVG sensitivity benchmark — report

**Status:** executed 2026-10-22/23. 60/60 units completed (5 CPU baselines +
7 real models, five cold-drug folds each). Decision: ADR 0028. Artifacts:
`stage2/03_bench_models/temp_env_all_5k/{output,results}/`.

This is a **sensitivity run**, not a replacement for the executed 18k
benchmark. It holds wells, Δ, controls, groups, folds, entrants, seeds and
budgets fixed and changes one thing: the gene axis, from all 18,164 corpus
genes to the 5,000 genes with the smallest sum of per-source HVG ranks.

## Executive summary

- **The 5k leaderboard changes the headline conclusion.** PRnet ranks first at
  **0.2735** macro AUPRC (k=200), ahead of `delta_mean` at **0.2491**. Its paired
  per-drug gap is **+0.0243 (95% CI +0.0113 to +0.0367)**, it wins in all five
  folds, and it also leads on the descriptive `atc-moa` subset (0.2261 vs
  0.2131). It is the only entrant above the mean baseline on this axis.
- **Raw AUPRC must not be compared across the two runs.** The 5k denominator
  raises random prevalence from 1.1% to **4.0%** at k=200, so every row rises;
  `no_change` moves from 0.0110 to 0.0400. The comparable quantity is the
  paired gap against `delta_mean` on the same folds.
- **The 18k best entrant does not survive as the best entrant.** PerturbNet falls
  from −0.0128 to **−0.0563** versus the baseline (change −0.0435, 95% CI
  −0.0495 to −0.0375), although its absolute score rises. CPA improves its gap
  slightly (+0.0147) but stays below; biolord, chemCPA, CMonge and Squidiff
  remain below.
- **Squidiff still fails to converge at 5,000 genes.** Cold-valid ε-MSE is
  0.889 at step 10k and then sits at ≈1.0 for the rest of the 10^5-step budget;
  test R² is −2.98e5 and RMSE 74.3. The gene space was not the cause of its
  failure.
- **Recommendation:** keep the 18k result as the primary benchmark and cite the
  5k run as evidence that the "no entrant beats the mean baseline" conclusion
  is **gene-space dependent**. PRnet's 5k advantage is consistent but modest
  (0.024 AUPRC; 62% of test drugs), and it does not overturn the 18k result on
  its own.

## The 5k gene axis

- Sources: `dtoxs_cm`, `screen_ctox`, `temp_env_all-smiles` (all 18,164 genes,
  identical `var` order).
- Per source: `scanpy.pp.highly_variable_genes(flavor="seurat",
  n_top_genes=None)` on `layers["log1p_cp10k"]` over all wells, no `batch_key`.
  Genes with mean in `[0.0125, 3]` are ranked 1..M by descending
  `dispersions_norm` (ties by corpus order); all other genes, including
  zero-mean/NaN, get rank 18,165. Sum the three ranks; keep the 5,000 smallest,
  ties by corpus order; output in corpus order.
- Eligible genes per source: 11,127 / 11,976 / 17,995. Low-expression genes in
  the final 5,000: 41 / 72 / 13. Rank-sum range of the selection: 64–21,875.
  Provenance: `data/gene_set_5k.manifest.json`; full ranks:
  `data/gene_set_5k_ranks.csv`.
- The selection is unsupervised and uses all wells, including wells of drugs
  that are test in some fold. This matches the authors' pipelines and keeps one
  axis across folds; it is declared in ADR 0028. No expression labels or Δ
  values enter the selection.
- Data parity: `build_data.py` verifies every rebuilt array against the 18k
  archive at the same gene columns, and verifies that `wells.parquet`,
  `groups.parquet` and `splits.parquet` are identical. The folds are therefore
  paired per drug and per fold with the 18k run.

## Dataset and split

- View: `temp_env_all-smiles` — 4,478 wells (4,206 perturbed / 272 control),
  460 drug keys, 16 plates, 3 doses (0.2/1/10 µM), 24 h, `ipsc_cm`.
- Δ = `log1p_cp10k(treated) − mean(log1p_cp10k(same-plate DMSO))`; group =
  `(drug_key, dose)`; 1,380 groups.
- Five cold-drug folds, 6:2:2: 276 train / 92 valid / 92 test drugs; every drug
  is test exactly once and valid exactly once. Test bundles are metadata-only.
  Split seed `20261021`; training seed `20261021`.

## Entrants and gene-space fidelity

| model | cold-drug path | 5k run setting | ideal_fit | gene-space note |
|---|---|---|---|---|
| PRnet | rFCFP4 × log-dose | 5,000 output | full | matches the widest author space (5,000 Sci-Plex); gene identity differs |
| CPA | frozen Morgan-2048 + unseen-drug dose scaler | 5,000 | full | width matches the authors' 5,000 HVG; different HVG statistic |
| chemCPA | RDKit2D molecule encoder | 5,000 from scratch | partial | 2.5× the widest author gene set; 978 pretrain still omitted |
| biolord | ordered RDKit2D + dose, transductive | 5,000 | partial | close to the 4,999 demo, 1.5× narrower than 7,561 |
| PerturbNet | ChemicalVAE → cINN | 5,000 | partial | matches sci-Plex (5,087), 5.1× LINCS |
| Squidiff | rFCFP4 × log-dose | 5,000 | full (form) | 25× the authors' 200 HVGs; still non-convergent |
| CMonge | RDKit2D + dose → Monge map | 5,000 | full | 5× the authors' 1,000 HVGs; 50-d OT bottleneck |
| 5 baselines | — | unchanged | none | train-mean, per-dose mean, ridge, nearest drug, zero |

Architecture, optimizer, schedule, seed and budget are unchanged from 18k; only
the output width and the gene-space ledger rows differ. All seven adapters
passed a fold-0 smoke at 5k before the queue (finite `(837, 5000)` outputs).

## Leaderboard (five cold-drug folds, k=200)

Headline = direction-macro, group-macro AUPRC at 200 of 5,000 genes.
`no_change` is the visible chance floor (0.0400), not a competitive entry.

| rank | model | AUPRC@200 | 95% CI | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 | ideal_fit |
|---:|---|---:|---|---:|---:|---:|---:|---:|---|
| 1 | PRnet | **0.2735** | 0.2590–0.2891 | 0.2872 | 0.2659 | 0.2664 | 0.2771 | 0.2707 | full |
| 2 | `delta_mean` | 0.2491 | 0.2321–0.2673 | 0.2662 | 0.2324 | 0.2417 | 0.2595 | 0.2459 | none |
| 3 | `delta_mean_dose` | 0.2448 | 0.2284–0.2625 | 0.2603 | 0.2290 | 0.2378 | 0.2560 | 0.2407 | none |
| 4 | PerturbNet | 0.1929 | 0.1790–0.2083 | 0.1964 | 0.1902 | 0.1863 | 0.2090 | 0.1824 | partial |
| 5 | CPA | 0.1559 | 0.1513–0.1605 | 0.1453 | 0.1642 | 0.1593 | 0.1681 | 0.1427 | full |
| 6 | biolord | 0.1467 | 0.1398–0.1539 | 0.1377 | 0.1473 | 0.1628 | 0.1364 | 0.1491 | partial |
| 7 | `ridge_fp` | 0.1387 | 0.1284–0.1501 | 0.1519 | 0.1354 | 0.1460 | 0.1295 | 0.1308 | none |
| 8 | `nearest_drug` | 0.1252 | 0.1155–0.1355 | 0.1288 | 0.1167 | 0.1211 | 0.1330 | 0.1264 | none |
| 9 | chemCPA | 0.1165 | 0.1129–0.1203 | 0.1141 | 0.1199 | 0.1226 | 0.1094 | 0.1164 | partial |
| 10 | CMonge | 0.1153 | 0.1108–0.1200 | 0.1128 | 0.1160 | 0.1238 | 0.1082 | 0.1155 | full |
| 11 | Squidiff | 0.0523 | 0.0517–0.0530 | 0.0527 | 0.0522 | 0.0537 | 0.0511 | 0.0519 | full (failure) |
| 12 | `no_change` | 0.0400 | 0.0400–0.0400 | 0.0400 | 0.0400 | 0.0400 | 0.0400 | 0.0400 | none |

Secondary metrics (headline, all folds):

| model | AUPRC@100 | AUPRC@500 | PCC | Spearman | R² | RMSE |
|---|---:|---:|---:|---:|---:|---:|
| PRnet | 0.2229 | 0.3571 | 0.3498 | 0.4051 | −0.2305 | 0.1731 |
| `delta_mean` | 0.2172 | 0.3101 | 0.2433 | 0.2009 | 0.0164 | 0.1657 |
| PerturbNet | 0.1555 | 0.2652 | 0.2261 | 0.2268 | −13.97 | 0.5312 |
| CPA | 0.1095 | 0.2535 | 0.2862 | 0.3289 | 0.0163 | 0.1669 |
| biolord | 0.1109 | 0.2176 | 0.0756 | 0.2019 | −7.01 | 0.3547 |
| Squidiff | 0.0315 | 0.1124 | 0.0019 | 0.0005 | −2.98e5 | 74.34 |

## Paired comparison against the 18k run

The two runs share wells, folds and test drugs, so the per-drug difference is
paired. `gap` = model AUPRC minus `delta_mean` AUPRC on the same drug, averaged
over 460 drugs. The last column is the change in that gap from 18k to 5k
(positive = the model improved relative to the baseline). CIs are a 2,000-draw
drug bootstrap, seed 20261022.

| model | 5k gap | 95% CI | 18k gap | change in gap | 95% CI |
|---|---:|---|---:|---:|---|
| PRnet | **+0.0243** | +0.0113 to +0.0367 | −0.0345 | **+0.0589** | +0.0516 to +0.0659 |
| PerturbNet | −0.0563 | −0.0688 to −0.0434 | −0.0128 | −0.0435 | −0.0495 to −0.0375 |
| CPA | −0.0932 | −0.1120 to −0.0746 | −0.1079 | +0.0147 | +0.0089 to +0.0204 |
| biolord | −0.1025 | −0.1214 to −0.0838 | −0.1014 | −0.0010 | −0.0050 to +0.0032 |
| chemCPA | −0.1327 | −0.1498 to −0.1152 | −0.1380 | +0.0054 | +0.0017 to +0.0092 |
| CMonge | −0.1339 | −0.1538 to −0.1150 | −0.1257 | −0.0082 | −0.0119 to −0.0046 |
| Squidiff | −0.1968 | −0.2150 to −0.1790 | −0.1776 | −0.0192 | −0.0216 to −0.0167 |

PRnet beats `delta_mean` on **62.2%** of the 460 test drugs (50.9% at 18k);
PerturbNet beats it on 45.2% (62.8% at 18k). The `atc-moa` subset agrees with
the full panel: PRnet 0.2261, `delta_mean` 0.2131, `delta_mean_dose` 0.2100,
PerturbNet 0.1610, CPA 0.1509.

## Squidiff convergence failure

Squidiff's cold-valid ε-MSE is 0.889 at the first 10k-step checkpoint and then
≈1.0 (random) through 100k steps on every fold. Its DDIM samples remain out of
distribution: test R² −2.98e5, RMSE 74.34, and AUPRC 0.0523, barely above the
0.0400 chance floor. Shrinking the decoder from 18,164 to 5,000 genes did not
fix the optimisation; the row stays a documented failure, not a ranking.

## Execution record

- 60/60 units completed; total measured wall time **12.7 h** (sum of unit wall
  clocks; one GPU unit at a time).
- CMonge is the outlier: 22,139 s total (6.15 h), CPU JAX, 1.1–1.6 GB RSS,
  versus 32 min/fold for the 18k run. CPA 10,545 s, Squidiff 9,409 s,
  PerturbNet 2,305 s. The other entrants and baselines are minutes.
- Peak RSS: PerturbNet 7.83 GB, biolord 3.72 GB, Squidiff 1.63 GB, CPA
  1.62 GB, CMonge 1.55 GB. The queue cap was raised from 8 GB to 10 GB after
  PerturbNet peaked near the old limit.
- CMonge folds 0 and 4 resumed from their own checkpoints after scheduler
  interruptions; resume is deterministic (same seed, same data) and the final
  statuses are `completed`.
- No test expression entered fit or checkpoint selection. Tests were run over
  raw fold predictions, not over adapter-reported scores.

## What this shows — and does not show

**Shows**

- The 18k conclusion is gene-space dependent. On a 5,000-gene HVG axis, PRnet
  beats the training-mean baseline by 0.024 AUPRC with a paired CI above zero,
  in every fold and on the pharmacology subset.
- The entrant ranking is not stable: PerturbNet moves from best entrant to
  4th overall, behind PRnet and the two mean baselines (still the second-best
  entrant); CPA's gap improves while PerturbNet's widens.
- Squidiff's failure is not caused by the 18k output width.

**Does not show**

- That PRnet is the better model for MIRI. The 5k run is one sensitivity axis;
  the 18k result remains primary, and PRnet's 5k margin is modest and
  gene-space-specific.
- That the 5k axis is the authors' preprocessing. It is one fixed Seurat HVG
  selection shared across three views, not a per-dataset or per-fold HVG fit.
- That the 5k numbers can be compared in absolute terms with 18k. Prevalence
  differs 3.6×; only paired gaps and ranks are meaningful.
- That the drug-level bootstrap is a population-level test: the 460 drugs are
  structurally clustered, so the CIs are optimistic relative to a
  cluster-aware bootstrap.

## Recommendation

1. Keep `temp_env_all_18k/` as the primary leaderboard and ADR 0027 as the
   contract.
2. Report PRnet's 5k result as a **sensitivity finding**, not a new winner: the
   mean baseline is beatable on an HVG axis, but not on the full transcriptome.
3. If a gene-space-robust choice is needed, prefer the model that is competitive
   on both axes; PRnet is the only entrant that is close to the baseline at 18k
   (−0.0345) and above it at 5k (+0.0243).
4. Do not re-tune on the 5k result. Any follow-up (per-fold HVG, batch-aware
   HVG, or a second 5k selection) is a new pre-registered axis and needs a new
   ADR.

## Reproduce

```bash
cd stage2/03_bench_models/temp_env_all_5k

# gene set + frozen data + parity check
python scripts/build_gene_set.py
python scripts/build_data.py
python scripts/build_data.py --check

# offline contracts
make test

# five CPU baselines on fold 0
bash scripts/run.sh --dry-run

# one fold-0 smoke per real adapter
BENCH_SMOKE=1 bash scripts/run.sh --model prnet --fold 0

# full 60-unit queue
bash scripts/run_queue.sh

# evaluation, run cards, fidelity, summary, 18k comparison
python scripts/evaluate.py --output-dir output --results-dir results \
  --registry configs/models.json --data-dir data
python scripts/summarize_run.py
python scripts/aggregate_fidelity.py
python scripts/make_summary.py
python scripts/compare_18k.py
python scripts/validate_run.py
```
