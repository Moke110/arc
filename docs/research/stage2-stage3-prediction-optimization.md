# Stage 2 预测优化 — what is recorded, what is withdrawn, and the live option space

**Status:** research note (read-only survey, 2026-11-13). Not an ADR; it changes no
contract. Every claim below is cited to a path in this tree, a Git object, or a
`tmp/` audit artifact. Where a number was never measured the text says
*not measured / not recorded* instead of estimating it.

**Question.** For Stage 2 (predicted drug-induced Δ transcriptome), how do we make
the prediction as accurate as possible — under three named sub-areas: model
selection / architecture, gene identity / chemistry / dose handling, and
pretraining–transfer.

**Current live state (verified).** `stage2/02_models_tenv/` is the only model step
(ADR 0063 → ADR 0064). Two entrants only: PRnet and the deterministic `Mean`
reference; accepted matrix 2 folds × 2 seeds. PRnet = upstream
`models.PRnet.PGM` imported unedited from `pert_model/PRnet/source`, fed
`chemical = Linear(2048→1024)(Morgan-2048) * log10(dose_µM + 1)`, Gaussian NLL on
`delta + baseline`, early stopping 5-stale-epoch / 20-above-min / 100-epoch cap.
The rest of Stage 2 was deleted on 2026-11-11 (commit `1543781`) and survives only
in Git.

---

## What has been tried (with outcomes)

### 1. The 18,164-gene cold-drug benchmark (12 entrants, ADR 0027)

Ran 2026-10-21/22, 60/60 units, **five** cold-drug folds, single seed `20261021`,
direction-macro / group-macro AUPRC@200. Read from
`git show f50955e^:stage2/03_bench_models/temp_env_all_18k/REPORT.md`.

| rank | entry | AUPRC@200 | PCC | Spearman | R² |
|---:|---|---:|---:|---:|---:|
| 1 | `delta_mean` (global mean Δ) | **0.1989** | 0.2291 | 0.1899 | 0.0198 |
| 2 | `delta_mean_dose` | 0.1948 | 0.2272 | 0.1860 | 0.0171 |
| 3 | `perturbnet` | 0.1861 | 0.2816 | 0.3288 | −1.3878 |
| 4 | `prnet` | 0.1643 | 0.2681 | **0.3840** | −0.7099 |
| 5 | `biolord` | 0.0974 | 0.1077 | 0.2371 | −5.4742 |
| 6 | `ridge_fp` | 0.0921 | 0.0761 | 0.0584 | −0.7552 |
| 7 | `cpa` | 0.0910 | **0.3151** | 0.3457 | **0.0376** |
| 8 | `nearest_drug` | 0.0807 | 0.0635 | 0.0502 | −1.5777 |
| 9 | `cmonge` | 0.0732 | 0.0366 | 0.1121 | −12.6827 |
| 10 | `chemcpa` | 0.0608 | 0.2704 | 0.2866 | −0.0863 |
| 11 | `squidiff` | 0.0212 | 0.0003 | 0.0023 | −4.9e5 |
| 12 | `no_change` | 0.0110 (= prevalence) | — | — | −0.0195 |

Findings that still stand: **(a) no entrant beat the global-mean baseline** on the
pre-registered metric (PerturbNet's paired CI was −0.0248…−0.0007, i.e. below);
**(b) models disagree across metric families** — CPA is 7th on AUPRC but best on
PCC/R²/RMSE, PRnet is best on Spearman; **(c) Squidiff did not converge** at
18,164 genes (valid ε-MSE 0.976 vs 1.0 random) so its row is a failure, not a
ranking; **(d) `delta_mean_dose` did not beat `delta_mean`** (0.1948 vs 0.1989) —
a per-dose mean predictor bought nothing over one global vector.

### 2. The 5,000-gene HVG sensitivity run (ADR 0028) — the ranking flips

Same wells, splits, entrants, seeds; only the gene axis changes (Seurat
`highly_variable_genes` on `log1p_cp10k`, finite mean ∈ [0.0125, 3], ranked by
`dispersions_norm`, sum of ranks over three CM views, 5,000 kept). Chance floor
rises 1.1% → 4.0% at k=200. Recorded outcome
(`git show f50955e^:stage2/03_bench_models/README.md`):

> **PRnet beats `delta_mean` (+0.0243, paired CI +0.0113…+0.0367), while the 18k
> best entrant PerturbNet falls behind.**

This is the single most important *methodological* fact in the repo's Stage 2
history: a 3.6× smaller gene axis re-orders the leaderboard. Any architecture
claim whose effect size is ≈0.02 AUPRC is inside that artefact band.

### 3. The rebuilt tenv shared-Morgan comparison (ADR 0042) — ran to completion, publication retired

`stage2/03_modes_tenv/` froze one Morgan-2048 representation, an adapter per
predictor, 3 cold-drug folds × seeds 17/42/2026, and admitted **five** predictors
through real full-axis train/save/reload/predict admission: PRnet, CPA, chemCPA,
CMonge, PerturbDiff. Final state recorded in
`git show 1543781^:stage2/03_modes_tenv/README.md`: **63/63 formal slots
completed** (CMonge fold 0 seed 17 needed the ADR 0042 supplemental attempt after
its original 14,400 s budget timed out; the supplement stopped by validation
patience at iteration 2,150).

Its leaderboard is **not in the working tree**: ADR 0059 retired the retrieval
scores, then ADR 0063/0064 replaced the metric family and the step was deleted.
The README states: *"Previous retrieval scores, metric figures, and historical
copies were removed; Git retains their version history."* So the five-model tenv
comparison exists, completed, and its numbers are reachable only through Git
history — **no accepted per-model tenv number survives today**.

### 4. PerturbDiff standalone on Tahoe and on tenv (ADR 0046)

Nine random-init fits + nine held-out predictions per dataset (folds 0–2, seeds
17/42/2026), two deterministic references per fold. Completed; two Tahoe
prediction retries after a power cut, recorded in
`execution/formal_record/tahoe/recovery.json` (deleted tree). The scores are
retired from the tree with the same ADR 0059 removal; only the training/terminal
provenance was kept.

Its one durable design finding is in ADR 0046: *PerturbDiff's nine recorded best
validation checkpoints sit at 6,000–13,000 steps and all nine runs stopped at
30,000*, i.e. the validation cadence + 20,000-step floor + patience-10 **was the
binding constraint**, not convergence. That is why `04` exists. Whether the
retuned schedule actually improved held-out scores: **not recorded** in the
working tree.

### 5. Tahoe → tenv transfer, four independent attempts, all negative

| attempt | contract | measured outcome |
|---|---|---|
| PRnet pretrain bring-up (ADR 0031) | 65,318 samples, plate split, 68 epochs, best valid NLL −2.628080 @ epoch 47, 956 s, 15,150 MiB | pipeline proof only; **declared deviation**: raw-count control input vs log target; no transfer claim (`git show 1543781^:stage2/05_prnet_pretrain/README.md`) |
| PRnet formal transfer (ADR 0032) | 5 drug folds, 3 seeds × 3 LRs × 2 inits = 90 candidates, 30 test units, direction AUPRC@200 | **negative**: Tahoe-init **−0.0063** vs random-init (95% CI −0.0074…−0.0052); both far below `delta_mean` **0.1645** (`git show cf3225e^:stage2/06_prnet_transfer/README.md`) |
| PerturbDiff transfer (ADR 0047) | warm-start / cold-cm-drug / cold-drug regimes, CM vs Tahoe vs Transfer vs Mix, 36 learned slots | **Transfer below CM in all nine paired seed comparisons**; mean paired difference −0.02822 (warm-start), −0.02262 (cold-cm-drug), −0.02566 (cold-drug) (`git show 72b4317:stage2/05_pdiff_transfer/README.md`) |
| PerturbDiff transfer optimization (ADR 0048/0049) | five one-at-a-time interventions, 15 target fits, cold-cm-drug fold 0 | random dose reset **+0.0014** AUPRC @1 µM; continuous log-dose encoder **+0.0017** @1 µM; denoiser-LR halved **+0.0013** @10 µM; adaLN-only worst on every dose; **original Transfer had the best mean PCC at every dose** (`git show 731f5db:stage2/06_pdiff_transfer_optimize/README.md`) |

Reading the four together: the deficit is ~0.02–0.03 AUPRC and the five
micro-interventions recover ~0.0015, i.e. **two orders of magnitude smaller than
the gap they were meant to close**. Transfer is negative on two different model
families (PRnet and PerturbDiff), so it is not one architecture's bug.

### 6. Transcriptome reconstruction (ADR 0052–0056) — measured, but not a predictor

`stage2/07_autoencoder/` compared reconstruct-from-itself over five matrix models
with a fixed 18,164 × 256 shared expression interface. Accepted numbers
(`git show 1543781^:stage2/07_autoencoder/README.md`):

| model | tenv PCC | tenv rel. sq. error | Tahoe PCC |
|---|---:|---:|---:|
| `mean` | 0.18690 | 1.00000 | 0.97104 |
| `PCA_1024` | **0.89992** | **0.11612** | **0.98939** |
| `MLP_1024` | 0.88008 | 0.13682 | 0.98331 |
| `AAEc_1024` | 0.70936 | 0.31154 | 0.96885 |
| `AAEcs_1024` | 0.47293 | 0.50106 | 0.96273 |
| `AAEcs_1024+64` | 0.45064 | 0.51289 | 0.95782 |

PCA beat every neural encoder; **both attention additions made reconstruction
worse**. This is reconstruction, not perturbation prediction — no Δ, no drug, no
cold-drug test — so it supports "the 18,164-gene axis is highly compressible"
(≈1,024 of 18,164 dims suffice for reconstruction) and nothing stronger.

### 7. XPert gene-axis / compression ablation (ADR 0057–0060) — never concluded

Round 1 was reduced on 2026-11-11 to two arms (G512, MLP512) × folds 0/1 × seeds
17/42 = eight runs, with a frozen plan and policy. The user then cancelled the
remaining runs, and the READMEs state an incomplete matrix must not publish final
figures. `git show 1543781^:stage2/09_xpert_optimize/README.md`:

> **"A failed run stops the queue… an incomplete matrix must not publish final
> figures or close #144."**

No accepted Round-1 result exists. The tree retains only the *design* (Δ head with
MSE+Pearson, Tahoe→tenv two-phase, native-axis continuous metrics, `AUPRC@full`
retrieval).

### 8. The current accepted comparison (ADR 0063 → ADR 0064)

`stage2/02_models_tenv/README.md`, `log/166-acceptance.json`,
`output/evaluation/per_run.csv`. 2 folds × 2 seeds, 918 test groups of which 169
are null, truth `padj < 0.05` with direction `sign(log2FC)` (the `|log2FC| ≥ 2`
gate was dropped by ADR 0064).

| entry | fold/seed | PCC | Spearman | RMSE | R² | Reliability@d* | changed-group Reliability | d* |
|---|---|---:|---:|---:|---:|---:|---:|---|
| PRnet | 0/17 | 0.1917 | 0.2864 | 0.2299 | −0.2940 | 0.19608 | **0.00000** | strictest k=100 |
| PRnet | 0/42 | 0.1371 | 0.2715 | 0.2609 | −0.8487 | 0.19608 | 0.00000 | k=100 |
| PRnet | 1/17 | 0.2144 | 0.2729 | 0.2205 | −0.1733 | 0.16558 | 0.00000 | k=100 |
| PRnet | 1/42 | 0.1659 | 0.2755 | 0.2491 | −0.6674 | 0.16558 | 0.00000 | k=100 |
| Mean | 0 | 0.1934 | 0.1717 | 0.2070 | **0.0321** | 0.00705 | 0.00877 | k=32 |
| Mean | 1 | 0.1961 | 0.1754 | 0.2067 | 0.0308 | 0.00883 | 0.01067 | k=31 |

Two facts matter. **PRnet's `Reliability@d` peaks at the strictest grid point with
a mean window of 0.03 genes** — the peak is an almost-empty-table solution, and its
changed-group Reliability is exactly 0. **PRnet beats `Mean` on Spearman
(0.27 vs 0.17) but loses badly on R²** (negative vs +0.03) and roughly ties on PCC.

### 9. The structural probe (the load-bearing measurement)

Measured in the Lead's session and **not yet committed anywhere**:
predicted Δ effective rank ≈ **4.0** vs **14.9** for measured Δ; cross-dose
direction correlation **0.945** for predictions vs **0.313** for measured Δ.

Independently re-checked here (read-only, `output/runs/prnet/fold_0/seed_17/
prediction.npy` + `output/inputs/fold_0/seed_17/test_cases.csv` +
`layers/delta` of `training-corpora/tenv.h5ad`, bounded 300-well slabs):

| statistic | predicted Δ | measured Δ |
|---|---:|---:|
| sample-space effective rank (participation ratio, 300 test wells) | **4.03** | **15.67** |
| group-space effective rank (459 test groups, all wells) | **4.02** (90% var in 10 dims, 99% in 18) | not measured here |
| cross-dose direction correlation, same drug, median over 218 within-drug dose pairs | **0.979** | **0.138** |

So: 459 cold-drug test conditions are represented by ≈4 effective directions, and
changing the dose of a drug barely changes the predicted *direction*. The
difference between 0.945 and 0.979 and between 0.313 and 0.138 is
sampling/aggregation (well-level vs group-level); both estimators agree on the
order of magnitude and on the conclusion.

### 10. The uncertainty / p-value audit (rejected route)

`tmp/prnet_uncertainty_audit/REPORT.md` (fold 0/seed 17, all 1,401 test wells):
the variance head is **well calibrated as a predictive sd** (coverage 0.9486 at
|z|≤1.96, 0.9450 over all wells), but the noise vector is not an uncertainty
source (median σ_draw 0.0015 vs median σ_head 0.1962, ratio ≈108; 40.4% of
draws exactly 0; z_draw coverage 0.0119), and the per-gene z-test yields **0 genes
at BH<0.05** for the median 3-well group, with 83.6% of all 81 significant-call
tails coming from **clipped predictions**. Verdict: a predicted-effect p-value is
not a viable substitute for the magnitude threshold. Its byproduct is the most
important architectural fact in this note:

> **40.4% of (well, gene) mean-channel outputs are exactly 0** because the decoder
> ends in a ReLU, so `Δ̂ = −baseline`, i.e. the model predicts a *complete
> knockdown* for 40% of its output.

---

## Withdrawn / untrustworthy conclusions

| conclusion | why it is withdrawn | authority |
|---|---|---|
| The ADR 0063 indirect family (`Precision@k`, `Spearman@k`, `max_k@P80S40`, `Recall@max_k`, `\|log2FC\|≥2` gate) | population chosen by the model being scored: `max_k@P80S40`/`Recall@max_k` defined for 0–4 groups per run, and the PRnet/`Mean` defined-group sets were **disjoint** in every fold/seed; `Precision@k` was NaN (not 0) on an empty gate; `reasons['Spearman@k']` mislabelled `k_below_3`. A metric whose population depends on the entrant cannot rank models | ADR 0064 "Why the @k family is retired" |
| Every Stage-2 retrieval/AUPRC score from steps 03, 04, 05, 06, 09 (incl. the 63/63 tenv matrix and all PerturbDiff transfer tables in the tree) | ADR 0059 collapsed retrieval to three AUPRCs and had the superseded scores and figures **removed from the working tree**; ADR 0063 then withdrew AUPRC reporting entirely (undefined on the 273 null groups, i.e. exactly the population the new contract must score) | ADR 0059, ADR 0063, ADR 0064 "AUPRC reporting is withdrawn" |
| `no_change` as a baseline row | a constant-zero Δ leaves the ranking family with no defined ordering (`Precision@k`=0, `Spearman@k`/`max_k` undefined) and the correlation family undefined; reporting it as a bar misrepresents an undefined entry as a bad score | ADR 0063 "Baselines" |
| "The 18k leaderboard ranks the model families" | the same entrants on the paired 5,000-gene axis re-order: PRnet moves above `delta_mean` (+0.0243, CI +0.0113…+0.0367) while the 18k-best PerturbNet falls behind. Differences ≈0.02 AUPRC are not robust to gene-axis choice | ADR 0028, `03_bench_models/README.md` |
| "PRnet's Reliability shows a 0.18 vs 0.007 win over Mean" | under the changed-group view the accepted record already reports PRnet **0.00000** on all four runs vs `Mean` 0.00877/0.01067; the earlier advantage was produced by the 169 null groups rewarding an empty table. Under the proposed A2 rescoring the two entries land at **0.0107 (PRnet) vs 0.0098 (Mean)** — Lead-measured this session, **not yet a committed record** | `per_run.csv`, README "Reading the result" |
| Predictions from the p-value / noise-draw route | head sd is a *predictive residual* sd, not the SE of the fitted mean; the noise vector is ~108× too small with 40% exact zeros; the median group yields 0 significant genes; 84% of call tails are ReLU-clipped zeros | `tmp/prnet_uncertainty_audit/REPORT.md` |
| The old v1 Tahoe training release | independently rejected: raw-source replay failed (`ValueError: tahoe raw source replay differs: X row48`). Cause: the old source aggregator dropped marker entries then indexed compacted arrays with the original offsets, misattributing and dropping counts. Disposition `REJECTED_SUPERSEDED`; the only usable Tahoe release is `training-corpora/tahoe_bulk.h5ad` | `stage2/01_bases/log/155-gates.md`, `log/records/tahoe_training_v1_rejection.json` |
| "Squidiff is a weak model" | 90× gene-space extrapolation (200 → 18,164) with the author's ε-MSE never leaving 0.976; a documented convergence failure, not a quality claim | `.../REPORT.md` §Squidiff |
| Any comparison of the ADR 0031 bring-up NLL (−2.628080) with another run's loss | Gaussian NLL is a log-likelihood: negative throughout, sign meaningless, not comparable across corpora/rows/splits; and the run used raw-count controls against log targets | ADR 0031 §2 and "Recorded deviation" |
| "Transfer works if you tune the dose bridge / learning rate" | the five ADR 0049 interventions move AUPRC by +0.0013…+0.0017 against a −0.023…−0.028 deficit, and leave original Transfer best on PCC at every dose | `06_pdiff_transfer_optimize` README |

**Untrustworthy by construction, not withdrawn:** the direct metrics in the current
step are single-well statistics over the 18,164-gene axis and carry TempO-Seq
well-level noise (README "Limitations"); the four run samples are a spread
indicator, not an interval; `d*` is an in-sample peak.

---

## Open options by sub-area

### Architecture

**What actually causes the collapse, per code:**

1. **Dose is one scalar multiplying one vector.**
   `stage2/02_models_tenv/execution/script/prnet_fit_predict.py` (`JointPRnet.forward`):
   `chemical = self.adapter(fingerprints) * torch.log10(dose[:, None] + 1.)`.
   For a given drug the adapter output direction is fixed, so the *entire* dose
   response is motion along a 1-D ray. With tenv's 0.2/1/10 µM the multiplier is
   0.079/0.301/1.041 — a 13× magnitude knob and **no shape knob**. This is a
   sufficient explanation of cross-dose direction correlation ≈0.98.
2. **A 64-dimensional conditioning code per sample.** `config.json` sets
   `drug_dimension: 1024, comb_dimension: 64, z_dimension: 64, noise_dimension: 10`;
   upstream `PGM` maps the 1024-d chemical to `c_dim=64` through `PAdaptor`, feeds
   `[x; c]` to the encoder and `[z(64); c(64); noise(10)]` to the decoder. Whatever
   the chemistry and dose say, the decoder sees a 64-d additive code — the drug
   effect is a global shift, not a gene-structured one.
3. **A hard ReLU on the predicted expression mean.** `pert_model/PRnet/source/
   models/PRnet.py` `PDecoder.forward` splits the 2·x_dim output and applies
   `self.relu(recon_x[:, :dim])` to the mean half. Measured consequence: 40.4% of
   (well, gene) outputs exactly 0 → predicted full knockdown, and 84% of
   significant-call tails come from clipped entries.
4. **A Gaussian NLL objective on 18,164 outputs dominated by baseline variance**,
   whose conditional mean under no cold-drug signal is the marginal (mean) effect —
   which is exactly what a 4-dim predicted Δ looks like.

**Options, cheapest first.** Cost anchor: a PRnet fit+predict run costs
**68.07–105.38 s CPU wall, ≤2000.8 MiB peak RSS** (no GPU was available this boot;
`nvidia-smi` blocked), so a full 2 folds × 2 seeds matrix is *minutes*, not hours.
Compute is not the constraint — design is.

| id | change | minimum edit | cost | what it can and cannot fix |
|---|---|---|---|---|
| A1 | **de-clip the output**: identity (or softplus on `baseline+Δ`) instead of ReLU on the mean half; optionally add a per-gene bias | 1 line in a *copied* seam module; upstream source stays untouched | ~free (same matrix cost) | removes 40% exact zeros and the "predict full knockdown" failure; should move R² and Spearman. Cannot raise the sample-space rank on its own |
| A2 | **dose as separate conditioning**: concat `[adapter(fp) ; dose_mlp(dose)]` (or `adapter(fp) + dose_mlp(dose)` in a wider code) instead of the rank-1 product | one forward function in `prnet_fit_predict.py` + config key | ~free | gives the model a *shape* knob per dose; the only option that can move cross-dose direction correlation. At 3 doses it is an identified 3-level effect, not a curve |
| A3 | **wider/channel-expanded conditioning**: raise `comb_dimension` 64→256 and/or `drug_dimension` | config only | ~free | tests whether the 64-d code is the rank ceiling. Confounded with A2 if run together |
| A4 | **low-rank gene decoder**: predict `Δ̂ = G·h(conditioning)` with `G` (18,164 × r), r ∈ {16, 64} | new head replacing the 18,164-way linear readout; ~1.2 M params at r=64 | cheap (CPU minutes) | injects gene-gene structure and shrinks the output parameter count; measured Δ's own rank is 15.7, so r ≈ 16–64 is the right order. It does not add conditioning information |
| A5 | **Δ-space loss**: train on Δ with MSE+Pearson instead of Gaussian NLL on `delta+baseline` | config/seam loss swap | cheap | aligns the objective with the scored quantity; precedent exists in-repo (Step 09 "Δ head with MSE and Pearson", ADR 0057). Note MSE still regresses to the conditional mean, so it cannot by itself create cold-drug signal |
| A6 | **gene-embedding output layer** (predict a coefficient against per-gene embeddings) | new head + embedding table | moderate | parameter sharing across genes; largely subsumed by A4 unless embeddings carry biology (see gene identity) |
| A7 | **SetTransformer / attention autoencoder backbone** (Latent XPert-style: 512 latent queries × 128 channels, frozen decoder, Δ backprop) | needs an AE trained first; GPU; two-phase training | expensive (the retired Step 07 fit was hours-scale and required a sealed split) | the only route that decouples "represent the transcriptome" from "predict the perturbation"; measured reconstruction says PCA/MLP already win at that job, so the AE buys little for reconstruction and its benefit for Δ is unmeasured |
| A8 | **non-architecture levers that were never tried**: capacity ablation (state latent 128 / hidden 256) and environment-compound ablation — both were written into ADR 0032 as tickets #10/#11 and explicitly **deferred, never run** | config | cheap | unknown; ADR 0032 says the round concluded negative before they were needed |

**What the repo cannot tell you yet:** whether the rank ceiling is *dose* (A2) or
*conditioning width* (A3) or *the loss pulling to the mean* (A5). A1/A2/A3/A5 all
fit in one small factorial and each costs minutes; the honest design is to run the
factorial rather than pick a favourite.

### Gene identity

**Current:** a fixed ordered 18,164-gene axis (`stage2/01_bases/output/
gene-identity/genes.csv`, `gene_index` 0–18163, ADR 0061), **no gene embeddings, no
gene-gene structure**; identity enters only as output position.

**Precedents in this repo, and what they are worth:**

- ADR 0057 gave XPert a **fixed random 256-float identity vector per gene** (seed
  20261111; HDGFL3 handled separately) — cheap, no biology, no measured benefit.
- ADR 0056 built a **learned 256-channel gene-position lookup** with a 1→64→256 GELU
  shared generator/restorer (33,281 frozen params) and it worked as a *scaffold*
  (uncompressed round-trip MSE 1.05e-07 Tahoe / 4.37e-07 tenv). That artifact was
  retired with Step 07: `output/shared-expression/` today holds only the
  **expression** matrices (4,478 × 18,164 and 67,018 × 18,164), not the lookup.
- `stage2/01_bases/output/api/bases/hvg.py` + `output/tenv-hvg/` provide the only
  other per-gene structure in use: Seurat normal-dispersion rankings, full-tenv
  descriptive plus **train-only** per-fold rankings (ADR 0061). These are selectable
  gene sets, not embeddings.
- **ESM2 protein embeddings are available as a tool, not as data.** Biomni exposes
  `generate_gene_embeddings_with_ESM_models` (Ensembl IDs → per-gene vectors), and
  the data lake is present at `/home/chang/.local/share/biomni/data/`, but **no
  precomputed ESM embedding table exists in this repo** (`find … -iname '*esm*'`
  empty). Using it would mean fetching 18,164 sequences and running ESM2 — a new
  external dependency and a new admission decision.
- Pathway/module tokens have raw material — `stage1/xpert-fda-miri-v3/input/
  MIRI_pathways.json` and MSigDB in the Biomni lake — but **no Stage 2 pathway
  feature has ever been built or scored**.

**What is justified at n ≈ 459 parents** (tenv: 459 `drug_parent`, 1,377
parent-dose groups, 3 doses, ~3 wells/group; verified from
`training-corpora/tenv.h5ad` metadata):

1. **A low-rank gene decoder (A4) is justified and cheap** — it is a
   parameter-sharing device, not a biology claim, and it matches the measured rank
   of real Δ (15.7).
2. **A gene-identity embedding layer (fixed random or learned) is justified only as
   parameter sharing**; with 459 conditions the model cannot learn 18,164
   independent gene responses, so sharing is the point.
3. **ESM2 / pathway embeddings are not justified yet**: the measured collapse is in
   the *sample* space (rank 4 across 459 conditions), not the gene space. Gene-side
   features cannot create between-condition variation. Cost is real and the
   expected effect is zero-on-the-binding-constraint.
4. **5,000-HVG re-run is not a new idea**: it ran historically (ADR 0028) and it
   *changed the ranking*. If the axis is revisited, that is a declared sensitivity
   axis with its own ADR, not an optimisation knob.

### Chemistry

**What is used now:** one shared representation for every entrant — RDKit Morgan,
radius 2, 2048 bits, chirality on, binary, cached per `drug_key` from the corpus's
normalized `parent_smiles` (`02_models_tenv/AGENTS.md`); `Linear(2048→1024)` adapter
trained jointly with the predictor on training responses only (ADR 0042; ADR 0063
"the same seam the retired adapter used"). One-hot and KPGT were considered and
**withdrawn** as the input policy; "a drug encoder trained independently outside
perturbation data is a separate model, excluded from this comparison" (ADR 0042).

**What actually exists in the tree** (verified, memory-safe header reads):

| representation | where | shape | coverage caveat |
|---|---|---|---|
| Morgan-2048 binary | derived per corpus, `output/inputs/*/morgan.npy` | 2048 | 459/459 tenv parents (corpus builds only with a resolvable structure) |
| HG biological embedding (HeteroGNN, contrastive drug–gene graph) | `stage1/xpert-fda-miri-v3/input/HG_drug_embeddings.npy`, `drug_emb_model/HG/` | 8,981 × 256 float32 | **L1000-indexed** (`drug_idx < 8981`); author reference, no training route left (library dormant, ADR 0008) |
| KPGT graph-transformer embedding | `drug_emb_model/KPGT/` (894 MB weights + 35 MB cached tensors for LINCS2020/GSE92742/FDA_drugs) | 512-d base output (ADR 0042 notes 2,304-d base) | LINCS/FDA-indexed; upstream commit unknown |
| **UniMol 3D atom features** | `stage1/xpert-fda-miri-v3/input/all_drugs_unimol_arr.npy` | **(11,104, 122, 514) float64, 5.57 GB** | L1000 `drug_idx`-indexed; a SMILES→idx join for the 459 tenv parents is **not measured** |
| ChemBERTa / MolFormer | **not in this repo** — only referenced inside third-party sources (`pert_model/ExPO/source/` `chem_encoder.py`, `pert_model/PerturbDiff/source/` `covariate_encoding.py`). `drug_emb_model/` contains **only HG and KPGT** | — | not an available asset |

**What Stage 1 XPert used:** the drug encoder consumes the **author HG 256-d
biological embedding** (byte-parity gated against `HG_drug_embeddings.npy`) *plus*
the frozen **UniMol feature array** indexed by `drug_idx`
(`stage1/xpert-fda-miri-v3/01_drug_embed/drug_embeddings.py` docstring; ADR 0007
"author HG reference byte-identical + precomputed UniMol array"). So Stage 1 is
chemistry-richer than Stage 2 by exactly the two representations Stage 2 withdrew.

**Decision-relevant trade-offs.** (i) Morgan is the only representation with
complete coverage of both corpora and no external dependency; (ii) HG/UniMol/KPGT
are all **L1000-row-indexed**, so adopting one requires a SMILES→L1000 join whose
coverage for 459 tenv parents and 352 Tahoe parents is **not measured** — that
measurement is one afternoon of table work and should precede any run; (iii) the
measured failure mode is *within-drug* (cross-dose correlation ≈0.98), so a
chemistry swap cannot fix it; (iv) ADR 0042 forbids changing the representation
inside the comparison without a new decision.

### Dose

**Current:** `chemical = adapter(Morgan) · log10(dose_µM + 1)` — one scalar. For
tenv that is 0.079 / 0.301 / 1.041 at 0.2 / 1 / 10 µM.

**Measured evidence that this is binding.** Cross-dose direction correlation is
**0.979** for predictions vs **0.138** for measured Δ (my 300-well re-check;
Lead's all-well figures 0.945 vs 0.313). The model cannot express a dose-dependent
*shape*, only a dose-dependent magnitude, and its predicted magnitudes are also
nearly dose-invariant. Independent corroboration from the historical benchmark:
`delta_mean_dose` (per-dose mean Δ) scored **0.1948** vs `delta_mean` **0.1989** —
even an oracle-style dose-stratified mean did not beat one global vector on the
cold-drug metric (this bounds how much *mean-level* dose structure is worth; it
does not measure dose shape).

**Alternatives and what the data can support:**

- **D1 — explicit dose embedding concatenated with the drug code** (A2 above). With
  three dose levels this is a 3-row table; identifiable, cheap, and it removes the
  rank-1 constraint. Preferred first move.
- **D2 — two- (or r-) component dose response per drug**: `Δ̂(d) = u_drug + s(d)·v_drug`,
  i.e. a dose-scaling of a *learned second* direction rather than of the whole code.
  Identifiable at 3 doses only as a 3-point family; a monotone prior in log-dose is
  a modelling assumption, not a measurement.
- **D3 — continuous log-dose encoder with a shared trainable projection** (the
  ADR 0049 `log` arm: `log10(d / 0.05 µM)` → linear → dose-embedding dim, separate
  vector at 0 µM). **This exact design was already run** on PerturbDiff, refit the
  Tahoe source under the same recipe, and gained **+0.0017** AUPRC at 1 µM on the
  cold-cm-drug validation split — no material improvement. It also fixes the
  cross-dataset mapping problem (Tahoe 0.05/0.5/5 vs tenv 0.2/1/10 µM, which the
  Step 05 bridge had handled by log-linear interpolation clamped at the ends, ADR
  0047 "Important limitation"). So the *known* continuous-dose encoder is not a
  proven win on this data; the untested variant is concat-with-full-Morgan (D1).
- **D4 — per-dose output calibration** (three calibration vectors applied to the
  prediction). Cheap, purely post-hoc, and would not change cross-dose direction
  correlation at all; useful only if the delivered table is used at one dose.
- **D5 — cross-source dose harmonisation** for any Tahoe-pretrained arm: the two
  corpora's doses do not overlap (0.05/0.5/5 vs 0.2/1/10 µM), so any transfer arm
  must carry a declared dose bridge; ADR 0047 records that as the priority fix and
  ADR 0049's `log` arm is the only implementation of it.

### Pretraining–transfer

**What the experiments concluded:** three independent negatives (§5 above). PRnet:
Tahoe pre-training makes cold-drug Δ *slightly worse* than random init (−0.0063)
and both are far below the mean baseline; PerturbDiff: Transfer is below CM-only in
all nine paired comparisons (−0.022…−0.028); the five optimisation interventions
recover ≈0.0015. The bring-up run itself (ADR 0031) explicitly claims only that the
pipeline runs, and its control input was raw counts while its target was
log1p-CP10K — a declared fidelity break that was never re-run at parity.

**What now exists** (verified from `training-corpora/manifest.json` and h5ad
metadata only): `tahoe_bulk.h5ad` = **65,318 treated samples × 18,164 genes**,
1,400 matched controls, **50 cell lines**, **352 drug parents**, **1,056
parent-dose groups**, doses **0.05/0.5/5 µM**, layers `control`/`delta`, obs carries
`cell_line_id`, `drug_parent`, `dose_um`, `matched_control_ids`. It is the only
usable Tahoe release (the v1 build is `REJECTED_SUPERSEDED`, §Withdrawn). Row-level
counts and exclusions: `tahoe_bulk_exclusions.csv` (300 dropped for missing
structure).

**Protocol options:**

1. **Pooled multi-source training ("Mix")** — implemented for PerturbDiff (ADR
   0047: pooled Tahoe+tenv training groups, union dose vocabulary, no source
   balancing or reweighting, checkpoint selected on tenv validation) and **never
   run for PRnet**. It is the cheapest transfer protocol and the only one whose
   arm design already exists in the repo's history. Its PerturbDiff outcome:
   **not recorded in the working tree** (scores removed by ADR 0059).
2. **Pretrain → finetune with an explicit dose parameterisation (D3/D1)** — the
   only transfer variant whose failure mode has been diagnosed (dose-bridge
   interpolation/clamping), and the only one where a fix was implemented (ADR 0049
   `log`). Measured gain: +0.0017 AUPRC on 1 µM.
3. **Source re-weighting / domain adaptation / cell-line conditioning** — never
   attempted. The Tahoe corpus carries `cell_line_id` (50 lines) and tenv has one
   cell identity (`hiPSC-CM`, ADR 0061), so a cell-line conditioning input is
   available and untested. **Not measured.**
4. **Transfer only the shape prior, not the signature**: e.g. pretrain a
   gene-factor (A4) basis on Tahoe and freeze it, then fit only the conditioning
   head on tenv. Never attempted; it is the natural pairing of A4 with the transfer
   question.
5. **Don't transfer.** The honest reading of ADR 0029/0031's domain gap — 50 cancer
   cell lines, one timepoint, 31/525 drug overlap by the corpus's desalted-parent
   rule (42 by the union of three identity routes) — plus three measured negatives.

**Also relevant and unmeasured:** the Tahoe corpus has **1,056 parent-dose groups vs
tenv's 1,377**, so Tahoe is *not* a larger condition-level corpus even though it has
15× the profiles; its extra size is cell-line × well replication. Any "more data
should help" argument has to survive that.

---

## What is blocked on the metric fix

The proposed **A2** fix is: per-group normalised threshold `λ = d / max|Δ̂_group|`
(scale-free per group instead of one run-wide `M`), macro average over **changed
groups only**, and `@λ` naming. It keeps ADR 0064's truth side (`padj < 0.05`,
direction from `sign(log2FC)`, 1,377 groups / 918 test groups / 169 null).

**Blocked on A2 (do not spend compute until it lands):**

- Any claim of the form "architecture X improves predicted-Δ accuracy" **measured
  by the `@d`/`Reliability` family**. Today the family's peak sits at the strictest
  grid point with a mean window of 0.03 genes, i.e. it rewards an empty table; a
  changed-group-only macro removes the 169 null groups that produce that incentive.
- Any **cross-entry** comparison of Reliability, any "delivered gene table" claim,
  and any statement about `d*`/table size as a product (README: "a consumer of a
  PRnet checkpoint must treat table size as part of the intervention").
- Any **model-selection decision** taken on the current peak `Reliability@d*` — it
  is an in-sample maximum by construction and, at `d* = M`, degenerate.
- The A2 rescore itself can be produced **without training**: the accepted
  `--accept-only` path already reloads the eight published runs and rescores twice
  byte-for-byte, so A2 is a driver change plus a replay, not a new experiment.

**Not blocked (can proceed now):**

- **Direct Δ-space metrics** (PCC, Spearman, RMSE, R²) — unchanged by A2. Any
  accuracy claim can be made on these today, and the current record already shows
  the key contrast (PRnet Spearman 0.27 vs `Mean` 0.17; PRnet R² negative vs
  `Mean` +0.03).
- **Structural read-outs**: predicted-Δ effective rank, cross-dose direction
  correlation, per-group window size. These are not part of the `@d` family and are
  reproducible from `prediction.npy` + the corpus Δ in seconds (done in §9).
- **Dose-parameterisation experiments**, whose primary read-out is the cross-dose
  correlation — a structural statistic, not a scored metric.
- **Chemistry-coverage measurements** (SMILES→L1000 join for HG/UniMol/KPGT) —
  pure table work, no metric dependency.
- **Pretraining-transfer paired contrasts**: pairing is by fold/seed on the same
  split, so the *difference* between arms is meaningful under whatever metric is
  current, including direct metrics only.
- **The gene-axis sensitivity question** (ADR 0028 logic) — governed by its own
  ADR, independent of the `@d` family.

---

## Recommended next experiment if only one can run

**One small factorial on the existing PRnet seam — not a new model, not the AE.**
Run the current 2 folds × 2 seeds matrix with exactly **one declared seam change**:

> **PRnet-D: replace the rank-1 dose product with an explicit dose embedding
> (`chemical = adapter(Morgan) ⊕ dose_mlp(dose)`) and de-clip the decoder output
> (identity instead of ReLU on the mean half).**

Why this one:

1. It attacks the two *measured* pathologies rather than a hypothesis: the rank-1
   dose product (§9, cross-dose corr ≈0.98 vs 0.14 measured) and the 40.4% clipped
   zero-expression outputs (§10).
2. It costs **minutes**: the accepted PRnet runs are 68.07–105.38 s CPU each, so
   the full matrix plus rescoring fits in well under an hour, with no GPU, no new
   dependency, no new data, and no change to `pert_model/PRnet/source` (the seam
   lives in `prnet_fit_predict.py`).
3. It is diagnostic whatever the outcome: if effective rank rises above ~6 and
   cross-dose correlation drops toward the measured 0.14, conditioning capacity is
   the constraint and the next spend is A3/A4/A5; if both stay flat (rank ≈4,
   corr ≈0.95+), then conditioning capacity is *not* the limit and no further
   architecture work should be funded before a predictability decomposition that
   asks which condition axis carries any cold-drug signal at all (baseline-only,
   dose-only, chemistry-only deterministic references, which ADR 0063 says needs a
   new decision because a reference cannot be added casually).
4. Pre-register the read-outs **before** the run: direct metrics (PCC/Spearman/
   RMSE/R², unchanged by A2), predicted-Δ effective rank (sample and group space),
   cross-dose direction correlation, and the A2 `@λ` curve only if A2 has landed.
   Do not select on peak `Reliability@d*`.

**Do not** spend this one slot on: a chemistry swap (cannot fix an within-drug
collapse, and needs a coverage measurement first); an ESM/gene-embedding layer (a
sample-space rank problem, not a gene-space one); an AE/SetTransformer backbone (the
retired Step 07 measurement says PCA/MLP beat attention at reconstruction, and the
fit is hours-scale); or a Tahoe transfer arm (three negatives, and the only
diagnosed fix worth +0.0017).

---

## Citations (paths + commits)

**ADRs (working tree, `docs/adr/`):** 0026 CM corpus/views · 0027 benchmark contract ·
0028 5k HVG sensitivity · 0029 Tahoe pseudobulk corpus · 0030 external-data layout ·
0031 Tahoe CM-schema view + PRnet pretrain bring-up · 0032 PRnet formal transfer ·
0033 dataset-corpus merge · 0042 tenv benchmark redesign (+ 2026-11-06/07 amendments) ·
0046 PerturbDiff standalone · 0047 Tahoe→tenv transfer · 0048 transfer diagnosis ·
0049 transfer optimization arms · 0051 Latent XPert · 0052–0056 Step 07 reconstruction ·
0057 XPert landmark reference · 0058 direction-matched PR · 0059 three-AUPRC retrieval ·
0060 Step 09 round-1 rerun · 0061 shared bases and indexed corpora ·
0062 accepted products/figures · 0063 direct+indirect metrics · 0064 @d threshold scan ·
0007/0008 XPert/HG finalisation and library purity.

**Live code and accepted products**

- `stage2/02_models_tenv/README.md`, `AGENTS.md` — entrants, matrix, accepted
  numbers, limitations, retirement boundary.
- `stage2/02_models_tenv/execution/script/prnet_fit_predict.py` — `JointPRnet.forward`
  (`chemical = adapter(fingerprints) * log10(dose+1)`), `predict_deltas`,
  `replay_prediction`.
- `stage2/02_models_tenv/execution/models/prnet/config.json`,
  `bounded_run.json` — `drug_dimension 1024`, `comb_dimension 64`, `z_dimension 64`,
  `noise_dimension 10`, hidden `[128]`, early-stop keys.
- `stage2/02_models_tenv/execution/script/mean_reference.py` — one global training-mean
  Δ vector, copied to every test well.
- `stage2/02_models_tenv/output/evaluation/per_run.csv` — per-run M, d*, Reliability,
  null/changed split, direct metrics.
- `stage2/02_models_tenv/log/166-acceptance.json` — acceptance-v2 record, replay proof,
  evaluation-driver commands, resources.
- `pert_model/PRnet/source/models/PRnet.py` — `PGM`, `PAdaptor`, `PDecoder` (final
  `nn.ReLU()` on the mean half at lines 216–231).
- `stage2/01_bases/output/training-corpora/{manifest.json, tahoe_bulk.h5ad, tenv.h5ad}`
  — 65,318 × 18,164 / 4,206 × 18,164, 1,056 / 1,377 parent-dose groups, dose sets.
- `stage2/01_bases/output/api/bases/hvg.py`, `output/tenv-hvg/manifest.json` — Seurat
  HVG rule, descriptive vs train-only rankings.
- `stage2/01_bases/output/gene-identity/genes.csv`, `output/shared-expression/` —
  frozen axis; expression releases only (no learned gene lookup survives Step 07).
- `stage2/01_bases/log/155-gates.md`, `log/records/tahoe_training_v1_rejection.json`,
  `log/155-inventory.md` — v1 Tahoe release rejection and compact failure evidence.

**Deleted steps, read from Git**

- `git show f50955e^:stage2/03_bench_models/temp_env_all_18k/REPORT.md` — 18k leaderboard,
  paired bootstrap, Squidiff failure, per-drug view.
- `git show f50955e^:stage2/03_bench_models/README.md` — 5k sensitivity result and the
  ranking flip.
- `git show 1543781^:stage2/03_modes_tenv/README.md` — 63/63 slots, admission registry,
  CMonge supplement, retirement of scores.
- `git show 1543781^:stage2/01_setup/README.md` — nine-model bring-up, "no model produced a
  whole-transcriptome prediction at setup time".
- `git show 1543781^:stage2/04_pdiff_hyperparameter/README.md` — PerturbDiff retune scope,
  nine fits per dataset.
- `git show 1543781^:stage2/05_pdiff_transfer/README.md` and
  `git show 72b4317:stage2/05_pdiff_transfer/README.md` — transfer contract and the
  pre-retirement paired differences.
- `git show 731f5db:stage2/06_pdiff_transfer_optimize/README.md` — five interventions and
  their ~0.001 magnitudes.
- `git show cf3225e^:stage2/06_prnet_transfer/README.md` — PRnet transfer negative result
  and deferred ablations #10/#11.
- `git show 1543781^:stage2/05_prnet_pretrain/README.md` — bring-up numbers and the
  raw-count-control deviation.
- `git show 1543781^:stage2/07_autoencoder/README.md` — reconstruction table.
- `git show 1543781^:stage2/09_xpert_optimize/README.md` — cancelled Round 1, Δ-head design.
- Commit `1543781` ("Retain only bases under Stage 2 per authorized scope") — the
  1,317-file deletion that removed these steps.

**Other evidence**

- `tmp/prnet_uncertainty_audit/REPORT.md` (+ `audit.json`) — variance-head calibration,
  noise-vector invalidity, ReLU clipping, rejected p-value route.
- `stage1/xpert-fda-miri-v3/01_drug_embed/drug_embeddings.py` and
  `stage1/xpert-fda-miri-v3/input/` — HG 8,981 × 256, UniMol (11,104, 122, 514)
  float64, idx→SMILES table.
- `drug_emb_model/README.md`, `drug_emb_model/AGENTS.md` — HG and KPGT only; KPGT
  894 MB weights + 35 MB cached tensors.
- `pert_model/README.md`, `pert_model/AGENTS.md` — the fourteen model 项, their licences,
  and the reference-only/excluded list.
- `stage2/AGENTS.md` — step contract, accepted products, acceptance/replacement rules.

**Numbers measured in this survey but not yet recorded anywhere else** (reproducible
from `prediction.npy` + `test_cases.csv` + `layers/delta`): predicted-Δ sample-space
rank 4.03 (300 test wells), group-space rank 4.02 (459 groups; 10 dims = 90% var, 18 dims = 99%),
measured-Δ rank 15.67, prediction cross-dose direction correlation 0.979 vs measured
0.138 over 218 within-drug dose pairs (300-well subsample, fold 0 seed 17). The
Lead's all-well figures (4.0 vs 14.9 rank; 0.945 vs 0.313 correlation) are the same
measurement at a different population.

**Explicitly not established by this survey:** the tenv five-model leaderboard
numbers (completed, then removed); the PerturbDiff standalone/transfer absolute
scores (removed); whether the ADR 0046 retune improved held-out scores; any
Tahoe→tenv dose-parameterisation effect on PRnet; chemistry-representation coverage
for HG/UniMol/KPGT over the 459 tenv / 352 Tahoe parents; whether any gene-embedding
or pathway feature helps predicted Δ; the A2 rescored values as a committed record.
