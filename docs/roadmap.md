# AICM Roadmap

Program-level milestone index for the AICM (AI cardiomyocyte) drug-repositioning effort.
Terminology lives in `CONTEXT.md`; decisions in `docs/adr/`; per-feature specs in
`.scratch/<slug>/`. A *Stage* is a program-level milestone — distinct from a workflow
*step* and from a *component* (see `CONTEXT.md` → Roadmap).

## Thesis (开题)

**利用虚拟细胞技术筛选缺血性心脏病治疗药物** — *Using virtual-cell technology to screen
drugs for ischemic heart disease.* This proposal frames the AICM effort below; the
Stage-1/2 screens target MIRI (myocardial ischemia-reperfusion injury) as the
ischemic-heart-disease transcriptional signature to reverse. Proposal deck:
`presentation/开题.pptx` (20 slides + speaker notes; the deck is now the single
proposal content source — the earlier Markdown outline was retired). Background
literature review: `presentation/my_review.md`.

Proposal framing (decided 2026-10-11): 5 chapters (1 立题依据 / 2 研究目标与内容含创新 /
3 现有研究结果 / 4 后续研究计划 / 5 参考文献); unified spine **模型→药效→闭环**; foregrounds
the **v2** screen as the representative Stage-1 result — **v1/v1.1/v3 are treated as immature
preliminary exploration and are not mentioned anywhere in the proposal**. Sample-scarcity
argument cites Amrute (Nature 2024) + Xie (Science 2026).

## Stage 1 — XPert + DToxS MIRI drug repositioning (complete)

**Method.** Finetune the author-pretrained XPert (L1000) on cardiomyocyte drug-perturbation
data, predict the 978-landmark perturbed transcriptome of each of the 1,004 L1000-contained
FDA drugs, score MIRI transcriptional reversal (WCS+RGES / Zhang / pathway RGES → integrated
`-0.3*RGES_zs - 0.2*Zhang_zs + 0.3*pathway_reversal_zs + 0.025*sum(critical efficacy zs)`),
then filter and rank.

**Screens.** All four share the same 11-step layout, screen the same 1,004 drugs, and are
anchored on the dtoxs-primary test sets. They differ only in finetune corpus:

| workflow | finetune corpus | screened / WCS<0 / ≥6-of-8 critical | role |
|---|---|---|---|
| `xpert-fda-miri-v1` | frozen DToxS-only (521 samples / 36 drugs) | 1004 / 244 / **42** | frozen baseline; canonical finalized screen; archived agentic report |
| `xpert-fda-miri-v1.1` | DToxS-LINCS (662 / 45) | 1004 / 95 / **20** | DToxS-LINCS testbed |
| `xpert-fda-miri-v2` | primary+iPSC union (958 / 45) | 1004 / 130 / **64** | primary+iPSC sibling |
| `xpert-fda-miri-v3` | all-CM union (~4,318 / ~180) | 1004 / 252 / **56** | current frontier |

The rebuilt workflows (v1.1/v2/v3) share byte-identical dtoxs-primary fold test sets under
the unified fold scheme (ADR 0014/0015/0017); v1 predates it. Per-plate normalization and
the union builds are ADR 0016/0017.

**Artifacts.** Per workflow: step 11 `output/` (`integrated_scores.csv`, `ranked_drugs.csv`,
`report.md`), step 05 benchmarks (`dtoxs-primary/`, `data-ablation/`), and `screen.ipynb`.

**Benchmark (step 05).** XPert finetune beats the L1000-pretrain zero-shot anchor; warm-start
> cold-drug. The data-ablation (v2/v3) shows test performance rising with finetune-corpus
scale and saturating — CM data still pays off at the top of the current corpus.

**Wet-lab validation (first experimental check).** Two candidates — trimethoprim (MP, 甲氧苄啶)
and mirtazapine (MT, 米氮平) — were taken into wet-lab experiments: CCK8 cardiomyocyte viability
(safety, and hypoxia/reoxygenation efficacy) and a mouse IR model (%EF/%FS/infarct, plus a non-IR
safety arm). Both were non-toxic in vitro and in vivo; **MT significantly reduced infarct area**
(IR+MT vs IR, p<0.01) while MP did not, and neither reached significance on echo. The first
experimental readout for the screen thus **partially supports MT** and **does not support MP**.
Record: `docs/wetlab-validation-stage1.md`; raw data `wet_data/` (gitignored).

**Known limitations.**
- XPert predicted perturbations are ~97.2% down-regulated DEGs with low inter-drug variance
  (documented in the v3 data-ablation spec) — a structural bias that limits drug discrimination.
- The screen's target is transcriptional MIRI reversal, not an experimental cardioprotection
  endpoint. The first wet-lab validation (above) is small-n, single-dose, and covers only two
  candidates, so the screen still lacks a systematic experimental ground truth.
- Cardiotoxicity scoring is deferred (`CONTEXT.md` → Cardiotoxicity score).

## Stage 2 — CM benchmark of open-source drug-perturbation models (planned)

**Objective.** Model selection: decide which open-source whole-transcriptome model should
produce the perturbomes that feed future MIRI drug scoring. Motivation: Stage 1's XPert is
landmark-only and biased (~97.2% down-regulated DEGs, low inter-drug variance).

**Paradigm.** Stage 2 is **whole-transcriptome end-to-end**: entrants are models that
predict the whole transcriptome; models that only predict the L1000 978 landmarks are
excluded (this excludes XPert and TranSiGen, the two incumbents). No 978-landmark subset is
used anywhere; Stage 1's evaluation is a *method reference only*. Stage 1's landmark line
remains its historical deliverable, but Stage 2 opens a new full-transcriptome line. (ADR
0020.)

**Versions & workflow.** Stage 2 v1 = `workflow/cm_bench_v1/` (core entrants). Stage 2 v2
extends v1 with the phase-2 model families. Repo conventions unchanged: model source lives
in `pert_model/`, CM data in `pert_data/` — both populated by download (`.git` removed) and
then untouched; all scripts and artifacts live in the workflow.

**Environment.** Each model runs in its own pinned conda env (dependency conflicts across
model families are near-certain), named `arc-stage2-<version>-<model>` (e.g.
`arc-stage2-v1-chemcpa`); the workflow drives each model via `subprocess`. The shared `arc`
env stays clean for Stage 1. This is a deliberate exception to the single-runtime convention
(ADR 0021).

**Workflow layout** (`cm_bench_v1`): `01_bench_data` (raw → per-plate Δ → common gene
universe → limma true sets → cold-drug/cold-dose split masks), `02_baselines`,
`03_<model>` (adapter: train on fold train → predict held-out drug/dose), `04_evaluate`
(direction-aware AUPRC, per-dataset/per-axis aggregation, leaderboard, report), plus
`_lib/` and `input/`.

**Entrants (v1 core).** chemCPA, biolord, CPA, scGen, State, plus mandatory baselines
(`no-change`, `mean-effect`, `ridge-fp`, `nearest-drug`). All support unseen-drug
(cold-drug) prediction and output the whole transcriptome. Stage 2 v2 adds generative/OT/KG
families (TxPert, scPRAM, CellOT, MCDP/CFM-GP/LDM-Drug, PRnet) and, in a later tier,
finetuned single-cell foundation models (scGPT, scFoundation).

**Single metric.** Direction-aware **AUPRC recovering the true significant DEG set**
(up- and down-regulated sets scored separately, then averaged). Rationale: the MIRI score
is a pathway-enrichment statistic over the ranked DEG vector, so what matters is whether a
model ranks a drug's *true* significant genes high (with the correct sign) — precision and
recall are both captured. Explicitly not PCC. Threshold-free on the predicted side (the
model only supplies a ranking). Random baseline = the true-set prevalence.

**Model I/O contract.** Input = (the dataset's baseline transcriptome for the cell context,
the drug, the **dose**); output = the perturbation Δ (whole transcriptome). Time is merged.
Because the baseline is part of the input and carries the dataset/context scale, each
dataset uses its own baseline and there is **no cross-dataset harmonization**. Per-plate
baseline = the mean of the control wells on the same plate (Stage 1's ADR 0016 method,
referenced only); raw per-dataset normalization (log-CPM/quantile) is applied before Δ.
Architectures whose native condition is a cell-line embedding rather than the
control-expression vector (chemCPA, CPA) need an adapter — the main engineering work of
`cm_bench_v1`.

**Evaluation axes.** Two generalization axes, both scored by the same metric: (1)
**cold-drug** — hold out drugs across all their doses; (2) **cold-dose** — hold out dose
levels. Cold-dose is feasible only on gse244740 (the only multi-dose dataset: 38 drugs × 8
doses × 3 wells); DToxS-primary/iPSC are single-dose. Per-(drug, dose) true sets in
gse244740 rest on 3 wells, so cold-dose true sets are statistically weaker than the
DToxS cold-drug true sets (16 wells/drug). The **selection criterion is cold-drug AUPRC**
(deployment predicts unseen FDA drugs); cold-dose is a mandatory robustness column, with
N/A for architectures that cannot take dose.

**Ground truth.** Per (dataset, drug, **dose**), time merged: limma moderated one-sample t
on the per-well deltas (treated − per-plate control) vs 0, BH-FDR < 0.05 (+ optional
`|mean Δ|` gate), sign of `mean Δ` = direction; drugs/doses with too few replicate wells are
marked unreliable and dropped from aggregation. Computed on the full-transcriptome gene
space, from the datasets' raw matrices (the 978 `cm_union.h5ad` is not used). Cold-drug
split = 5-fold CV over the drug union (each drug tested exactly once); eligible drugs have
≥4 treated wells in at least one dataset. Cold-dose split = leave-dose-level-out on
gse244740, holding out one interior (interpolation) and one extreme (extrapolation) level.

**Gene universe.** A single common universe = the **intersection** of the participating
datasets' gene sets (gene-symbol/Ensembl mapped). True DEG sets and AUPRC are computed only
in this universe, so the AUPRC random baseline (true-set prevalence) is comparable across
models. No 978-landmark subset anywhere.

**Datasets & drug scope.** No in-HG restriction. The full-transcriptome CM datasets:
dtoxs_primary (~12k expressed genes, 3′-DGE), dtoxs_ipsc (38k, RNA-seq), and gse244740
(44k, ScreenSeq) — together ~**79 drug-like compounds with SMILES** (46 DToxS small
molecules + 42 gse244740 − 9 overlapping). gse262419 (TempO-Seq panel, ~102 mostly
environmental compounds) is excluded unless confirmed otherwise. gse244740 SMILES are empty
in the union build and must be recovered from the raw series matrix.

**Model preparation (compute-constrained).** No entrant is pretrained or trained from
scratch: each starts from a **public pretrained checkpoint** and is adapted to CM data by
**finetuning** on the fold's train drugs. A pre-flight checkpoint-availability audit records
each entrant's public weights and CM-finetune recipe; entrants without a usable public
checkpoint are dropped from v1 and documented, and zero-shot is used only where a checkpoint
exists but no finetune recipe does. A model's benchmark score is its cold-drug AUPRC in that
best-fit form.

**Aggregation.** Primary metric = **macro** AUPRC (per-dataset, datasets equally weighted;
within a dataset, averaged over its (drug, dose) points). Micro (pooled) and per-drug-equal
views are reported alongside as robustness diagnostics.

**Exit condition.** A model leaderboard + a recommended winner + a reproducible benchmark.
Stage 2 does **not** produce a MIRI drug ranking; the full-transcriptome MIRI scorer
(rebuilt 08/09 on the whole transcriptome) is deferred to a follow-on stage informed by the
winner.

**Status.** Spec published at `.scratch/cm_bench_v1/spec.md`; 12 execution tickets published
at `.scratch/cm_bench_v1/issues/` (ready-for-agent).
