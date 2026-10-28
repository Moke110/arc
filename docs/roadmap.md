# AICM Roadmap

Program-level milestone index for the AICM (AI cardiomyocyte) drug-repositioning effort.
Terminology lives in `GLOSSARY.md`; decisions in `docs/adr/`; per-feature specs in
`.scratch/<slug>/`. A *Stage* is a program-level milestone — distinct from a workflow
*step* and from a *component* (see `GLOSSARY.md` → Roadmap).

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
- Cardiotoxicity scoring is deferred (`GLOSSARY.md` → Cardiotoxicity score).

## Stage 2 — model setup and selection (active)

**Objective.** Select the open-source whole-transcriptome drug-perturbation model that
should produce the perturbomes feeding future MIRI drug scoring. The objective is
unchanged; the plan is now staged so that model setup is not blocked by benchmark
design.

**Setup (current).** Prove that the nine deployed non-preprint models can run on
their native assets on this machine. Per model: a pinned environment, a load/train
run, one prediction artifact, a run log, and a one-page run card. No cross-model
metrics, no folds, no CM adaptation, and no benchmark design yet.

| model | native asset | smoke target | caveat |
|---|---|---|---|
| PRnet | LINCS checkpoint + 978-gene demo | `test_demo.py` | **setup-only**; 978 is quarantined from Stage 2 evaluation |
| CPA | Kang checkpoint + Kang h5ad | load and predict held-out Kang cells | IFN-β stimulation, not a chemical screen |
| chemCPA | public `ComPert` checkpoints | load and forward | public weights are L1000 977 / 2,000 genes, 188 drugs — **not** 5,000/20,570; real-data inference needs retraining |
| biolord | none published | train a small model on `GSM_new` | per-drug categories, not a cold-drug structural model |
| PerturbNet | LINCS/sci-Plex cINN + ChemicalVAE | load and latent forward | full inference needs Python 3.7 / TF 1.15 |
| scVIDR | none published | synthetic 2,000-gene train + dose predict | dose is log-linear latent interpolation |
| Squidiff | none published | 2,000-gene drug-structure forward | dose is `rFCFP × log10(dose+1)` |
| CMonge | bundled sciplex demo | official demo train + eval | 100-gene demo; categorical dose conditions |
| DECANT | none published | 2,000-gene synthetic forward | no processed cell bags or checkpoint |

Preprint models (X-Pert, DORA, BaiZe, PerturbLDM) are excluded and their code is
not vendored. WAVE and ExPO are reference-only (no dose / 978-only).

**978 quarantine.** PRnet's native demo uses the 978-landmark LINCS matrix. It is
allowed only under `stage2/01_setup/prnet/` and its outputs are labelled
`setup-only`. No later Stage 2 evaluation reads a 978 artifact; the
whole-transcriptome guard still applies to Stage 2 itself.

**Environment.** One dedicated environment per model; the shared `arc` env is not
modified. Existing four keep `arc-stage2-bringup-<model>`; new models use the same
prefix; biolord and CMonge keep their dedicated venv overlays. The exact commands
and resolved package versions live in each run card.

**Layout.** `stage2/PLAN.md`, `stage2/01_setup/<model>/`, and
`stage2/audit/`. The old `cm_bench_v1` benchmark is archived at
`stage2/_archive/cm_bench_v1/` and is not a dependency.

**Exit condition.** 9/9 run cards with a successful command and a saved artifact; any
blocked model recorded with its exact error; no upstream model code edited; a short
setup summary states which models are ready for the selection stage.

**Selection stage — deferred.** After setup, decide which model(s) deserve CM adaptation,
which gene space (978 / 5,000 / 20,570) to evaluate, whether CPA/biolord need chemical
retraining, and what metric/folds to use. No benchmark design is fixed now.

**Status.** Setup completed 2026-10-19; results and selection-stage implications at
`stage2/01_setup/SUMMARY.md`. Plan at `stage2/PLAN.md`; audit at
`stage2/audit/entrant_audit.md`. `cm_bench_v1` is archived, not deleted.
