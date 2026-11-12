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

## Stage 2 — shared foundations

Stage 2 retains only `stage2/01_bases/`: shared gene identity, separate Tahoe
and tenv training corpora, tenv measured truth, training-only HVG selection,
and shared retrieval metrics/reporting. Build and acceptance commands are in
`stage2/01_bases/README.md` (agent routing: its `AGENTS.md`). The other Stage 2
steps and their local assets were deleted at the user's request on 2026-11-11.
Prior benchmark records can be inspected in Git; their ignored predictions and
checkpoints are no longer local inputs. The original consumer-migration scope
is not an acceptance claim for this retained foundation.

The Stage 2 exit is now fixed: **the predicted Δ on a selected gene axis, as
accurate as possible**. Narrowing the frozen 18,164-gene axis to that selected
axis — from tenv HVG plus public myocardial-disease datasets — is being decided
and delivered under the [Stage 2 exit → Stage 3 exit decision
map](https://github.com/Moke110/arc/issues/171). Model, chemistry, dose and
pretraining-transfer optimisation is recorded on that map only, as a decision
card with its supporting research in `docs/research/`.

## Stage 3 — cardiomyocyte phenotype direction prediction (first-wave corpus audited)

Predict the direction of change (+1 / 0 / −1) of one phenotype readout under one
specified experimental condition from the drug-induced perturbome predicted by a
Stage 2 model.

Both stage boundaries are fixed (2026-11-13): the Stage 2 exit and Stage 3 entry
are **the predicted Δ on the selected gene axis**, and the Stage 3 exit is a
**phenotype change direction a simple non-high-throughput wet-lab readout can
detect** (for example cell viability), not an in-vivo endpoint. Input
representation and model design remain open; they are recorded, with their
evidence, on the [Stage 2 exit → Stage 3 exit decision
map](https://github.com/Moke110/arc/issues/171).

The first effort plans the drug–phenotype experimental database: 24–48 h drug
exposure followed by MTT/CCK8, the same exposure followed by hypoxia/reoxygenation
and MTT/CCK8, and other single readouts with at least 500 non-repeat drug–readout
observations under comparable conditions. It does not require paired measured
transcriptomes or process Stage 2 predictions. Scope and counting boundaries are
recorded in [ADR 0034](adr/0034-stage3-phenotype-data-boundary.md).

The canonical decision index is the
[Stage 3 phenotype-database wayfinder map](https://github.com/Moke110/arc/issues/21).
The current source-first strategy ([ADR 0044](adr/0044-stage3-source-first-corpora.md))
selects Heliste 2021 for separate baseline and H-R chemical-perturbation CTG
`cell_viability` corpora at `stage3/01_corpus/output/`. The public
[build/verify workflow](../stage3/01_corpus/README.md) preserves supplied normalized
values, all four doses and source-cell identity. Chemical structures require
frozen resolution evidence; unresolved treatments remain in coverage reports.
Statistics and endpoint labels are outside this corpus contract.

The earlier endpoint-gated first wave is superseded. Its unchanged
[handoff](../stage3/01_corpus/execution/historical/handoff.md) and coverage evidence
remain under `execution/historical/`; the retired implementation remains in Git.
Current implementation is tracked by [spec #52](https://github.com/Moke110/arc/issues/52)
and [ticket #53](https://github.com/Moke110/arc/issues/53). Model views and training
remain later work.
