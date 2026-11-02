# 06_prnet_transfer — PRnet: Tahoe pretraining and temp_env CM finetuning

This step answers the question Stage 2 still has open after `05_prnet_pretrain`:
**does Tahoe-100M pretraining improve prediction of drug responses in human
iPSC-cardiomyocytes (iPSC-CM) for drugs never measured there?** The treatment
target is the general CM drug response; the held-out unit is the drug, not a
gene set. The study is pre-registered by an approved spec: the decisions, the
metric, the budgets and the negative-result clause were fixed before execution
(`.scratch/stage2-prnet-formal-transfer/spec.md`, gitignored; ticket graph in
`.scratch/stage2-prnet-formal-transfer/ticket-plan.md`).

Status: **in progress** (2026-11-02). The two audits are built and merged
(#1 eligibility rule v1.3 / handoff v2, #2 per-observation QC) and the frozen input
(#3) is built and verified (118 main-set drugs, 951 wells, 125/125 verify); no
training has been executed. Each section below flips to executed only when its
ticket closes with artifacts and evidence.

## The first-round design (fixed)

| | |
|---|---|
| question | prediction of the transcriptomic response of unseen drugs in iPSC-CM; paired Tahoe-init vs random-init comparison |
| scope | 24 h only; Tahoe + temp_env only (no ScreenSeq, no DToxS branch, no GSE247541) |
| CM data | temp_env (`GSE262419`) restricted to the eligibility ∩ QC set; low-response and toxicity/stress wells are kept and stratified |
| gene axis | the frozen 18,164 genes of `02_cm_data` corpus/v1, order preserved |
| model | the author PRnet `PGM`: hidden 128, drug latent 64, state latent 64, noise 10; FCFP4 × log10(dose+1); full-parameter finetune; Gaussian NLL loss |
| split | five folds by `drug_key` (~6:2:2); every dose/replicate of a drug stays together; test and valid once per drug |
| control | CM train uses same-plate/same-vehicle matched controls; validation, test and inference get one fixed training CM consensus `c*`; experimental truth Δ uses matched controls and never enters a model input |
| primary metric | direction AUPRC@200, macro-averaged over drug×dose conditions after within-condition replicate aggregation |
| selection | validation AUPRC@200 (max) drives checkpointing, early stopping and schedule; test never touches selection |
| repetition | 3 seeds; 3 independent Tahoe pretrains reused across 5 folds; 2 inits × 5 folds × 3 seeds × 3 LRs = 90 CM candidates; 30 CM test units |
| baselines | `no_change`, training-mean Δ, plus dose-mean Δ / ridge-fingerprint as auxiliary |

## Layout

Per ADR 0030 §3 this container splits into eleven independently-concluded
substeps, one per published ticket, each carrying `AGENTS.md`, `README.md`,
`execution/` and `output/`: `06.1_eligibility/`, `06.2_qc/`, `06.3_input/`,
`06.4_pipeline/`, `06.5_tahoe_pretrain/`, `06.6_cm_random/`, `06.7_cm_tahoe/`,
`06.8_baselines/`, `06.9_evaluation/`, `06.10_capacity_ablation/`,
`06.11_env_compound_ablation/`. See [`AGENTS.md`](AGENTS.md) for what each owns
and [`docs/adr/`](../../docs/adr/) for the contracts. There is exactly one
eligibility audit, one QC audit, one frozen input build, and one evaluator;
every experiment consumes them rather than re-deriving a split.

## Status table

| substep | ticket | status | evidence |
|---|---|---|---|
| `06.1_eligibility/` eligibility audit | #1 | built and verified — rule v1.3, handoff v2 (118 include / 28 uncertain / 314 exclude) | [`README`](06.1_eligibility/README.md), `output/current.json` |
| `06.2_qc/` per-observation QC | #2 | built and verified — locked rule set | [`README`](06.2_qc/README.md) |
| `06.3_input/` frozen inputs, folds, consensus | #3 | frozen and verified — 118 drugs / 951 wells, 125/125 verify | [`README`](06.3_input/README.md), `output/manifest.json` |
| `06.4_pipeline/` pipeline + smoke | #4 | Phase A built; Phase B real-data smoke green (17/17 steps, GPU peak 2.9 GB, host peak 8.5 GiB) | [`README`](06.4_pipeline/README.md), `output/smoke/smoke_record.json` |
| `06.5_tahoe_pretrain/` 3 Tahoe pretrains | #5 | wrapper built and profiling-verified (75.9 s/full epoch, 11.3 GiB peak); formal sweep running | [`README`](06.5_tahoe_pretrain/README.md) |
| `06.6_cm_random/` 45 random-init CM candidates | #6 | executed — 45/45 (40 ok / 5 capped), LR 1e-3 in all folds, 15 checkpoints selected | [`README`](06.6_cm_random/README.md), `output/record/record.json` |
| `06.7_cm_tahoe/` 45 Tahoe-init CM candidates | #7 | built and wrapper-profile-verified (`#5` seed 20261021, 2-epoch profile); 45-unit matrix pending the Lead-run sweep | [`README`](06.7_cm_tahoe/README.md), `output/profile/profile_record.json` |
| `06.8_baselines/` five-fold baselines | #8 | built and verified — real five-fold record, `delta_mean` 0.1645 AUPRC@200 is the primary bar | [`README`](06.8_baselines/README.md), `output/summary_table.csv` |
| `06.9_evaluation/` 30-unit evaluation + report | #9 | planned | — |
| `06.10_capacity_ablation/` capacity ablation | #10 | planned (post-first-round) | — |
| `06.11_env_compound_ablation/` env-compound ablation | #11 | planned (post-first-round) | — |

## Boundaries

- The CM responses are general iPSC-CM responses at 24 h; nothing here becomes a
  hypoxia/reoxygenation efficacy claim, and no dose extrapolation is claimed.
- The checkpoint from `05_prnet_pretrain` is a pipeline bring-up artifact
  (raw-count controls vs log targets). It is **not** the formal pretrain and is
  not reused as one; formal pretrains come from `runs/`.
- A negative result — no pretraining gain, or no advantage over simple baselines
  — is a valid first-round outcome and is reported as such.
