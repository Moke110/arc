# 06_prnet_transfer — PRnet: Tahoe pretraining and temp_env CM finetuning

This step answers the question Stage 2 still has open after `05_prnet_pretrain`:
**does Tahoe-100M pretraining improve prediction of drug responses in human
iPSC-cardiomyocytes (iPSC-CM) for drugs never measured there?** The treatment
target is the general CM drug response; the held-out unit is the drug, not a
gene set. The study is pre-registered by an approved spec: the decisions, the
metric, the budgets and the negative-result clause were fixed before execution
(`.scratch/stage2-prnet-formal-transfer/spec.md`, gitignored; ticket graph in
`.scratch/stage2-prnet-formal-transfer/ticket-plan.md`).

Status: **planned** (2026-11-02). No audit, no frozen input and no new training
has been executed. Each section below flips to executed only when its ticket
closes with artifacts and evidence.

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

`eligibility/` · `qc/` · `input/` · `pipeline/` · `runs/` · `baselines/` ·
`evaluation/` · `ablations/capacity/` · `ablations/env_compounds/` — see
[`AGENTS.md`](AGENTS.md) for what each owns and
[`docs/adr/`](../../docs/adr/) for the contracts. There is exactly one
eligibility audit, one QC audit, one frozen input build, and one evaluator;
every experiment consumes them rather than re-deriving a split.

## Status table

| area | ticket | status | evidence |
|---|---|---|---|
| eligibility audit | #1 | planned | — |
| QC audit | #2 | planned | — |
| frozen inputs / folds / consensus | #3 | planned | — |
| pipeline + smoke | #4 | planned | — |
| 3 Tahoe pretrains | #5 | planned | — |
| 45 random-init CM candidates | #6 | planned | — |
| 45 Tahoe-init CM candidates | #7 | planned | — |
| baselines | #8 | planned | — |
| evaluation + report | #9 | planned | — |
| capacity ablation | #10 | planned (post-first-round) | — |
| env-compound ablation | #11 | planned (post-first-round) | — |

## Boundaries

- The CM responses are general iPSC-CM responses at 24 h; nothing here becomes a
  hypoxia/reoxygenation efficacy claim, and no dose extrapolation is claimed.
- The checkpoint from `05_prnet_pretrain` is a pipeline bring-up artifact
  (raw-count controls vs log targets). It is **not** the formal pretrain and is
  not reused as one; formal pretrains come from `runs/`.
- A negative result — no pretraining gain, or no advantage over simple baselines
  — is a valid first-round outcome and is reported as such.
