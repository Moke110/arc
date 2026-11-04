# 06_prnet_transfer — PRnet: Tahoe pretraining and temp_env CM finetuning

This step answers the question Stage 2 still has open after `05_prnet_pretrain`:
**does Tahoe-100M pretraining improve prediction of drug responses in human
iPSC-cardiomyocytes (iPSC-CM) for drugs never measured there?** The treatment
target is the general CM drug response; the held-out unit is the drug, not a
gene set. The study is pre-registered by an approved spec: the decisions, the
metric, the budgets and the negative-result clause were fixed before execution
(`.scratch/stage2-prnet-formal-transfer/spec.md`, gitignored; ticket graph in
`.scratch/stage2-prnet-formal-transfer/ticket-plan.md`).

Status: **first round executed and concluded (2026-11-02)** — audits, frozen input,
pipeline, three seeded pretrains, both 45-candidate matrices, baselines and the
30-unit evaluation are built, recorded and merged; the pre-registered result is
**negative** (Tahoe-init does not beat random-init and neither beats the simple
baselines). The two post-first-round ablations (#10/#11) are deferred with the
concrete remaining work recorded on their tickets. No claim below extends beyond
this scope.

**Post-merge path reconnection (2026-11-04).** The 02+04 merge (ADR 0033) moved the
CM v1 evidence to `stage2/_frozen/cm-data-v1/` and the Tahoe CM-schema view to
`stage2/05_prnet_pretrain/output/pretrain_view/`. `06.2`, `06.3` and `06.4` data
roots were reconnected to those locations and to
`stage2/02_dataset_corpus/output/corpus/v2/tahoe.h5ad`: paths were reconnected to
the post-merge locations; the recorded results, hashes and digest keys are
unchanged.

## The first-round design (fixed)

| | |
|---|---|
| question | prediction of the transcriptomic response of unseen drugs in iPSC-CM; paired Tahoe-init vs random-init comparison |
| scope | 24 h only; Tahoe + temp_env only (no ScreenSeq, no DToxS branch, no GSE247541) |
| CM data | temp_env (`GSE262419`) restricted to the eligibility ∩ QC set; low-response and toxicity/stress wells are kept and stratified |
| gene axis | the frozen 18,164 genes of the CM `corpus/v1` tables (now `stage2/_frozen/cm-data-v1/corpus/v1`), order preserved |
| model | the author PRnet `PGM`: hidden 128, drug latent 64, state latent 64, noise 10; FCFP4 × log10(dose+1); full-parameter finetune; Gaussian NLL loss |
| split | five folds by `drug_key` (~6:2:2); every dose/replicate of a drug stays together; test and valid once per drug |
| control | CM train uses same-plate/same-vehicle matched controls; validation, test and inference get one fixed training CM consensus `c*`; experimental truth Δ uses matched controls and never enters a model input |
| primary metric | direction AUPRC@200, macro-averaged over drug×dose conditions after within-condition replicate aggregation |
| selection | validation AUPRC@200 (max) drives checkpointing, early stopping and schedule; test never touches selection |
| repetition | 3 seeds; 3 independent Tahoe pretrains reused across 5 folds; 2 inits × 5 folds × 3 seeds × 3 LRs = 90 CM candidates; 30 CM test units |
| baselines | `no_change`, training-mean Δ, plus dose-mean Δ / ridge-fingerprint as auxiliary |

## Layout

Per ADR 0030 §3 this container splits into substeps, one per published ticket.
**Nine are built** — `06.1_eligibility/`, `06.2_qc/`, `06.3_input/`,
`06.4_pipeline/`, `06.5_tahoe_pretrain/`, `06.6_cm_random/`, `06.7_cm_tahoe/`,
`06.8_baselines/`, `06.9_evaluation/` — each carrying `AGENTS.md`, `README.md`,
`execution/` and `output/`. **Two are deferred ablations and have no directory**
(the post-first-round work was not needed once the round concluded negative):
`06.10_capacity_ablation/` (#10) and `06.11_env_compound_ablation/` (#11). See
[`AGENTS.md`](AGENTS.md) for what each owns and [`docs/adr/`](../../docs/adr/)
for the contracts. There is exactly one eligibility audit, one QC audit, one
frozen input build, and one evaluator; every experiment consumes them rather than
re-deriving a split.

## Status table

| substep | ticket | status | evidence |
|---|---|---|---|
| `06.1_eligibility/` eligibility audit | #1 | built and verified — rule v1.3, handoff v2 (118 include / 28 uncertain / 314 exclude) | [`README`](06.1_eligibility/README.md), `output/current.json` |
| `06.2_qc/` per-observation QC | #2 | built and verified — locked rule set | [`README`](06.2_qc/README.md) |
| `06.3_input/` frozen inputs, folds, consensus | #3 | frozen and verified — 118 drugs / 951 wells, 125/125 verify | [`README`](06.3_input/README.md), `output/manifest.json` |
| `06.4_pipeline/` pipeline + smoke | #4 | Phase A built; Phase B real-data smoke green (17/17 steps, GPU peak 2.9 GB, host peak 8.5 GiB) | [`README`](06.4_pipeline/README.md), `output/smoke/smoke_record.json` |
| `06.5_tahoe_pretrain/` 3 Tahoe pretrains | #5 | executed — 3/3 seeds (29/95/53 epochs, best valid AUPRC@200 0.2008/0.2222/0.2313) | [`README`](06.5_tahoe_pretrain/README.md), `output/formal/tahoe_pretrain_record.json` |
| `06.6_cm_random/` 45 random-init CM candidates | #6 | executed — 45/45 (40 ok / 5 capped), LR 1e-3 in all folds, 15 checkpoints selected | [`README`](06.6_cm_random/README.md), `output/record/record.json` |
| `06.7_cm_tahoe/` 45 Tahoe-init CM candidates | #7 | executed — 45/45 (42 ok / 3 capped), per-fold LR 1e-4/1e-4/1e-3/1e-3/1e-4, 15 checkpoints selected, warm start 14/14 keys per seed | [`README`](06.7_cm_tahoe/README.md), `output/record/record.json` |
| `06.8_baselines/` five-fold baselines | #8 | built and verified — real five-fold record, `delta_mean` 0.1645 AUPRC@200 is the primary bar | [`README`](06.8_baselines/README.md), `output/summary_table.csv` |
| `06.9_evaluation/` 30-unit evaluation + report | #9 | executed — 30/30 test units scored, 0 failed; **negative first round**: Tahoe-init −0.0063 vs random-init (95% CI [−0.0074, −0.0052]) and both far below `delta_mean` 0.1645 | [`README`](06.9_evaluation/README.md), `output/record/summary.json` |
| `06.10_capacity_ablation/` capacity ablation | #10 | deferred (post-first-round, not a gate) — plan and remaining work recorded on the ticket | — |
| `06.11_env_compound_ablation/` env-compound ablation | #11 | deferred (post-first-round, not a gate) — plan and remaining work recorded on the ticket | — |

## Boundaries

- The CM responses are general iPSC-CM responses at 24 h; nothing here becomes a
  hypoxia/reoxygenation efficacy claim, and no dose extrapolation is claimed.
- **Record provenance.** Locked audit records under `06.1`, `06.4` (smoke cards),
  `06.5`, `06.6` and `06.7` embed the absolute worktree path where the run was
  produced (e.g. `/home/chang/arc-wt/t05-pretrain`). ADR 0030 §4 keeps frozen
  evidence path strings byte-for-byte, so these historical strings are deliberately
  not rewritten; the 06.9 record, regenerated after the code review, resolves every
  path repo-relatively.
- The checkpoint from `05_prnet_pretrain` is a pipeline bring-up artifact
  (raw-count controls vs log targets). It is **not** the formal pretrain and is
  not reused as one; formal pretrains come from `runs/`.
- A negative result — no pretraining gain, or no advantage over simple baselines
  — is a valid first-round outcome and is reported as such.
