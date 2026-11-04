# 0032 — PRnet formal transfer: paired Tahoe-pretrained vs random-init CM study

Status: accepted (user-approved 2026-11-02; the ticket graph is published as
GitHub issues #1–#11 in `Moke110/arc`)

**Update (2026-11-04, ADR 0033):** step labels in this record were rewritten by the
dataset-corpus merge; the decision itself is unchanged. "The frozen 02/04/05
outputs" now means: the `02_dataset_corpus` corpus tables, the `05_prnet_pretrain`
retired-view/pretrain evidence, and the `03.1` benchmark input — the retired v1 CM
views they name live read-only in `stage2/_frozen/cm-data-v1/`, and no result,
hash or digest key in this study changed.

## Context

`05_prnet_pretrain` proved PRnet can be trained on the Tahoe CM-schema view
(ADR 0031) but it selected on validation loss with raw-count controls against
log-scale targets, and it never finetuned on the CM corpus. The open question is
whether Tahoe pretraining improves prediction of **drug responses in human
iPSC-CM for drugs never measured there**. The full decision record — interview
questions Q1–Q24, the rejected options and the fixed technical defaults — is
`.scratch/stage2-prnet-formal-transfer/spec.md` (gitignored) and
`stage2/audit/prnet_training_design_review_2026-11-02.md`. This ADR records the
architecture-level decisions only; the concrete contract (columns, fold seed,
consensus membership) is frozen by the input-contract ticket in `06/input/`.

## Decision

1. **One new step, one frozen input.** `stage2/06_prnet_transfer/` builds a new
   two-source input (temp_env CM + Tahoe view) on the frozen 18,164-gene axis
   with explicit eligibility and QC layers. The frozen 02/04/05 outputs and the
   03.1 benchmark input are **not** regenerated or reused as the new folds.
2. **Drug-level five-fold split.** Folds partition `drug_key` (~6:2:2) so all
   doses and replicates of a drug stay together; every drug appears once as test
   and once as valid; both initialisations share the splits.
3. **Paired initialisations.** A random-init branch and a Tahoe-init branch share
   folds, seeds, learning-rate candidates and budgets; 3 seeds, 3 LRs, 90 CM
   candidates, 30 test units (folds × seeds × inits), no seed cherries.
4. **Fixed consensus input at evaluation.** CM training uses same-plate
   same-vehicle matched controls; validation, test and inference receive one
   consensus `c*` built only from permitted training-branch DMSO wells (per-plate
   log mean, then equal weight across plates). Experimental truth Δ uses matched
   controls and is computed only at scoring time. The consensus is never built
   from a test drug's response.
5. **Direction AUPRC@200 selects and scores.** Both stages train Gaussian NLL but
   select checkpoints and early-stop on validation direction AUPRC@200 (max
   mode); the test split is only scored after candidate selection and never
   informs LR, epoch, seed or capacity.
6. **Two acceptance questions, both reported.** (a) Does Tahoe pretraining beat
   random init? (b) Does either beat simple baselines (`no_change`,
   training-mean Δ)? A negative answer is a valid first-round outcome; no tuning
   until positive.
7. **Ablations follow, not gate.** Capacity (state latent 128, hidden 256) and
   extra-environment-compound training are separate post-first-round tickets.

## Consequences

- New pretrains are required: the 05 bring-up checkpoint is pipeline evidence
  only and is not the formal pretrain.
- Reusing 03.1's `cases.py`/`controls.npz` unchanged would leak per-plate test
  controls into the model input; this ADR rejects that path.
- `uns['control_mean']` keeps its raw-count meaning; any log-scale control cache
  gets a new name and explicit row membership.
- The study predicts general 24 h iPSC-CM responses on the fixed consensus
  background; it is not a hypoxia/reoxygenation efficacy claim and no
  dose-extrapolation claim is made.

## Rejected alternatives

- Matched test-plate controls as model input (leaks test-plate control
  information into validation/test inference).
- Selecting checkpoints on validation NLL or on test metrics.
- Screening all 525 corpus keys or the old 113/252 automatic tiers as if they
  were the pharmacology-qualified set.
- Modifying the frozen corpus, views, benchmark input, or `05` evidence in place.
