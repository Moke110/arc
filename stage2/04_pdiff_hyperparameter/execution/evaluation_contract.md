# Held-out evaluation artifacts

The active retrieval contract is [ADR 0059](../../../docs/adr/0059-stage2-significance-retrieval.md)
and the [shared evaluator](../../evaluation/README.md). Only tenv retrieval is
reported, using `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` on full-CM
measured significance. Run under the step's scheduling and memory limits:

```bash
python stage2/evaluation/publish.py --step 4
python stage2/evaluation/evaluate.py replay --step 4
```

The following producer-plan and prediction schemas describe frozen training
provenance. They do not authorize Tahoe PR/AUPRC or the retired positional
scoring interface.

## Freeze the plan before test inference

Use schema `step04-evaluation-plan-v1` and set `stage` to `formal`. The
`folds` array lists the three folds of each of `tahoe` and `tenv`. Each entry
contains `request` (path to the frozen Step 04 predictor request relative to
the plan directory), `request_sha256`, and `recipe_sha256_by_seed` (a mapping
from strings `"17"`, `"42"`, and `"2026"` to the selected, normalized
recipe hashes that will appear in `run.json`). `stage: "fixture"` supports
bounded synthetic inventories. The request binds the source release,
predictor files, fold membership and shared gene axis. Both release manifests
pin training Δ, the ordered held-out cases, and evaluator-only group truth
and directional labels.

For an independent formal dataset evaluation, add `dataset: "tenv"` or
`dataset: "tahoe"`, keep exactly that dataset's three folds and all three
seeds. The default formal plan covers both datasets. Every formal evaluation
requires the externally recorded plan digest in `--plan-sha256`. Plans written
by `formal.py freeze` also bind the selected recipe files, inference batching
and runtime identities. The worker requires that same plan hash before
accessing predictor inputs.

## Record the outcome of each planned slot

Use schema `step04-evaluation-outcomes-v1` with an `outcomes` array. Each
entry has `dataset`, integer `fold`, integer `seed`, and `status`. For a
completed slot include `prediction` (path to its prediction manifest
relative to the outcomes directory) and `prediction_manifest_sha256`. For
a failed, timed out, incomplete, interrupted or cancelled slot, use
`failed`, `timeout`, `incomplete`, `interrupted`, or `cancelled` and
include `reason`. An omitted slot appears in the output as `missing`.

A completed slot's directory contains the following files.

- `run.json` and `best.pt` from Step 04 training. The producer verifies the
  saved EMA state through the checkpoint reload path before predicting.
- `prediction.npy` with float32 Δ, one row per held-out case in frozen
  `test_cases.csv` order and one column per gene on the frozen axis.
- `cases.json` with the ordered array of `[case_id, group_id]` pairs used
  for those rows.
- `prediction.json` with schema `step04-heldout-prediction-v1`.
  Copy `identity`, `recipe_sha256`, `checkpoint_sha256`, `best_step`,
  `best_validation_loss`, and `model_state_sha256` from `run.json`.
  Include `run_sha256`, `cases_sha256`, and `prediction_sha256` for the
  respective files.
- For every formal plan, `prediction.json` also includes `plan_sha256`
  and `replay_sha256`. The replay file `replay.json` binds the exact
  validation loss and probe result to the same plan and EMA state. The
  outcomes inventory carries that plan hash too. The evaluator verifies
  these bindings before opening any held-out truth or labels.

The evaluator verifies the checkpoint bytes against the training record,
matches the selected recipe and exact case/group order, and checks the
selected step against the strict minimum in the validation history.
Paths relative to the plan, outcomes, and source manifest stay under their
own directories. The frozen predictor request retains its absolute
`input_root`.

## Read the outputs

`output/retrieval/` holds the generated run/group/dose tables, summary, pooled
PR points, PR curve and AUPRC PNG/SVG figures, and hash manifest. Exact group AP
and fixed-group pooled PR follow the shared contract. Undefined scores retain
counts and reasons. Continuous metrics remain available. Deterministic
reference seed copies are retained as outcomes and counted once per fold in
uncertainty summaries. Frozen producer bindings and measured label hashes are
checked again before final publication.
