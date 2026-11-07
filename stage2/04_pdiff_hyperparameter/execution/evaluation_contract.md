# Held-out evaluation artifacts

Run the evaluator with the `arc` interpreter. Repeating the command against
the same frozen inputs and a fresh output directory replays the evaluation.
Use the scheduling rules in the step README for full-size work.

```bash
python stage2/04_pdiff_hyperparameter/execution/script/evaluate.py \
  <plan.json> <outcomes.json> <fresh-output-dir>
```

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

The evaluator verifies the checkpoint bytes against the training record,
matches the selected recipe and exact case/group order, and checks the
selected step against the strict minimum in the validation history.
Paths relative to the plan, outcomes, and source manifest stay under their
own directories. The frozen predictor request retains its absolute
`input_root`.

## Read the outputs

`per_group.csv` scores the means of the matching wells, with undefined
metric reasons. `per_run.csv` includes all expected slots, equal-weight
group means, valid-group counts and status. `summary.csv` gives run means,
sample standard deviations, valid-run counts and completion counts for
each dataset and model. Undefined metric cells are blank.

Each reference directory binds one frozen group list to two single Δ vectors.
The vectors apply to every listed group in that fold. `delta_mean` gives
equal weight to the mean Δ of each training group. `no_change` is zero Δ.
Reference rows have a blank seed and contribute one result per fold.

`evaluation.json` binds the plan, outcomes and output hashes, records
whether every expected result is complete, and discloses the prior Step 03
scoring of the formal tenv test memberships.
