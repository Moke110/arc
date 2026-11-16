# Single-dose tenv 5k model comparison

Implementation of approved [#226](https://github.com/Moke110/arc/issues/226).
The new step uses the released Bases G_5000 view, 1 µM treated wells and
matched controls. No normalization or axis is refitted here. Stage 3 is frozen.
The previous step stays accepted until the full successor matrix passes #231.

`execution/prepare.py --view <accepted view-5k> --candidate <empty directory>`
freezes bounded reader inputs, population, original-control mapping, Morgan2048,
and two shared 5-fold splits. Split seed 226 is registered before training;
validation units round up with `ceil(0.1 * outer_training_units)` (parents for
cold-drug, wells for warm-start). Warm wells are shuffled per parent and assigned
round-robin with a seeded starting fold. Validation retains at least one training
well per parent, explicitly failing if coverage cannot be met.

The public implementation is `output/api/tenv5k/protocol.py`: train-only Mean
and Mean_cell include each matched original control once; Mean_drug uses training
treated wells only and is unavailable for cold-drug. All methods predict treated
expression, with Δ obtained by subtracting the same matched control.

Tests run through the approved split, prediction and scoring boundaries:
`python -m unittest discover -s execution/tests -v`.

Shared learner interface: add `output/api` to the Python path and open
`tenv5k.inputs.Inputs(prepared_directory)`. It verifies every frozen digest before
returning memmapped `expression`, `control`, `delta`; treated rows are first and
indexed by `obs`. `rows(scenario, fold, partition)` returns shared integer rows;
`fingerprints(rows)` maps their normalized parent identity to float32 Morgan2048.
`genes` preserves the released gene_index order; `config` freezes seeds and stopping.
Predictors must output treated expression in the requested row/gene order.
`score_prediction(..., prediction_ids=...)` rejects reordered identities.

Run `execution/means.py --inputs <prepared directory> --candidate <empty directory>`
for all deterministic fold baselines. Each applicable run writes `prediction.npy`,
`identity.csv`, `per_well.csv`, `per_drug.csv`; `status.csv` explicitly marks the
five cold-drug Mean_drug runs unavailable. Payloads are candidates until #231's
complete matrix acceptance. Reconstruct inputs with prepare rather than copying
large accepted view payloads into a worktree.

The first verified freeze is recorded in `log/frozen-inputs/` with input/output
hashes and shared identities. The real 30-entry Mean slice completed with 25
applicable runs and five N/A entries. See `log/acceptance.md` for bounded validation
and resource evidence. This establishes shared inputs and baselines, not the final
90-run model comparison.

PRnet native slice (#228): `execution/prnet.py --inputs <prepared> --candidate
<empty directory> --source <immutable pert_model/PRnet/source> --scenario
cold-drug --fold 0 --seed 17 --device cpu` trains the registered Sci-Plex5k
architecture from random weights. Its sole chemical interface change is shared
binary Morgan2048 replacing upstream feature Morgan1024; the constant
`log10(1 + 1 µM)` multiplier remains explicit. GaussianNLL training and the native
noise-plus-output-Normal draw are retained. Batch512, Adam1e-3/weight_decay1e-8,
and ReduceLROnPlateau(.5, patience5, min_lr1e-8) follow Sci-Plex; validation uses
fixed seed22617, ordered complete validation, strict MSE improvement, patience50,
max1000 and best restoration. A singleton final training batch joins the previous
batch for BatchNorm. Source digests, adaptations, input identity, order and seeds
are frozen in each checkpoint configuration.

The CLI then launches `--replay-only` in a fresh process; this entry never trains.
Completion requires exact native sampled prediction replay on the recorded device
and Python/Torch/NumPy environment. Candidate outputs include checkpoint,
prediction, identity, per-well/per-drug scores, learning curve, configuration,
environment, resource/terminal records and hashes. Failed candidates have terminal
failure and no accepted publication. `--smoke-rows 8 --max-epochs 2` is an explicitly
marked bounded smoke, never a formal matrix run. Shared learner lifecycle lives
in `tenv5k.learner`; model adapters own native loss, initialization and sampling.
