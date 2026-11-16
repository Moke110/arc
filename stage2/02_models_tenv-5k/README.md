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

TranSiGen uses the ordinary upstream `train_TranSiGen.py` configuration:
5000→1200 encoder, latent100, decoder800, Morgan2048→400 embedding,
beta0.1, dropout0.1, batch64, Adam lr0.001/weight_decay0.00001, no scheduler.
The 5000-gene dimension, shared Morgan fingerprints, random initialization and
shared split/stopping replace script inputs/defaults; native selective Xavier
initialization, summed VAE objective and CPU latent sampling remain unchanged.
Singleton final training batches are skipped as upstream; evaluation retains them.
Inference directly calls native `forward(x1, features)[3]`, so treated labels and
the target reconstruction branch never consume prediction RNG. Fixed seed,
row order and batch64 apply identically to validation, prediction and reload.
Source code is hash-checked before construction; no pretrained asset is loaded.

`execution/transigen.py --inputs <freeze> --source <TranSiGen/source>
--candidate <empty directory> --scenario cold-drug --fold 0 --seed 17`
builds a candidate and launches a fresh Python process for exact sampled replay.
Use `--smoke-rows 8 --max-epochs 2` for bounded integration checks; smoke status
is explicit and cannot substitute for a formal matrix entry. `--replay-only`
performs no training. Native-source tests accept `TRANSIGEN_SOURCE` when the
immutable payload is materialized in another checkout.
