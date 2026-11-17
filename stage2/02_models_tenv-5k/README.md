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

TranSiGen uses the ordinary upstream `train_TranSiGen.py` configuration:
5000→1200 encoder, latent100, decoder800, Morgan2048→400 embedding,
beta0.1, dropout0.1, batch64, Adam lr0.001/weight_decay0.00001, no scheduler.
The 5000-gene dimension, shared Morgan fingerprints, random initialization and
shared split/stopping replace script inputs/defaults; native selective Xavier
initialization and summed VAE objective remain unchanged. The approved
`transigen-native-f64-v1` numerical variant initializes in float32, then uses
float64 parameters, floating buffers, inputs, forward/loss/backward, Adam moments,
validation and inference. Original float32 CPU latent normal draws are cast to
float64 before scale/add; their call order stays native. Predictions serialize
as float64 with a finite check after conversion. The author cuDNN policy remains
part of the execution identity. Replay rejects old or changed precision identities.
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

MLP_Morgan is the XPert-defined chemical baseline with this project's
preregistered configuration: Morgan2048→Linear128→ReLU→Linear5000, a linear
output and no control, dose, dropout or pretraining input. It uses mean MSE,
Adam lr0.001/weight_decay0.00001, batch64 (including singleton batches), no
scheduler, and the shared strict validation-MSE stopping/best-restore lifecycle.
Complete baseline training code was not found in the official GitHub/Figshare;
this configuration is not claimed as a reproduction of author training code.

`execution/mlp.py --inputs <freeze> --candidate <empty directory> --scenario
cold-drug --fold 0 --seed 17` trains then invokes exact checkpoint replay in a
fresh process. `--smoke-rows 8 --max-epochs 2` is explicitly a smoke candidate;
`--replay-only` verifies without training. The candidate records the project
implementation digest, configuration, identities, predictions, matched-control
Δ metrics, learning curve, environment and measured process peak RSS.

Matrix orchestration (#231) uses `execution/run_matrix.py --plan-only` without
loading models or input payloads. `--all` requires explicit `--inputs`, `--means`,
`--prnet-source`, `--transigen-source`, `--candidate`, and `--device` paths/identity.
`--resume` skips only complete, exact-replayed formal identities after checking
all artifact hashes and the frozen plan; changed code/configuration is rejected.
`--limit-runs 1` profiles the first full registered run; it does not shorten its
1,000-epoch cap or patience 50. Continue it with the same arguments and `--resume`.
`--accept-only` evaluates the complete matrix without training, and `--publish`
requests publication only after complete OOF metric and CSV figure replay.

Comparison first reconstructs every well's out-of-fold prediction across all five
folds for each seed, then averages repeated-well Δ vectors within each parent,
then computes metrics, then averages the three seed metrics. This matters in
warm-start, where one parent spans several test folds. Paired parent bootstrap
uses 10,000 common draws with frozen statistical seed 2262026. Cold-drug Δ PCC
alone selects a winner; otherwise both candidates remain. Figures read only the
three recorded CSV plot-value tables and the pure shared palette implementation.

The historical native-float32 90-run attempt stopped at 18 complete runs: the
19th registered identity has a reproducible native float32 non-finite loss.
See [the terminal-block record](log/231-terminal-block.md) for diagnosis and
preserved evidence. The author CUDA backend correction does not resolve that
failure. Previous accepted results and the old step remain available; the new
comparison has not been accepted or published.

On 2026-11-17 the user approved the numerical-method amendment in ADR 0074.
All 30 TranSiGen identities must run from zero as `transigen-native-f64-v1`;
the three historical float32 successes remain evidence only. Exactly 15 completed
cold-drug PRnet entries may be imported only after immutable original
runtime551a3c9 identity receipts, fresh-process exact sampled replay and scientific
equivalence checks. The other 45 learner entries require training; deterministic
Mean reuse requires manifest/scoring compatibility. The complete 90-run acceptance
remains outstanding. No result is inferred from bounded variant smoke checks.

Interrupted same-identity `--resume` now preserves the entire prior run under
candidate `attempts/` before a clean restart. A trained run awaiting replay uses
only the replay entry point, with pre-replay terminal evidence retained; a
complete run with corrupt digests still fails closed. Failed replay recovery is
recorded and can subsequently restart the same registered identity.
Publication preflights every destination, stages digest-verified copies, and
journals installation renames with rollback. A terminated installation is rolled
back from its journal before another publication attempt. The heavy-command
boundary refuses to start below 12 GiB available memory. Native imports and their
sibling imports disable bytecode additions in immutable source trees.

Tracked metric CSVs include all six methods, explicitly including cold-drug
Mean_drug as N/A. Truth-constant parent identities, undefined denominators and
constant prediction counts are recorded with plot values and displayed as figure
annotations. Prediction-constant counts count seed-parent vectors; deterministic
Mean methods have a single vector per parent. True constants retain error scores.

The v2 matrix plan separates shared input/evaluation identity from each method's
scientific identity and implementation hashes. CSV report rows include
`method_identity` (paired tables include candidate/reference identities), so the
stable `TranSiGen` method name is explicitly identified as
`transigen-native-f64-v1`. Acceptance records all 90 execution origins, including
imported versus newly trained evidence; it still refuses an incomplete matrix.

Use `execution/run_matrix.py --import-cold-prnet --inputs <freeze> --means <means>
--prnet-source <source> --transigen-source <source> --candidate <new candidate>
--device cuda --import-original <original matrix-231> --import-runtime
<preserved runtime551a3c9 worktree> --compatibility-audit <reviewed JSON>` before
continuing `--all --resume` with the same paths/code. Import is restricted to the
15 registered cold PRnet identities and the approved original plan digest. The
audit binds `source_plan_sha256`, exact `target_method` from the new plan,
`source_runtime: "551a3c9"`, `reviewed_by`, and `changed_code` rationales keyed by
changed method dependency paths. Only the explicitly reviewed original→new shared learner, PRnet entrypoint and
new attempt-identity helper digest transitions are allowed; PRnet adapter, inputs and scoring must be unchanged.
A rationale alone cannot authorize other code or scientific changes.

Every import runs a read-only fresh process using the frozen runtime, original
source, recorded device and exact software environment. It checks checkpoint
configuration and exact prediction dtype/values. Original run files, including
configuration, terminal, checkpoint and manifest, are copied byte for byte into
the new candidate, accompanied by `import-receipt.json` and `import-provenance/`
(original plan, reviewed audit, frozen Python code, timed replay evidence).
Original artifacts are never rewritten or relabelled as fresh training. These
portable receipts are rechecked on resume, acceptance and publication. New runs
must have a preregistered execution identity and a separate execution receipt.

Each new training attempt uses an OS-generated nonce, separate from model RNG.
The registration's canonical JSON SHA256 and nonce are passed to the training and
fresh replay child via `--registration`, persisted in configuration, checkpoint
and manifest, and checked on record/resume/acceptance. A historical complete run
cannot become fresh evidence even in a previously registered slot. Retries use a
new nonce and retain the former registration and interrupted artifacts.

Runtime identity scans owned direct `execution/*.py` entrypoints/helpers and
recursive `output/api` Python source. Nested execution candidates/checkpoints and
their archived runtime provenance are products, so imports cannot change their
own frozen code identity. A candidate may use another execution subdirectory;
source-root ancestors, direct execution source files and locations within the
public API source tree are rejected to keep runtime code in the identity.
