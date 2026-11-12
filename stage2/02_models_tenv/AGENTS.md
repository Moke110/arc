# AGENTS.md — stage2/02_models_tenv

Routing for the Stage 2 tenv model-comparison step (#163, #164, #165). ADR 0063
is the authority for the metric families, run matrix, palette and figures; the
step contract is [`../AGENTS.md`](../AGENTS.md); the shared readers, evaluator,
palette and drawing API are owned by [`../01_bases/`](../01_bases/README.md) and
are never re-implemented here. Human-readable results are in
[`README.md`](README.md).

## The accepted matrix and its one command (#165)

`execution/run_matrix.py` is the one-command publish/replay. It builds each
fold's complete predictor-input release, trains PRnet for the declared
2 folds x 2 seeds **serialized** through `execution/capped.sh` under
`/usr/bin/time -v`, writes the deterministic `Mean` reference, scores the whole
run set with the shared #164 driver, verifies the #162 plot-value contract,
publishes the stable products, and renders and replay-verifies both figures.

```bash
python execution/run_matrix.py --plan-only   # print the declared identities
python execution/run_matrix.py --all         # candidate -> acceptance -> publish
python execution/run_matrix.py --all --resume
```

Candidate state lives in `execution/checkpoints/matrix-165/`, never in
`output/`; a run that is not `complete` aborts before any stable path is
touched. `figure/script/render_figures.py --figure-data figure/data --out-dir
figure` rebuilds both PNGs from the tracked plot values alone and never opens
the step's `output/` tree.

## Environments (two, recorded)

- **Training/prediction:** the pinned isolated interpreter
  `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python`. The shared
  Python 3.13 environment cannot import the upstream trainer package. The
  lockfile is `log/env.prnet.lock`; every run re-records the interpreter,
  torch/numpy versions, device and the imported upstream module digest in its
  `identity.json`.
- **Input preparation:** `/home/chang/miniconda3/envs/arc/bin/python`, because
  the accepted `bases` readers need that environment's newer `anndata`.
  `log/env.prepare.json` records the preparation interpreter and package
  versions.
- The pinned environment needs `NUMBA_CACHE_DIR` (numba 0.60 cannot write its
  cache beside the read-only package); `prnet_fit_predict.py` sets it, and every
  environment record includes it.

## Inputs (read, never re-derived)

- Training corpus `../01_bases/output/training-corpora/tenv.h5ad`
  (4,206 treated wells × 18,164 genes; `layers/delta` is the target,
  `layers/control` the same-plate baseline, `X` the treated expression).
  Read metadata through `bases.corpora.training_header` and expression only in
  bounded `h5py` row slabs (`float32`, 64-row slabs). Never build a dense
  wells × genes matrix and never open the corpus with
  `anndata.read_h5ad(backed="r")` (it materialises every layer).
- Fold/seed identities `../01_bases/output/fold_n_seed/` through
  `bases.fold_seed.FoldSeedIndex.load`; the step declares 2 folds × 2 seeds and
  selects `folds.csv`/`seeds.csv` rows by prefix, never by name.
- Morgan-2048 comes from the corpus's normalized `parent_smiles` with the
  retired adapter input product's representation: RDKit Morgan, radius 2,
  2048 bits, chirality on, binary, keyed by `drug_key` and ordered by
  `drug_key`. Do not invent another chemical representation.
- `pert_model/PRnet/source` is immutable: only `models.PRnet.PGM` is imported,
  and the source tree is never edited or copied into this step.

## Commands (run from the step directory)

```bash
# Bounded predictor inputs (shared environment; corpus/fold-product paths must
# be passed explicitly from a linked worktree, whose ignored payloads are absent)
execution/capped.sh /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python \
  execution/script/prepare_inputs.py --corpus <tenv.h5ad> \
  --fold-product <fold_n_seed> --output output/bounded/inputs \
  --fold-row 0 --seed-row 0 --train-wells 512 --validation-wells 128 --test-wells 64

# Bounded fit/predict (pinned environment)
execution/capped.sh /usr/bin/time -v \
  /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python \
  execution/script/prnet_fit_predict.py --request <request.json> \
  --config execution/models/prnet/bounded_run.json
```

The bounded run publishes `prediction.npy` (all 18,164 genes, finite, declared
`(n_cases, 18,164)` shape), case identities, the selected checkpoint, the
effective config, `identity.json` (hashes), `checkpoint.json` (history,
early-stopping record, behavior probes, resource peak) and `terminal.json`.
A failed, interrupted or over-budget run writes `terminal.json` with that
status and no prediction. The bounded run is #163 admission only; the accepted
2 × 2 matrix is #165 and is recorded in `log/165-acceptance.json`.

## Early stopping and selection (ADR 0063)

`execution/script/early_stopping.py` owns both rules: stop after 5 consecutive
epochs without a strict validation-loss decrease, independently stop after
20 consecutive epochs above the running minimum, never exceed 100 epochs, keep
the lowest-validation-loss checkpoint. The 5-epoch rule dominates the 20-epoch
rule (recorded in `checkpoint.json` and `terminal.json` evidence); the
20-epoch rule is still implemented and checked on its own counter.

## Accepted products (ADR 0062, #165)

| path | content | tracked? |
|---|---|---|
| `output/inputs/{fold_X}/{seed_Y}/` | complete predictor-input release per matrix cell and its manifest | ignored payload |
| `output/runs/{entry}/{fold_X}/{seed_Y}/` | `prediction.npy`, `case_ids.json`, `best.pt`, `config.json`, `checkpoint.json`, `identity.json`, `terminal.json`, `progress.json`, `resources.json` | ignored payload |
| `output/evaluation/` | `per_well.csv`, `per_group.csv`, `per_run.csv`, `evaluation.json` from the shared driver | ignored payload |
| `output/figure_data/` | candidate copies of the two plot-value tables (identical to `figure/data/`) | ignored payload |
| `figure/metrics.png`, `figure/learning_curve.png` | the report figures | tracked |
| `figure/data/metrics_values.csv`, `figure/data/learning_curve_values.csv` | complete plotted values in the shared `PLOT_VALUE_COLUMNS` schema | tracked |
| `figure/script/render_figures.py` | CSV-only renderer (no matrices, checkpoints or truth) | tracked |
| `log/165-acceptance.json`, `log/resources/165-*` | acceptance record with hashes/statuses/peaks and the raw `/usr/bin/time -v` files | tracked |

## Tests

Run each file individually with the repository Python 3.13; never the full
suite beside another heavy job:

```bash
/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python -m unittest discover \
  -s execution/tests -p test_early_stopping.py -v
/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python -m unittest discover \
  -s execution/tests -p test_prnet.py -v
/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python -m unittest discover \
  -s execution/tests -p test_inputs.py -v
```

`test_prnet.py` is the ADR 0063 seam-3 test: it runs the real command on a
synthetic release and asserts reload-exact prediction, case identity/order,
finite full-axis output, the early-stopping record, and honest terminal status
for failed and over-budget runs.

Seam-2/4 tests added by #165 run with `/home/chang/miniconda3/bin/python`:
`test_run_matrix.py` (declared 2 x 2 plan, `/usr/bin/time -v` parsing,
fail-closed terminal status and candidate-only publish, recorded figure layout),
`test_mean_reference.py` (the retired `delta_mean` arithmetic, cold-drug
semantics, driver-validated artifacts) and `test_figures.py` (both PNGs rebuild
byte-identically from tracked `figure/data` with the step `output/` unreadable).

## Records and resources

`log/` holds the environment locks, the bounded-run acceptance record
(`log/163-bounded-run.json`), raw `/usr/bin/time -v` measurements
(`log/resources/`) and the per-phase `progress.json` beside the run's published
artifacts. Disposable run payload
lives under `output/` and stays untracked. Heavy commands run one at a time
through `execution/capped.sh` (systemd user scope `MemoryMax=10G
MemorySwapMax=8G`, or an equivalent 10 GiB `RLIMIT_AS`) with 2 threads and
`/usr/bin/time -v`; record the measured peak.
