# Stage 2 whole-transcriptome benchmark contract

**Status:** accepted (2026-10-21)

Stage 2 must choose one open-source model to supply whole-transcriptome
perturbomes to future MIRI drug scoring. Setup (ADR 0021, `01_setup/README.md`)
proved that nine models run in isolated environments, but it did not compare
them: the authors' benchmarks use different gene spaces, doses, cell contexts and
splits, so cross-paper numbers are not evidence. `stage2/03_models_cmbench/`
therefore fixes a single cold-drug benchmark **before** any model is scored.
The contract below is frozen; the executed leaderboard and report are derived
from it.

## Decisions

### 1. Data view and target

- **Dataset:** `stage2/02_cm_data/output/datasets-1/temp_env_all-smiles.h5ad`, the
  largest curated view in the ADR-0026 bundle: 4,478 wells (4,206 perturbed /
  272 control), 460 drug keys, 18,164 genes, 16 plates, three measured doses
  (0.2 / 1 / 10 µM), one 24 h timepoint, `ipsc_cm`. The `atc-moa` views are
  smaller and would make the comparison depend on a pharmacology subset; the
  `all-smiles` view is the broadest reproducible panel.
- **Δ space:** `Δ = log1p_cp10k(treated) − mean(log1p_cp10k(same-plate DMSO))`.
  `cp10k` is counts per 10k; `log1p` is the frozen convenience transform. The
  subtraction is always by the well's own `control_group` (physical plate), never
  a global control mean.
- **Group:** one `(drug_key, dose)` pair; 1,380 groups. The target for a group is
  the mean Δ over its 3–6 wells. The model emits per-well Δ; the runner averages
  within `group_id`. This keeps replicate noise visible to the adapter and makes
  the scored object the group mean.
- **Drug identity:** the desalted parent InChIKey with stereochemistry retained
  (ADR 0026), with one representative parent SMILES per drug key. Salts and
  free bases are the same drug and cannot straddle a cold-drug split.

### 2. Cold-drug split

- Five folds, **6:2:2** per fold: 276 train / 92 valid / 92 test drugs (ticket
  06). Every drug is test exactly once and valid exactly once; `valid(k) =
  test(k−1)`.
- Assignment is by **connectivity-block cluster** (458 blocks) with a
  well-count-stratified capacity-constrained greedy, seed `20261021`, so
  structural analogues cannot leak between roles.
- The long `splits.parquet` has one row per treated well (21,030) and is
  fail-closed on ten gates: 276/92/92 drugs, disjoint roles, no test drug in
  train or valid, test wells 837/846/837/846/840, and complete coverage.
- Controls are never held out and never enter the test bundle: they are shared
  once via `output/input/controls.npz` (16 groups × 12 cells) and the adapter receives
  only the control inputs it needs.

### 3. Entrant set and fidelity labels

The admission gate is **cold-drug capability**: the model must consume a
structure-derived drug representation so an unseen drug can be scored. Seven
models pass:

| model | cold-drug path | `ideal_fit` |
|---|---|---|
| PRnet | rFCFP4 × log-dose | full |
| CPA | frozen Morgan-2048 pert encoder + unseen drug dose scaler | full |
| chemCPA | RDKit2D molecule encoder | partial (978 pretrain omitted) |
| biolord | ordered RDKit2D + dose, transductive | partial |
| PerturbNet | ChemicalVAE → cINN | partial (no control transcriptome) |
| Squidiff | rFCFP4 × log-dose | full |
| CMonge | RDKit2D + dose → Monge map | full |

- **scVIDR and DECANT are excluded from Stage 2 entirely.** scVIDR has no drug
  input (its prediction is latent arithmetic between observed control and treated
  cells of the same drug), so an unseen drug has no rows to interpolate. DECANT's
  forward pass also has no drug input and its target is the mean of the same
  treated bag it is handed; with 3 wells per group it sees its own answer. A
  cold-drug version of either is a new architecture, not an adapter.
- There is **one leaderboard**. `ideal_fit` is a report column (`full` /
  `partial` / `none` for baselines), never a second table. The partial labels are
  not excuses: every deviation is also listed per model in `FIDELITY.md`.

### 4. Training rules

- **Random init only.** No released checkpoint is loaded anywhere. The
  author-mandated multi-stage recipes are reproduced on each fold's training
  drugs/cells; where an author stage needs external data, the substitution or
  omission is declared (PerturbNet's ZINC ChemicalVAE → train-fold SMILES;
  chemCPA's L1000 978 pretrain → omitted).
- **Author defaults** for architecture, optimizer, schedule, loss and
  preprocessing; validation selects the epoch/checkpoint only. No per-model HPO.
  Fixed-step models evaluate a periodic checkpoint grid (Squidiff 10% of 10^5;
  biolord every 50 of 500; CMonge every 1,000 of 10,000); early-stopping models
  use the author's patience.
- **18,164 genes for every model**, in corpus order, with no HVG selection, gene
  projection, or zero-fill. This is a config/data-driven dimension, not an
  architecture change. Models whose author space was smaller carry that
  extrapolation in the ledger.
- **Budget:** 15.5 h soft / 16 h hard per model across five folds; 6 h per fold.
  Folds run sequentially on one GPU. A timed-out fit predicts from its last
  checkpoint and is flagged; a crash is retried once only for environmental
  failures; a fold with no usable checkpoint is `failed`, never imputed.
- **Single seed** `20261021`; deterministic flags where the framework allows.
  No second seed: it would double the model budget.

### 5. Metric

- **Primary:** direction-aware AUPRC at **k = 200**. For each group, the true
  positive set is the top-200 genes by true mean Δ (up) plus the bottom-200
  (down); the prediction is ranked the same way. AUPRC is computed per direction
  and macro-averaged, then macro-averaged over groups, then averaged over folds.
  A random predictor scores 1.1% per direction.
- **Secondary:** PCC, Spearman, R² about the true mean, raw RMSE in Δ units, and
  NRMSE; k = 100 and 500 are reported as sensitivity, never as a second ranking.
- **No micro AUPRC.** Pooling groups would let large groups dominate; the group
  macro is the contract.
- **Uncertainty:** stratified 92-drug cluster bootstrap, 2,000 resamples, seed
  `20261022`, paired adjacent differences with percentile 95% CI. A model with
  fewer than five completed folds is reported as partial and cannot be named
  winner.
- **Sanity anchors:** constant predictions get NaN PCC/Spearman (not 0);
  `no_change` reproduces the 1.1% prevalence; perfect and sign-flipped
  predictions are checked in the offline tests.

### 6. Adapter and artifact contract

- Every adapter is a **subprocess** in its own environment (ADR 0021), with a
  two-subcommand CLI: `fit --train --valid --config --run-dir --controls` and
  `predict --cases --checkpoint --out --controls`. No benchmark module imports
  model internals.
- Case bundles are shared per fold: train/valid carry expression and Δ; the test
  bundle is **metadata only** (no `X`, `log1p`, Δ, or expression row index) and
  `cases.verify_cases` fails closed on a leak. Test well IDs are visible; adapters
  must not join them back to global expression.
- Per unit: `config.json`, `run.log`, `checkpoint/`, `pred_valid.npy`,
  `pred_test.npy`, group predictions/truth, `metrics.json`, `valid_metrics.json`,
  `status.json`, and artifact checksums. A unit is resumed only when
  `status == completed` and all artifacts exist.
- The leaderboard, `metrics.csv`, `per_drug.csv`, `bootstrap.npz`, and
  `run_cards.csv` are generated from the raw fold predictions by
  `evaluate.py`/`summarize_run.py`, never hand-edited.

### 7. Quarantine

No 978-landmark artifact is read in training or evaluation (ADR 0020/0024).
PRnet's 978 demo remains setup-only. No Stage 1 code or data is imported.

## Alternatives rejected

- **Warm-drug table:** would reward models that memorise drug identity and would
  not answer the MIRI question (new drugs). One cold-drug leaderboard instead.
- **978 landmarks:** quarantined and too narrow for whole-transcriptome MIRI
  perturbomes.
- **Micro/pooled AUPRC:** dominated by large groups and by the mean-expression
  baseline.
- **Per-model HPO:** would confound model quality with tuning effort and break the
  16 h budget.
- **Released checkpoints:** would make the comparison a checkpoint-selection
  exercise and leak author pretraining data across the cold-drug split.
- **A second "ideal-form" table:** encourages reading a partial model as a
  winner; `ideal_fit` is a column.

## Consequences

- The leaderboard is a **comparative, contract-bound** result, not a reproduction
  of any author's reported score. Every model carries a fidelity ledger; the
  report must state what the benchmark does and does not show.
- The two exclusions are documented, not rows: scVIDR and DECANT have no
  cold-drug path, and a new architecture is out of scope.
- The benchmark fixes the evidence for choosing the Stage 2 model; it does not
  prove clinical efficacy, and it does not start the MIRI adaptation. That is the
  next effort.
