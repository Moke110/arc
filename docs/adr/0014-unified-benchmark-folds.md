# Unified benchmark fold scheme across v1 / v1.1 / v2

The three MIRI workflows (`workflow/xpert-fda-miri-v1`, `-v1.1`, `-v2`)
benchmarked XPert with inconsistent fold schemes, so their benchmark numbers
were not comparable and the evaluation's meaning was underspecified. This ADR
defines **one fold principle** and applies it identically (same logic, same
seed 2024, same ratio) to all three workflows, regenerating every fold column
and making each workflow self-contained.

## The fold principle (per workflow, on that workflow's finetune h5ad)

- **warm-start (10 folds)** — *drug-overlap*: every drug in the test set also
  appears in the train set. Implemented as a per-drug sample split: each
  drug's samples are randomly assigned `max(1, round(0.2·n))` to test, the
  rest to train (single-sample drugs go to test); any sample never in test is
  rescued into a random fold. Column values are only `train` / `test` (plus
  `''` for excluded). No `valid`.
- **cold-drug (10 folds)** — *drug-disjoint*: no drug in the test set appears
  in the train set. Each fold sample-weight-randomly selects drugs whose total
  samples ≈ 20% of all valid samples as test; every other drug is train. A
  coverage rescue guarantees every drug is in test ≥ 1 fold. Column values are
  only `train` / `test` (plus `''`). No `valid`.
- **No category or data-source splits** — `cold-cate` (v1's 6-fold
  leave-one-drug-category-out on the DToxS KI MOA) is dropped everywhere; no
  fold is stratified on `moa`/drug category or `dataset_id`.
- **Ratio** — per fold, `N(test) : N(train) ≈ 1:4` (target `test_frac =
  0.2`), inside the allowed band `1:9 – 3:7`. The split scripts assert and log
  the per-fold test fraction.
- **No validation set** — finetune epochs are fixed, so fold columns carry
  only `train`/`test`; the XPert loader's `valid → test` fallback is relied
  on (no library change).
- **Column lifecycle** — regenerating splits drops every prior fold/full
  column first, so each h5ad ends with exactly `warm-start_fold-1..10`,
  `cold-drug_fold-1..10`, and the all-train full column
  (`dtoxs_full` / `cm_full`).
- **Deterministic** — fold generation is reproducible under `seed = 2024`.
  Values are full strings (`train`/`test`/`''`), never a fixed-width numpy
  dtype (which previously truncated `'test'` → `'t'`).
- **Invariant gate** — each workflow's split script (the single test seam)
  asserts after writing: value set ⊆ {`train`, `test`, `''`}; warm-start drug
  overlap = 1.0; cold-drug drug overlap = 0.0; per-fold test fraction ∈
  [0.10, 0.30]; every valid drug (and warm-start sample) in test ≥ 1 fold; no
  stale fold columns. Any failure exits 1.

## v1 re-baselining

v1's frozen `input/dtoxs.h5ad` fold columns were regenerated **in place**
(folds only — expression, schema, and `dtoxs_full` are unchanged), which
re-bases v1's archived benchmark: the previously archived 26-fold table
(warm-start 10 + cold-drug 10 + cold-cate 6) is replaced by the 20-fold
warm-start + cold-drug table. Downstream steps are unaffected because they
consume the all-train model / predictions, not the fold columns. v1's
`create_dtoxs_splits.py` also fixes the legacy `dtype=str`/`<U1` truncation
bug and now generates cold-drug as well, so v1's folds are reproducible by
code rather than frozen by hand.

## v2 self-containment

v2's union h5ad had **no split columns on disk** (step 03 was rebuilt after
step 05 ran). Its split script now regenerates the unified columns (fixing
that stale state), and all cross-workflow comparison tooling is removed:
`eval_dtoxs_folds.py` (union warm-start weights evaluated on v1's DToxS folds),
`compare_v1_benchmark.py`, their eval config, and their outputs. Each workflow
now stands alone; the step-03 schema/parity gates (data integrity vs v1's
h5ad, not benchmarks) are preserved. v2's benchmark also aligns to **4
finetune epochs** (was 3), matching v1/v1.1.

## Consequences

- v1, v1.1, and v2 each bench `warm-start` + `cold-drug` (20 folds), so the
  benchmark table has the same shape everywhere and cold generalization is
  measured identically.
- `benchmark.png` in all three workflows overlays each fold as a black dot
  (jittered) on the corresponding mean ± std bar.
- GPU steps re-run the XPert finetune/benchmark; new benchmark numbers are
  environment-specific (GPU-fp-nondeterministic), exactly as before.
  Compare-to-history should cite this ADR.
