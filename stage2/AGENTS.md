# AGENTS.md — Stage 2

Before working in a step, read its `AGENTS.md` and `README.md`, ADR 0062,
and the scientific ADRs it names. Apply the repository execution-environment rules.

## Active steps

- `01_bases/` — shared foundation: frozen identity, corpora, evaluator truth,
  fold/seed products, metrics, palette and figure APIs. Every later step imports
  this implementation instead of re-deriving it.
  - Gene axis (ADR 0066): the step releases two axis versions, `G_2000` and
    `G_5000`, from the per-well tenv-HVG dispersion and cardiomyocyte-derived
    Control-vs-Disease DEGs fused at two levels (datasets inside a disease category,
    then categories), selected by a 50/50 tenv-HVG / pN quota. `output/` is being
    rebuilt as `raw_counts/` → `corpus/` (integer raw counts, no normalization, no
    Δ) → `view-2k/`/`view-5k/` (full-axis normalization, then subset; `X`, `control`,
    `Δ`), with `disease-deg/` and `gene-axis/`; the measured truth and `log/` are
    retired pending an evaluation redesign. Routing: `stage2/01_bases/AGENTS.md`.
- `02_models_tenv/` — tenv model comparison (ADR 0064, superseding the indirect
  and figure portions of ADR 0063): the PRnet entrant and the `Mean` reference,
  trained and evaluated on `tenv`. Bounded fit/predict admission is #163; the
  shared evaluation driver is #164; the accepted 2 folds x 2 seeds matrix and its
  `execution/run_matrix.py` publish/replay are #165; the @d reliability family is
  #167, the driver/tables #168, the three figures #169, and the no-training
  acceptance/replay plus retirement #170. The three report figures are
  `figure/metrics.png` (four direct-metric panels), `figure/metrics_indirect.png`
  (two-by-four reliability panels) and `figure/learning_curve.png`; each has one
  tracked `figure/data/*.csv` plot-value table.

## Step contract
- `output/`: minimal accepted products reused by other steps: corpora, gene/HVG tables,
  evaluator truth, public metric/reader code and reusable model checkpoints. Each has
  one stable path and one owner; public implementations live here and execution imports them.
- `figure/*.png`: human-facing report figures.
- `figure/data/{figure}_{panel}.csv`: complete plotted values for each panel, including
  curve coordinates, groups, uncertainty and reference values where applicable.
- `figure/script/`: reproduce PNGs from those CSVs alone, without matrices, checkpoints
  or recomputing scientific metrics. Record provenance when a CSV derives from output.
- `execution/`: reproducible build/audit code, configurations, tests and disposable
  working intermediates. Keep necessary reconstruction dependencies available.
- `log/`: execution/audit records, environments, seeds, acceptance results, hashes,
  resource measurements, migration mappings and compact failure evidence.
- `AGENTS.md`: agent execution/reproduction workflow; `README.md`: human-readable
  accepted results, limitations and links to figures and reusable products.

## Acceptance and replacement

Build candidates outside formal output, verify structure/data constraints, scientific
criteria, downstream reads and CSV-only figure reproduction, then replace stable paths.
Record metric tolerances for stochastic results and hashes for file integrity. A failed
candidate leaves the accepted products available. Retain only accepted final products;
remove obsolete payloads after acceptance and reference migration, preserving required
provenance and audit evidence. Historical sealed records keep their original bytes and
paths; log a migration mapping. External `source/` trees and Stage 1 remain untouched.
