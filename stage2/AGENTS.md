# AGENTS.md — Stage 2

Before working in a step, read its `AGENTS.md` and `README.md`, ADR 0062,
and the scientific ADRs it names. Apply the repository execution-environment rules.

## Active steps

- `01_bases/` — shared foundation: frozen identity, corpora, evaluator truth,
  fold/seed products, metrics, palette and figure APIs. Every later step imports
  this implementation instead of re-deriving it.
- `02_models_tenv/` — tenv model comparison (ADR 0063): the PRnet entrant and the
  `Mean` reference, trained and evaluated on `tenv`. Bounded fit/predict
  admission is #163; the 2 x 2 run matrix is #165; report figures are #164.

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
