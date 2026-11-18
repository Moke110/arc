# AGENTS.md — Stage 2

Before working in a step, read its `AGENTS.md` and `README.md`, ADR 0062,
and the scientific ADRs it names. Apply the repository execution-environment rules.

## Active steps

- `01_bases/` — shared foundation (ADR 0066 integer corpus; ADR 0073 disease/axis
  contract): frozen identity, raw counts, corpus, paired 5k expression views,
  complete disease source evidence, tier-first DEG1000 fusion, 2k/3k/4k/5k nested
  gene axes, fold/seed products, palette and the single public API. Later steps
  import it rather than re-deriving it. Axes use half tenv-HVG, half unused
  f1000-first/rank_sum-ordered disease candidates, with explicit HVG fallback.
  ADR0077 views normalize on selected G_5000 integer totals and retain matched
  control/delta semantics. Old02/03 experiments keep their full-axis scale. The human rule DOCX is at the step root.
  Routing: `stage2/01_bases/AGENTS.md`.
- `02_models_tenv-5k/` — single-dose tenv native model comparison (ADR 0074):
  shared splits and Mean baselines, random-weight PRnet/TranSiGen/MLP_Morgan,
  parent-level Δ evaluation, checkpoint replay and CSV-only report figures.
  Before changing or replaying it, read its `AGENTS.md` and `README.md`.

- `03_tsg_dose/` — TranSiGen dose-input optimization round1 (ADR0076 amends
  ADR0075):126 cold-drug/held-dose fits, validation-first selection and three
  report PNGs. Historical final fits/weights are retained for audit; this step
  supplies no downstream weights. Read its `AGENTS.md` and `README.md`.

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
