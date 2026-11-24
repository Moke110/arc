# AGENTS.md — Stage 2

Before working in a step, read its `AGENTS.md` and `README.md`, ADR 0062,
and the scientific ADRs it names. Apply the repository execution-environment rules.

## Active steps

- [util](util/AGENTS.md) owns gene and drug identities, representations, ECFP4,
  seed selection, generic vector metrics, the palette and CSV-only plotting.
  It also owns release publication and locking, shared DEG and ranking kernels,
  R worker resolution, and derivation manifest verification. Read the
  [public API inventory](util/README.md) and reuse compatible operations.
  Keep fold validation, metric aggregation and report tables with the study.
- [01_data_raw](01_data_raw/AGENTS.md) publishes source-native integer counts, coverage
  and match status.
- [02_gene_axis](02_gene_axis/AGENTS.md) owns the four ordered HVG-DEG prediction
  axes and their selection evidence.
- [03_view](03_view/AGENTS.md) publishes the full raw-first 5k expression release
  for downstream inputs.
- [04_model_view](04_model_view/AGENTS.md) compares single-dose tenv models.
- [05_tsg_loss](05_tsg_loss/AGENTS.md) reserves the loss comparison.
- [06_tsg_drug_emb](06_tsg_drug_emb/AGENTS.md) compares drug representations.
- [07_tsg_dose](07_tsg_dose/AGENTS.md) compares dose methods in Tahoe and tenv.
- [08_tsg_width_dose](08_tsg_width_dose/AGENTS.md) compares width and dose encodings.
- [09_tsg_transfer](09_tsg_transfer/AGENTS.md) compares Tahoe-to-tenv transfer.

Use each study's `execution/raw_first.py` for new inputs and fixed-prediction
evaluation. Use `figure/raw_first/render.py` or the documented `render` command
for CSV-only rendering with util. The [historical evidence index](03_view/03.1_hvg-deg-5k/historical/README.md) locates sealed records and their original provenance.

## Step contract
- `output/`: minimal accepted products reused by other steps: corpora, gene/HVG tables,
  evaluator truth, public metric/reader code and reusable model checkpoints. Each has
  one stable path and one owner; public implementations live here and execution imports them.
- `figure/*.png`: human-facing report figures.
- `figure/data/{figure}_{panel}.csv`: complete plotted values for each panel, including
  curve coordinates, groups, uncertainty and reference values where applicable.
- `figure/script/`: reproduce PNGs from those CSVs alone, without matrices, checkpoints
  or recomputing scientific metrics. Record provenance when a CSV derives from output.
  Use a named `figure/<successor>/` subtree when a receipt seals the predecessor's
  paths. Bind the successor inventory to its own acceptance receipt.
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
