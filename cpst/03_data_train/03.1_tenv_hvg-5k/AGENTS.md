# 03.1 tenv HVG-DEG 5k training

Read the [release contract and commands](README.md), [Step 03 routing](../AGENTS.md)
and [ADR 0062](../../../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md).
03.1 supplies `hvg-5k` and the reference population for the other tenv axes.

1. Resolve the accepted [native tenv corpus](../../01_data_raw/README.md),
   ordered [HVG-DEG 5k axis](../../02_gene_axis/02.1_hvg-deg/output/release/manifest.json),
   and original [Step 07 fold receipt](../output/frozen-folds/247-axis-local-smoke.json)
   with its [fold manifest](../output/frozen-folds/manifest.json). Check those
   pins in the [03.1 manifest](output/manifest.json). In a linked worktree, use
   `--raw` for the installed native corpus and keep the accepted axis available
   in the executing checkout for fold validation.
2. Run [build_tenv_training.py](../execution/build_tenv_training.py) with
   `publish --axis hvg-5k` to a scratch candidate. Apply the
   [execution environment limits](../../../docs/agents/execution-environment.md)
   and run build and replay serially. Pool compatible DMSO counts by plate
   and treatments by plate, canonical drug parent and measured dose. Preserve
   one cell line per plate. Project to the ordered axis before summing integer
   counts, then apply CP10K and natural `log1p` using each pooled row's selected-axis
   total. Delta uses the normalized matched plate control, and `C` gives each
   normalized plate control equal weight.
3. Run `verify` on the candidate and require byte-exact replay. Preserve the
   member and exclusion ledgers and require agreement of member fold labels
   within each pooled treatment. In cold-dose-drug folds, the held parent and
   held dose intersect in test, train and valid come from outside their union,
   and each remaining held strip has an empty assignment with its reason.
   Publish to stable `output/` with `--expected-sha256` from the successful
   verification. The builder uses the [shared publication API](../../util/README.md).

Use the [shared tenv reader](../output/api/training_data.py),
`open_tenv_training(root=None)`, for downstream checks. Verify `control_row`
against the control table and read controls or treatments in windows of at
most 256 rows. The [acceptance receipt](../log/312-source-specific-releases.json)
records the formal digest, partition checks and replay.
