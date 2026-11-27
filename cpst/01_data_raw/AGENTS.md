# CarPS-Transcriptome raw data

Read the [step README](README.md) and the ownership decisions in
[ADR 0083](../../docs/adr/0083-stage2-corpus-axis-view-boundaries.md) and
[ADR 0085](../../docs/adr/0085-stage2-raw-and-training-data.md).

Preserve nonnegative integer counts, observation identities and native gene
columns at the finest available source granularity. Retain vehicles, blanks,
unmatched observations and their match metadata. tenv and Tahoe bulk require
complete common-axis coverage. For DToxS-primary, align shared counts by native
gene ID and preserve every original table membership and control context.
[Step 02](../02_gene_axis/AGENTS.md) owns gene selection.
[Step 03](../03_data_train/AGENTS.md) owns pooling, normalization, Δ and training
exclusions.

## Read the accepted release

Use `open_corpus(...)` from [`native_corpus.py`](output/api/native_corpus.py) with
the materialized `output/` directory. Its `manifest.json` pins the seven
dataset directories. The reader checks each dataset manifest and count payload
before reading it. Use `counts(start, stop)` for windows of at most 256 rows or
`batches()` to iterate. Use `native_columns()` to resolve ordered common-gene
indices and `matched_control_rows()` to resolve vehicles. Unmeasured genes raise
an error and unmatched treatments return no vehicle rows.
[Accepted release digests](log/migration-303.json)

## Build and verify

- Run [`build_corpus.py`](execution/build_corpus.py) from the repository root in
  `arc`. For a fresh build, use pinned count and annotation inputs in the
  [perturbation data lake](../../pert_data/README.md), including the declared local
  compound annotation, and the tracked [gene identity release](../util/output/gene-identity/manifest.json).
  Set `--materialized-root` to the checkout holding those data-lake inputs.
- Run `verify` against the installed migrated release to check its inventory,
  digests and coverage. Match its printed digest to the acceptance receipt.
  Default `publish` on that installation also checks the
  active training consumers. Use an unused `--output` for a new candidate and
  verify the corresponding Step 03 releases before switching generations.
  Pass an explicit sealed v1 `--predecessor` to `verify` for equivalence replay.
  Pass the same input and a separate scratch output to `migrate`.
- For Tahoe scans, follow the [execution limits](../../docs/agents/execution-environment.md),
  run one heavy job at a time and use bounded reads. The builder checks source
  pins before reusing plate checkpoints and independently replays any staged
  `--tahoe-built` input. Keep source files read-only.
