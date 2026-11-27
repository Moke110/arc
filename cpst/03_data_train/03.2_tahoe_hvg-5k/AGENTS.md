# Tahoe HVG-DEG 5k well-level data

Read the [Step 03 routing](../AGENTS.md) and this step's [release report](README.md).

## Read and verify the accepted release

- Use [`tahoe_expression.open_tahoe(path=None)`](output/api/tahoe_expression.py)
  as a context manager. Put this step's `output/api` and `cpst/util/output/api`
  on `PYTHONPATH`. Pass the materialized release directory when the active
  worktree contains only the tracked manifest. Use `read` or `batches` for
  `X`, `control`, or `delta` in windows of at most 256 rows.
- Preserve the ordered [HVG-DEG 5k axis](../../02_gene_axis/02.1_hvg-deg/README.md)
  and the sealed `tahoe_bulk` full-view recipe. Keep treatment rows at well
  granularity, sum matched vehicle counts once per source observation, apply
  CP10K and `log1p` separately on the selected axis, and retain `delta = X - control`.
  The [release manifest](output/release/manifest.json) pins the axis and corpus
  generations and the H5AD, membership ledger and exclusion ledger.
- From the repository root, run
  `python cpst/03_data_train/03.2_tahoe_hvg-5k/execution/publish_tahoe.py verify`.
  Use `--output` for a materialized release in another checkout. Require `PASS`
  and the manifest and asset digests recorded in the
  [acceptance receipt](../log/312-source-specific-releases.json).

## Reconstruct from a sealed combined release

Supply the original combined full-view manifest and its three Tahoe assets.
Match its manifest digest to `sealed_manifest_sha256` in the accepted Tahoe
manifest. With [`publish_tahoe.py`](execution/publish_tahoe.py), run `publish`
with `--sealed` and separate scratch `--output` and `--staging` directories.
Verify the candidate with `verify --sealed` and its `--output`. Require `PASS`
and equality with the accepted manifest and its three asset digests before
targeting the stable output. The builder uses the
[shared release and publication APIs](../../util/README.md).
