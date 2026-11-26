# Tahoe HVG-DEG 5k well-level data

This release retains the accepted Tahoe well-level H5AD, membership and
exclusion ledgers. Its manifest records the original combined manifest digest
and axis and corpus generations.

Add this substep's `output/api` and `stage2_CarPS-VAE/util/output/api` to `PYTHONPATH` and
call `tahoe_expression.open_tahoe()`. The reader returns at most 256 rows of
well-level `X`, `control` or `delta` per call.

```bash
BUILDER=stage2_CarPS-VAE/03_data_train/03.2_tahoe_hvg-5k/execution/publish_tahoe.py
python "$BUILDER" verify
```

For historical reconstruction, pass `--sealed <combined-release-directory>` to
`publish`. The builder hard-links the sealed assets on the same filesystem and
copies them across filesystems. The three assets keep their original digests.
