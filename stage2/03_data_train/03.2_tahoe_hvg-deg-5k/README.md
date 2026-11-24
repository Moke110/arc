# Tahoe HVG-DEG 5k well-level data

This source-specific route exposes the accepted Tahoe well-level H5AD,
membership and exclusion ledgers from the [sealed combined full view](../03.1_hvg-deg-5k/README.md).
The three assets keep their original digests. The new manifest records the
combined manifest digest and the original axis and corpus generations.

Add `stage2/03_data_train/03.1_hvg-deg-5k/output/api` to `PYTHONPATH` and
call `expression_view.open_view()` with this substep's `output/release`
directory and `source="tahoe_bulk"`. This reader returns at most 256 rows
of well-level `X`, `control` or `delta` per call.

```bash
BUILDER=stage2/03_data_train/03.2_tahoe_hvg-deg-5k/execution/publish_tahoe.py
python "$BUILDER" publish
python "$BUILDER" verify
```

The builder hard-links the immutable archive assets on the same filesystem
and copies them across filesystems. Pass `--sealed` when the archive payloads
are materialized in a different checkout.
