# tenv Denoise-DEG 5k training

Read the [Step 03 routing](../AGENTS.md) and the [release summary](README.md).

1. Build a scratch candidate using `publish --axis denoise-5k`. Use the
   [source-native tenv counts](../../01_data_raw/README.md), the
   [ordered Denoise-DEG 5k axis](../../02_gene_axis/02.2_denoise-deg/output/release/manifest.json), the
   [frozen tenv folds](../output/frozen-folds/tenv/manifest.json) and
   [Step 07 fold receipt](../output/frozen-folds/247-axis-local-smoke.json), and the
   [accepted HVG-DEG 5k tenv reference](../03.1_tenv_hvg-5k/output/manifest.json).
   [Builder](../execution/build_tenv_training.py) · [Build and replay procedure](../README.md#tenv-training-axes)
2. Project native integer counts to the ordered axis before pooling. Normalize
   each pooled profile with CP10K on its selected-axis total and natural `log1p`.
   Compute Δ against the normalized DMSO control from the same plate.
3. Run `verify --axis denoise-5k` for byte-exact replay. Require equality with
   HVG-DEG 5k for the ordered plate and treatment keys, matched control rows,
   membership and exclusions, and all six fold assignments and strip reasons.
   Pass the verified candidate manifest digest as `--expected-sha256` when
   publishing to the stable release. [Acceptance receipt](../log/313-six-tenv-axes.json) · [Release manifest](output/manifest.json)

Read the materialized release with `tenv_axes.open_tenv_axis_training("denoise-5k")`.
Pass the Step 03 root to use a release in another checkout. [Shared reader](../output/api/tenv_axes.py)
