# tenv HVG-DEG 2k training

Read the [Step 03 routing](../AGENTS.md) and the [release summary](README.md).

1. Build a scratch candidate with `publish --axis hvg-2k`. Use the
   [source-native tenv counts](../../01_data_raw/README.md), the
   [ordered HVG-DEG 2k axis](../../02_gene_axis/02.1_hvg-deg/output/release/manifest.json),
   the [frozen tenv folds](../output/frozen-folds/tenv/manifest.json) and
   [Step 07 fold receipt](../output/frozen-folds/247-axis-local-smoke.json), and the
   [accepted HVG-DEG 5k tenv reference](../03.1_tenv_hvg-5k/output/manifest.json).
   [Builder](../execution/build_tenv_training.py) · [Build and replay procedure](../README.md#tenv-training-axes)
2. Project the native integer counts onto the ordered 2k axis before pooling.
   Apply CP10K on each pool's selected-axis total followed by natural `log1p`.
   Compute Δ against the normalized DMSO control from the same plate.
3. Run `verify --axis hvg-2k` for byte-exact replay. Require equality with
   HVG-DEG 5k for the ordered control and treatment keys, matched control rows,
   membership and exclusions, and all six fold assignments and strip reasons.
   Pass the verified candidate digest as `--expected-sha256` when publishing
   to the stable output. [Release manifest](output/manifest.json)
4. After both 2k releases are available, run the independent audit. Require
   exact agreement for ordered genes, pooled library totals, normalized profiles,
   matched Δ and the overall control vector. Confirm the six established releases
   still load. [Audit command](../execution/audit_tenv_2k.py) · [Acceptance receipt](../log/366-two-tenv-axes.json)

Read the materialized release with `tenv_axes.open_tenv_axis_training("hvg-2k")`.
Pass the Step 03 root to read a release in another checkout. [Shared reader](../output/api/tenv_axes.py)
