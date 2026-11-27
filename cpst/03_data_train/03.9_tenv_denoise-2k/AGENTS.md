# tenv Denoise-DEG 2k training

Read the [Step 03 routing](../AGENTS.md) and the [release summary](README.md).

1. Build a scratch candidate with `publish --axis denoise-2k`. Use the
   installed [native tenv counts](../../01_data_raw/README.md), the
   [ordered Denoise-DEG 2k axis](../../02_gene_axis/02.2_denoise-deg/output/release/manifest.json),
   the [frozen tenv folds](../output/frozen-folds/tenv/manifest.json) and
   [Step 07 fold receipt](../output/frozen-folds/247-axis-local-smoke.json), and the
   [accepted HVG-DEG 5k tenv reference](../03.1_tenv_hvg-5k/output/manifest.json).
   The fold provenance stays bound to HVG-DEG 5k.
   [Builder](../execution/build_tenv_training.py) · [Build and replay procedure](../README.md#tenv-training-axes)
2. Project the native integer counts onto the ordered 2k axis before pooling.
   Use each pool's selected-axis total for CP10K, then apply natural `log1p`.
   Compute Δ against the normalized DMSO control from the same plate.
3. Run `verify --axis denoise-2k` for byte-exact replay. Require equality
   with HVG-DEG 5k for ordered control and treatment keys, matched control
   rows, membership, exclusions, all six fold assignments and strip reasons.
   Pass the verified manifest digest as `--expected-sha256` when publishing
   to the stable output. [Release manifest](output/manifest.json)
4. After both 2k releases are available, run the independent audit against
   native counts. Require exact agreement for ordered genes, pooled library
   totals, normalized profiles, matched Δ and the overall control vector.
   Confirm the six established releases still load.
   [Audit command](../execution/audit_tenv_2k.py) · [Acceptance receipt](../log/366-two-tenv-axes.json)

Use `tenv_axes.open_tenv_axis_training("denoise-2k")` to read the materialized
release. Pass the Step 03 root to read it in another checkout.
[Shared reader](../output/api/tenv_axes.py)
