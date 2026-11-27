# 03.5 tenv Denoise-DEG 3k training

Read the [Step 03 operating guide](../AGENTS.md) and the [release summary](README.md).

## Build and verify

1. Run the [tenv builder](../execution/build_tenv_training.py) with `--axis denoise-3k`
   using the [Step 03 build procedure](../README.md#tenv-training-axes). Use
   [native tenv integer counts](../../01_data_raw/README.md), the published
   [ordered 3k axis](../../02_gene_axis/02.2_denoise-deg/output/release/manifest.json),
   the [frozen Step 07 tenv folds](../output/frozen-folds/tenv/manifest.json), and
   the accepted [HVG-DEG 5k population reference](../03.1_tenv_hvg-5k/output/manifest.json).
2. Project counts to 3k in the published order. Pool DMSO by plate and cell line,
   and treatments by the same keys plus canonical parent and measured dose. Apply
   CP10K and natural `log1p` using each pooled profile's 3k count total. Compute Δ
   by subtracting the normalized matched DMSO profile from the normalized treatment.
3. Build outside `output/` and run `verify` for byte-exact replay. Require equality
   with the 5k reference for ordered control and treatment keys, control links,
   raw membership and exclusions, all six fold assignments and strip reasons.
   Publish to `output/` with `--expected-sha256` set to the verified manifest digest.
   Compare the [release manifest](output/manifest.json) with the
   [six-axis acceptance receipt](../log/313-six-tenv-axes.json).

Open a materialized release with `open_tenv_axis_training("denoise-3k")` from
the [shared tenv axis API](../output/api/tenv_axes.py). Read control and treatment
windows through the [bounded reader](../03.1_tenv_hvg-5k/README.md#downstream-use-and-verification).
