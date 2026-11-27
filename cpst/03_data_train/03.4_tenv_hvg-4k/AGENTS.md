# 03.4 tenv HVG-DEG 4k training

Read the [Step 03 operating guide](../AGENTS.md) and the [release summary](README.md).

1. Use [the tenv builder](../execution/build_tenv_training.py) with `--axis hvg-4k`
   and the [Step 03 build procedure](../README.md#tenv-training-axes). Supply the
   [native tenv counts](../../01_data_raw/README.md), the published
   [ordered HVG-DEG axes](../../02_gene_axis/02.1_hvg-deg/output/release/manifest.json), the
   [frozen Step 07 tenv folds](../output/frozen-folds/tenv/manifest.json), and the
   accepted [HVG-DEG 5k reference](../03.1_tenv_hvg-5k/output/manifest.json).
   Project integer counts to 4k before pooling by plate for DMSO and by plate,
   canonical parent and measured dose for treatments. Compute CP10K and natural
   `log1p` on each pooled 4k profile, then subtract the normalized same-plate
   DMSO profile for Δ.
2. Build outside `output/` and run `verify` for byte-exact replay. Check that the
   ordered control and treatment keys, raw membership and exclusion ledgers, and
   all six fold assignments and strip reasons match the 5k reference. Publish to
   `output/` with `--expected-sha256` set to the verified candidate digest. The
   [acceptance receipt](../log/313-six-tenv-axes.json) records the checks for the
   [published manifest](output/manifest.json).

Read a materialized release with `open_tenv_axis_training("hvg-4k")` from
[the shared tenv axis API](../output/api/tenv_axes.py). Use the
[bounded read methods](../03.1_tenv_hvg-5k/README.md#downstream-use-and-verification)
for control and treatment slices.
