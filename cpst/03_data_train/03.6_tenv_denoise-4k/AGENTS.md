# 03.6 tenv Denoise-DEG 4k training

Read the [Step 03 operating guide](../AGENTS.md) and [release summary](README.md).

1. Use the [tenv builder](../execution/build_tenv_training.py) with
   `--axis denoise-4k` and follow the [Step 03 build procedure](../README.md#tenv-training-axes).
   Resolve the [native tenv integer corpus](../../01_data_raw/README.md), the
   [ordered Denoise-DEG 4k axis](../../02_gene_axis/02.2_denoise-deg/output/release/manifest.json),
   the [frozen tenv folds](../output/frozen-folds/tenv/manifest.json) with their
   [original Step 07 receipt](../output/frozen-folds/247-axis-local-smoke.json),
   and the [HVG-DEG 5k population reference](../03.1_tenv_hvg-5k/output/manifest.json).
2. Project the native counts in published gene order before pooling. Pool DMSO
   by plate and cell line, and treatments by the same keys plus canonical drug
   parent and measured dose. Apply CP10K and natural `log1p` using each pooled
   row's 4k count total. Compute Δ from the normalized treatment and its
   matched normalized DMSO row.
3. Build a scratch candidate and run `verify` for byte-exact replay. Require
   agreement with the 5k reference on ordered plate and treatment keys, control
   links, membership and exclusion ledgers, all six fold assignments and strip
   reasons. To reproduce the accepted release, match the digest of its
   [manifest](output/manifest.json) and the digest recorded in the
   [six-axis acceptance receipt](../log/313-six-tenv-axes.json). Publish to
   stable `output/` with `--expected-sha256` set to the verified candidate digest.

Open the installed release with `open_tenv_axis_training("denoise-4k")` from the
[shared tenv axis API](../output/api/tenv_axes.py). Use its
[bounded reader](../03.1_tenv_hvg-5k/README.md#downstream-use-and-verification)
for control, treatment and Δ windows.
