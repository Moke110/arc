# 03.3 tenv HVG-DEG 3k training

Read the [local report](README.md) and [Step 03 routing](../AGENTS.md).
This release uses the Step 03 tenv builder with `--axis hvg-3k`.

1. Resolve the [native tenv counts](../../01_data_raw/README.md), the ordered
   [HVG-DEG 3k axis](../../02_gene_axis/02.1_hvg-deg/output/release/manifest.json),
   the [frozen tenv folds](../output/frozen-folds/tenv/manifest.json) and
   [original Step 07 receipt](../output/frozen-folds/247-axis-local-smoke.json).
   The fold ledger retains its HVG-DEG 5k provenance. Use the accepted
   [03.1 release](../03.1_tenv_hvg-5k/output/manifest.json) as the population
   reference and check the input pins in the [03.3 manifest](output/manifest.json).
2. Run [build_tenv_training.py](../execution/build_tenv_training.py) with
   `publish --axis hvg-3k` to a scratch output. Follow the
   [shared build procedure](../README.md#tenv-training-axes) for staging and
   linked-worktree inputs. Pool compatible DMSO counts by plate and treatments
   by plate, canonical drug parent and measured dose. Project the integer
   counts to the ordered 3k axis before pooling. Apply CP10K and natural
   `log1p` to each pooled profile using its 3k library total, then compute Delta
   against the normalized matched plate control.
3. Run `verify --axis hvg-3k` on the candidate. Require byte-exact replay and
   equality with 03.1 for ordered control and treatment keys, matched control
   rows, membership and exclusion ledgers, and all six fold assignments and
   strip reasons. Promote to stable `output/` with `--expected-sha256` from
   the successful verification. The builder uses the
   [shared publication API](../../util/README.md). See the
   [six-axis acceptance receipt](../log/313-six-tenv-axes.json) for the
   recorded replay and population checks.

For downstream reads, use `open_tenv_axis_training("hvg-3k")` from the
[shared axis reader](../output/api/tenv_axes.py). Pass the Step 03 root to read
an installed release in another checkout. Use the
[tenv reader](../output/api/training_data.py) for bounded matrix windows and
read each treatment's matched control through its `control_row`.
