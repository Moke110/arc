# Step 03 training data workflow

- **tenv plate-pooled release.** Pick an axis ID from the [tenv release index](README.md#tenv-training-axes)
  and read that release's guide. Use [03.1](03.1_tenv_hvg-5k/AGENTS.md) for the
  `hvg-5k` population reference and the shared build procedure below for any axis.
- **Tahoe well-level release.** Use [03.2](03.2_tahoe_hvg-5k/AGENTS.md) to verify
  the installed HVG-DEG 5k view or reconstruct it from the sealed combined view.
- **Archive independence.** Use the active-input smoke procedure in the
  [Step 04 guide](../04_tsg_axis/AGENTS.md) to verify the materialized inputs
  with the retired archive unavailable. [Smoke entry point](../04_tsg_axis/execution/verify_active_inputs.py)

## Build and replay a tenv release

1. Resolve the [Step 01 native corpus](../01_data_raw/AGENTS.md), the selected
   [Step 02 ordered axis](../02_gene_axis/AGENTS.md), the [frozen tenv fold manifest](output/frozen-folds/tenv/manifest.json)
   and the original [Step 07 receipt](output/frozen-folds/247-axis-local-smoke.json).
   Check these pins against the target release manifest. The builder validates
   folds against the accepted HVG-DEG 5k axis in the executing checkout. For
   each other axis, supply the accepted [03.1 population reference](03.1_tenv_hvg-5k/output/manifest.json).
2. Project native integer counts onto the ordered axis before pooling DMSO
   by plate and treatments by plate, canonical drug parent and measured dose.
   Preserve one cell line per plate and match controls through `control_row`.
   Apply CP10K on each pooled row's selected-axis total and natural `log1p`.
   Subtract the normalized matched control to obtain Δ. Give each normalized
   plate control equal weight in `C`. Keep all six frozen fold assignments and
   their strip reasons. Preserve the ordered keys, control links, membership
   and exclusion ledgers across axes. [Data boundary](../../docs/adr/0085-stage2-raw-and-training-data.md)
3. Build and replay serially with the repository RAM cap and single-thread
   limits. Use separate scratch candidate and staging directories. Run this
   example from the repository root, substituting the axis ID and stable
   destination from the leaf guide as needed.
   [Execution environment](../../docs/agents/execution-environment.md) · [Builder](execution/build_tenv_training.py)

   ```bash
   step03=$PWD/cpst/03_data_train
   step03_work=$PWD/tmp/step03-denoise-3k
   export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
   python "$step03/execution/build_tenv_training.py" publish --axis denoise-3k \
     --output "$step03_work/candidate" --staging "$step03_work/staging" &&
   python "$step03/execution/build_tenv_training.py" verify --axis denoise-3k \
     --output "$step03_work/candidate" --staging "$step03_work/staging"
   ```

   Use `--raw`, `--axis-root`, `--fold-root`, `--fold-receipt` or `--reference`
   to resolve inputs in another checkout. `--axis-root` selects the release
   directory for the requested axis. Keep the accepted fold axis available
   through the executing checkout's default Step 02 path.
4. Require `PASS` from byte-exact replay and population agreement with 03.1
   for alternate axes. To reproduce an accepted release, match its manifest
   digest to the leaf receipt. When promoting, publish with `--expected-sha256`
   set to `manifest_sha256` from the successful `verify`. The builder validates
   the digest before installing through the shared publication API.
   [Shared API](../util/README.md#shared-input-and-release-apis)

   ```bash
   python "$step03/execution/build_tenv_training.py" publish --axis denoise-3k \
     --output "$step03/03.5_tenv_denoise-3k/output" --staging "$step03_work/staging" \
     --expected-sha256 "${VERIFIED_SHA256:?Set VERIFIED_SHA256 from the successful replay}"
   ```

## Audit the two 2k releases

After both 2k releases are materialized, run
`python cpst/03_data_train/execution/audit_tenv_2k.py` from the repository root.
Use `--raw` for an alternate native corpus directory, `--training` for the
Step 03 root containing the two 2k releases, and `--established` for the Step 03
root containing the six established releases. Require exact agreement with
native counts for ordered genes, pooled totals, normalized profiles, matched Δ
and `C`. Check that all six established releases load with their published axis
digests. [Audit entry point](execution/audit_tenv_2k.py) · [Receipt](log/366-two-tenv-axes.json) · [Axis decision](../../docs/adr/0101-carps-vae-prediction-axes.md)

## Read accepted inputs in bounded windows

Add `cpst/03_data_train/output/api` to `PYTHONPATH` and open releases as context
managers. Use `tenv_axes.open_tenv_axis_training(axis_name, root=None)` for any
of the eight axes. Its `root` is a Step 03 directory. Use
`training_data.open_tenv_training(root=None)` for an explicit release directory
or the default `hvg-5k`. Preserve the six-axis lookup in `training_data` for
sealed consumers. Read controls and treatments or Δ in windows of at most 256
rows, and resolve each treatment's matched control with `control_row`.
[Axis lookup](output/api/tenv_axes.py) · [Matrix reader](output/api/training_data.py)

For Tahoe, follow the [03.2 reader setup](03.2_tahoe_hvg-5k/AGENTS.md#read-and-verify-the-accepted-release)
and use `tahoe_expression.open_tahoe(path=None)` with `read` or `batches` for
`X`, `control` or `delta` in windows of at most 256 rows. Pass the materialized
release directory when this checkout contains the tracked manifest alone.
[Tahoe reader](03.2_tahoe_hvg-5k/output/api/tahoe_expression.py)
