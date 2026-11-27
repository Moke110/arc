# Step 07 CarPS-Transcriptome alpha

Read the release report and the decision that freezes the Step 06.4 `XΔ-XΔ`
arm as alpha. For model, fit and replay changes, use the shared training and
evaluation APIs. [Release report](README.md) · [ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) · [Shared APIs](../util/README.md#transigen-training-and-evaluation)

## Verify the alpha release

1. Run from the materialized repository root in the `arc` environment with CUDA
   and the pinned Git revisions available. Use the sealed Step 06.4 formal
   matrix and report, the ordered tenv Denoise-DEG 5k training release, the frozen
   tenv KPGT vectors at `cpst/util/output/representations`, and the author
   TranSiGen source at `pert_model/TranSiGen/source`. Keep the Step 06.4
   checkpoint directories read-only. [Formal receipt](../06_tsg_shape/06.4_round4/log/359-formal.json) · [Training manifest](../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) · [Representation reader and manifest pins](../util/output/api/cpst_util/representations.py) · [Author source manifest](../../pert_model/TranSiGen/manifest.sha256)
2. Check available memory before launching the capped verifier. Follow the
   execution limits and sample GPU, WSL and Windows host memory every ten
   seconds with the coordinator sampler. [Execution limits](../../docs/agents/execution-environment.md) · [Coordinator sampler](../06_tsg_shape/06.1_round1/execution/sample_resources.py)

   ```bash
   free -g
   systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
     env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
     python cpst/07_cpst_alpha/execution/verify_alpha.py
   ```

   For identical inputs at alternate locations, append `--round-four`,
   `--training`, `--representations` or `--source` with the matching directory.
   For standalone replay of a registration through the original Step 06.4 runner,
   follow its pinned-checkout procedure. [Alpha verifier](execution/verify_alpha.py) · [Step 06.4 guide](../06_tsg_shape/06.4_round4/AGENTS.md)
3. Require the verifier's `PASS` before accepting alpha or a study that depends
   on it. It checks the pinned Git source, the sealed ledgers and all nine
   fold–seed coordinates. Validation predictions and objective must match
   exactly at sampling seed `22617`, and held-out expression and signed Δ arrays
   must match exactly at `100042`. The full result must equal the accepted
   receipt. [Release manifest](output/release.json) · [Parity receipt](log/362-alpha.json) · [Resource record](log/362-resources.json)

## Call the fixed model

1. Add `cpst/07_cpst_alpha/output/api` to the Python import path and import
   `AlphaAdapter` from `cpst_alpha`. Pass the author source directory and the
   dose transform fitted on the training fold. Call
   `adapter.build(training_seed, device)` to initialize the model. For inference
   from a fitted model, restore the corresponding checkpoint and set the
   prediction sampling seed using the study's replay path. [Public module](output/api/cpst_alpha.py) · [Shared checkpoint and inference code](../util/training/transigen.py)
2. Call `adapter.predict(model, inputs, rows)` with measured same-plate DMSO
   controls, raw positive doses and frozen chemical vectors in the registered
   condition order. The adapter applies the training-fold log-dose transform
   and returns sampled treated expression in float32 on the ordered gene
   axis. The alpha release uses 5,000 genes and predicts treated expression
   by adding the decoded signed Δ to the measured control. [Step 04 input API](../04_tsg_axis/output/api/tsg_axis/inputs.py) · [Shared native adapter](../util/training/native.py)
3. Use the native 2,304-input KPGT path for the alpha default. A successor
   supplies its own row-aligned vectors and `chemical_width` for other
   representations. To customize the affine projection, pass a
   `projection_factory` that accepts `(native_width, 400)` and returns a CPU
   `torch.nn.Linear` of that shape. The adapter replaces the first projection
   under isolated RNG state. Keep the variant's representation, fit registration
   and evidence with its owning experiment. [Alpha implementation](../util/training/cpst_model.py) · [Step 08 guide](../08_cpst_chem-emb/AGENTS.md)

For constructor or projection changes, run the focused test in `arc` with
the author source materialized. It checks native state and prediction parity
plus preservation of downstream initialization and RNG state for a variant.

```bash
python -m unittest discover -s cpst/07_cpst_alpha/execution/tests -p test_alpha.py
```
