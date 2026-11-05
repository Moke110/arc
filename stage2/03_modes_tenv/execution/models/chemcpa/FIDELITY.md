# chemCPA — native full-axis qualification adapter

Governing spec: #40; model qualification: #49. Upstream:
`pert_model/chemCPA/source/`, pinned
`43e830eb0958c54e4aa64442c17ec0fed19b3f15`; formal publication:
[NeurIPS 2022](https://neurips.cc/virtual/2022/poster/53227).
Interpreter: `/home/chang/miniconda3/envs/arc-stage2-bench-chemcpa/bin/python`;
package snapshot: `env.actual.lock`. Native setup evidence remains in
`stage2/01_setup/execution/chemcpa/`; its 977-gene checkpoint forward is separate
from this study's full-axis runtime gate. No public model weights are loaded.

## Preserved model and author optimization

The adapter imports immutable native `chemCPA.model.ComPert`. Its native basal
expression encoder, decoder, chemical embedding MLP, drug adversary and
amortized doser remain intact. It has 18,164 input genes and 36,328 native
mean/variance outputs. The expression output is the native nonnegative mean;
variance is native softplus. Prediction feeds the approved same-plate DMSO mean
into the encoder, then subtracts this baseline from the native mean output.

The source setup's `output/run.log` records a released RDKit checkpoint chemical
table of `(188,194)`. The approved added module is `Linear(2048,194)`; its output
is supplied through the existing embedding-table `.weight` interface. This table
contains one common Morgan vector per drug, not independent learned drug rows.
The native chemical MLP and amortized chemical-plus-dose network consume the
adapted vectors. The adapter is included in the autoencoder optimizer and is
randomly initialized jointly with every predictor parameter. It adds 397,506
parameters. No independent chemical learning or external pretraining occurs.

`config.json` copies the fixed author `config/model/hparams/default.yaml` recipe,
not a random `set_hparams_` search: latent32, AE width256/depth4, chemical MLP
width128/depth4, adversary width128/depth3, doser width64/depth3, batch256.
Maximum201epochs comes from `config/training/default.yaml`. Three Adam optimizers
retain author learning rates/weight decays, StepLR every6epochs with gamma.9,
and native norm1 clipping for autoencoder and doser updates.

Training follows the author's self-reconstruction method: allowed training Δ
plus its approved baseline reconstructs treated expression; this expression is
the native encoder input and Gaussian reconstruction target. This is not a new
paired baseline-to-treated training objective. Counterfactual prediction encodes
the measured baseline, as the native evaluation code does. On minibatches
whose index modulo3 is0, drug BCE plus the native latent-gradient squared
penalty updates the adversary; other minibatches update autoencoder, native
chemical network, joint adapter and doser with
`GaussianNLL - 9.100951626404369 * drug_BCE`. Native penalty scale is
`0.4550475813202185`. All components and update counts are persisted.

## Explicit adaptations and limits

- The native `lightning_module.validation_step` is a placeholder, while its
  validation-end hook reads test-treated responses and R². This hook is never
  called. To satisfy ADR0042, each epoch evaluates the same native combined AE
  objective on validation responses, and minimum loss selects the checkpoint.
  The fixed nominal patience50 is made effective at this epoch cadence.
  This differs from the upstream checkpoint-frequency50/R² hook and is disclosed.
- No variable cell-type covariate is supplied by the frozen tenv case contract;
  `num_covariates=[0]` is a native option. There is no invented plate embedding or
  artificial covariate adversary. The measured plate control remains the
  prediction baseline input.
- The adversary output dimension covers the frozen drug metadata table;
  held-out drugs have no positive training responses. Validation BCE on unseen
  drugs can be difficult to interpret, so reconstruction and adversary losses
  remain separately visible. No test metric chooses a checkpoint.
- The source MLP ignores its configured `dropout` hyperparameter; this adapter
  does not invent missing dropout layers. Native BatchNorm is retained and
  switched to eval mode for validation/inference. A singleton final training
  batch joins the preceding minibatch.
- The source Lightning wrapper calls this alternating scheme through
  `global_step`; the adapter counts minibatches explicitly so the configured
  `adversary_steps=3` has a reproducible meaning with multiple optimizers.
- float32 mmap batches replace AnnData bulk loading. All 18,164 genes participate;
  there is no HVG bottleneck, appended pretrained-gene bridge, response lookup,
  or output padding. Only allowed predictor files are opened.
- SIGTERM/SIGINT and batch-checked training/inference budgets leave partial
  evidence and an interrupted checkpoint, never a completed prediction.
  Scheduler cumulative four-hour/retry accounting remains authoritative.

The full-axis subprocess test is intentionally separate from real-data
admission. It checks chemical/native/doser/adversary parameter changes, both
optimization branches and a positive gradient penalty, finite full-axis cold
predictions, baseline/dose sensitivity and exact checkpoint reload. An actual
tenv pilot and scheduler resource evidence are required before formal admission.

The synthetic test passed (8training/2validation/2test cases, 18,164 genes,
2epochs, batch2, CPU) in 4.19s wall time at 1,359,644KiB peak RSS under the
10G RAM/8G swap scope; see `synthetic_test.log` and `synthetic_test_record.json`.
The scheduler's `TENV_BATCH_SCALE` is honored for the single permitted OOM
retry; requested/effective batch sizes and scale are recorded in checkpoint config.

Reproduce the synthetic seam test under the repository memory cap:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
  /usr/bin/time -v /home/chang/miniconda3/envs/arc-stage2-bench-chemcpa/bin/python \
  -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_chemcpa.py -v
```

Formal command: the isolated interpreter plus `execution/script/chemcpa_adapter.py`.
The artifact CLI appends `--request <request.json>`; pilot-only recipe overrides
use an explicit `--config <json>` included in the command identity.
