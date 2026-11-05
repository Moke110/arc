# PRnet tenv fidelity and qualification

Source: `pert_model/PRnet/source/`, pinned
`f19174bde3ed2633f54c7831799cc38c4ffc7a0d`; publication
[Qi et al., Nature Communications 2024](https://doi.org/10.1038/s41467-024-53457-1).
The adapter imports the unchanged native `models.PRnet.PGM`. Environment:
`/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python`, package record
`env.actual.lock`. No external checkpoint or pretrained drug encoder is loaded.

## Architecture, inputs and objective

Each run initializes `Linear(2048,1024)` and the PGM randomly from its training
seed. The common binary Morgan vector is projected, multiplied by the native
`log10(dose_um+1)` scaling, then passed into the native `CombAdaptor`. The measured
same-plate baseline enters the PGM encoder; native independent 10-dimensional
normal noise enters the decoder. PGM encoder, latent, decoder and adaptor are
unchanged, including the author trainer's normal(1,.02) BatchNorm scale initialization. The expression axis is the complete frozen ordered 18,164 genes.

The request's training/validation arrays are Δ; the minibatch loader adds each
approved baseline to recover the expression target. Native positive expression
mean and softplus variance feed `torch.nn.GaussianNLLLoss`, the author's `GUSS`
objective. The default architecture and optimization follow `train_sciplex.py`:
hidden `[128]`, z=64, chemical condition=64, dropout=.05, batch512, 100epochs,
Adam lr=.001, weight decay=1e-8; ReduceLROnPlateau factor=.5, patience5, min_lr1e-8.
The joint Morgan adapter adds 2,098,176 trainable parameters. No HPO is performed.

## Required deviations and explicit conventions

- The agreed Morgan representation replaces native feature-based rFCFP1024.
  The linear projection is trained jointly with all native predictor parameters.
- The upstream trainer optimizes Gaussian NLL but selects/schedules on sampled
  validation MSE (`PRnetTrainer.py:351`). ADR0042 requires validation to match
  the author's training loss, so this adapter selects/schedules on Gaussian NLL.
  It uses eval mode and repeats one fixed validation noise stream each epoch.
- Early stopping uses the accepted nominal 20 consecutive nonimproving epochs.
  Upstream's `patient <= 20` then `else` actually stops on the 22nd failed epoch;
  the adapter makes its 20-epoch policy explicit in `config.json`.
- Prediction uses the native conditional Gaussian mean for one independently
  seeded normal-noise stream, then subtracts the same approved baseline.
  It does not add observation sampling noise from the learned variance.
- Source AnnData/DataLoader plumbing is replaced by request allowlists and
  float32 mmap minibatches. No entire well×gene dense copy is made. A final
  singleton training batch is joined to the preceding batch for native BatchNorm.
- Validation and inference do not update BatchNorm/dropout. Small fixed probe
  checks compare the chosen checkpoint with its reloaded forward, baseline
  sensitivity and adapter/native parameter changes; no test-treated expression
  is opened. The probe uses only allowed test chemical/dose/control inputs.
- Four-hour training and bounded inference limits are checked between batches;
  scheduler cumulative budgeting/retries remain authoritative. SIGTERM/SIGINT
  request a safe stop with `interrupted.pt`, not completed prediction status.
  Retry initialization is fresh; the adapter does not silently reuse a stale
  partial state. Completed run resume is validated by the artifact CLI.

## Evidence so far

TDD subprocess test first failed for missing config/adapter, then passed on a
compact full-axis fixture (6train/2validation/2test, 2epochs, batch2, CPU). It
asserts finite 2×18,164 predictions, exact case order, matching request identity,
nonzero parameter changes for both modules, baseline sensitivity, identical
selected-checkpoint reload, Gaussian validation objective and blob checksums.
The capped run took 2.15s wall time and peaked at 1,036,336KiB RSS under
`MemoryMax=10G`, `MemorySwapMax=8G`; see `synthetic_test_record.json`.

The separately scheduled real-tenv pilot also passed: fold 0 / seed 17, two
epochs, all 18,164 genes and 1,401 finite test-well predictions. Exact reload,
trained adapter/predictor changes and baseline sensitivity passed. Peak RSS was
1,423,696 KiB; GPU reservation peaked at 864,026,624 bytes. See `admission.json`
and `pilot_record/`. On explicit resume, `resume_audit.json` records revalidation
of current identities, all input hashes, prediction and checkpoint hashes, and
the actual environment lock. This qualifies execution, not convergence or
formal test performance.

Run the synthetic seam test from the repository root:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
  /usr/bin/time -v /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python \
  -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_prnet.py -v
```

The artifact command is the isolated interpreter plus
`execution/script/prnet_adapter.py --request <request.json>`. An optional
`--config <json>` is hashed by the benchmark command identity and is required
for the explicitly labelled two-epoch pilot override. Formal runs use the
default `config.json` unchanged.
