# CPA whole-gene Morgan adapter

The immutable source is `pert_model/CPA/source/`, upstream commit
`fbd7c0250edc23eff003a10c99655579c53afd63`, BSD-3-Clause. The isolated interpreter
is `arc-stage2-bench-cpa`; package pins are recorded in `stage2/_envs/cpa.lock`.
Native setup evidence is `stage2/01_setup/execution/cpa/README.md` and establishes
checkpoint loading on the Kang data, separately from this 18,164-gene admission.

The adapter invokes native `CPAModule` with Gaussian reconstruction and native
`CPATrainingPlan.training_step`, retaining its reconstruction NLL, adversarial
cross entropy, gradient penalty, alternating optimizer schedule and learning-rate
schedule. Inputs are streamed from `.npy` memmaps; training reconstruction uses
observed treated log-expression and approved measured controls. The predictor
encodes same-plate baseline expression and adds the native perturbation latent;
its Gaussian mean is converted back to Δ. Test-treated expression is absent from
the request.

The author's existing `drug_embeddings` interface receives the frozen common
Morgan table. The native `pert_transformation = Linear(2048, n_latent)` is the
agreed joint adapter. Native `configure_optimizers()` omits this projection from
its optimizer parameter lists: this local harness explicitly adds its parameters
to the autoencoder optimizer and verifies a nonzero trained change. Upstream
source bytes are unchanged. This is a disclosed harness correction.

The native `doser_type='linear'` option passes actual µM dose without fitting a
per-drug dose scaler. This replaces the usual learned logsigmoid option, whose
unseen-drug parameters would remain untrained. It is an existing native dose path,
not a replacement backbone. The fixed 18,164-gene axis, shared Morgan extraction,
streamed training host and single-context/no varying-covariate configuration are
the other recipe deviations. The 512×3 encoder/decoder, 128-dimensional latent,
2000-epoch cap, batch128, native three optimizers, 45-epoch StepLR cadence and
patience10 are fixed without test-score tuning. Two-epoch pilot overrides prove
execution and do not establish convergence.

`checkpoint.json` records selected validation reconstruction NLL, all deviations,
parameter counts, behavioral adapter/predictor/adversary changes, baseline
sensitivity, exact reload equality, timings and peaks. Model and optimizer
states are saved together in checksum-declared `best.pt`. Only successful full-axis
real train/save/reload/predict evidence establishes runtime admission; synthetic
fixtures alone do not establish it.

The scheduler-owned real-data pilot completed on fold0/seed17 using all
18,164 genes, 2,388 treated training wells plus 16 measured controls, 417
validation wells and 1,401 full-axis test predictions. Two complete epochs used
9.976 seconds training; inference/checkpoint reload used 0.687 seconds. Adapter,
predictor and adversary parameters changed, baseline mean sensitivity was
nonzero, and reload predictions agreed exactly. Peak host RSS was 1,876,412 KiB;
GPU allocation/reservation peaked at 1,359,595,008/1,503,657,984 bytes. The
checkpoint and telemetry evidence copies are under `pilot_record/`; full binary
artifacts remain local under `execution/pilots/cpa/`. `admission.json` and
`formal_profile.json` record this runtime admission and measured resource
profile. No test metrics were calculated at admission.
