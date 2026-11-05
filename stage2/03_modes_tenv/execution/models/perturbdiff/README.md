# PerturbDiff qualification

Ticket #51, part of #40. **Qualified** for bounded full18,164-gene execution and
measured formal-recipe feasibility. The native setup is separately recorded in
`stage2/01_setup/execution/perturbdiff/`; source and deviations are in `FIDELITY.md`.
The effective isolated environment is pinned without changing shared arc.

The current fresh profile is `gpu_native_profile.queue.json`, with output under
`execution/pilots/perturbdiff_gpu_native_checkpoint_profile/`. It preserves native
8×512 Cross_DiT,1,000 training diffusion steps,100 reverse DDIM steps, native
MSE+energyMMD and GPU EMA0.99/every10steps+epoch. All1,401 well predictions over
18,164genes are finite and correctly ordered, exact reload passes, and joint
adapter/predictor updates plus baseline/dose output participation are nonzero.
The100-step qualifier uses checkpoint interval100; formal `config.json` uses the
author scratch200,000 steps and native checkpoint interval20,000. Every native
epoch validates; off-epoch eligible checkpoints get fresh validation. Selection
uses minimum fresh author validation loss among eligible checkpoint states.

**Measured current resource profile:**100steps plus validation/checkpoint8.385s,
inference88.799s, time-v RSS2,501,780KiB, conservative process tree2.461GiB,
CUDA allocated5,156,252,160bytes/reserved5,305,794,560bytes. Predictor, Lightning
metadata and EMA parameter devices all report`cuda:0`.

Steady training steps average0.0574673s. The200,000-step projection includes
epoch validation, native epoch EMA, up to10checkpoint writes and progress I/O:
central3.294h/run, conservative observed-overhead estimate3.314h, or3.645h with
a disclosed10% runtime planning cushion. Nine serial runs plus measured inference
project roughly30.05h. These are100-step extrapolations, not completed formal
fits. The hard4h cumulative training ceiling remains enforced; a budget-stop is
incomplete and cannot become successful truncated training. No test metrics or
convergence claim is part of qualification.

The old CPU-EMA profile and its5.533h estimate are explicitly **superseded** in
`resource_assessment_cpu_ema_superseded.json`: parent `nn.Module.to(cuda)` moved
native tensors but bypassed the Lightning child `to()` metadata update, causing
native EMA to be allocated on CPU. The adapter now explicitly calls native.to,
asserts device agreement and passes the GPU artifact regression. Native EMA math,
decay and schedule are unchanged. The first GPU profile also used unnecessarily
dense checkpoint writes; native20k cadence is restored and independently tested.

`admission.json`, `resource_assessment.json`, `gpu_native_profile_*.json` and the
telemetry/time-v log preserve current proof. Large checkpoints and prediction
arrays stay local-only. Earlier `pilot_*` and `profile_*` records are historical
setup/diagnostic evidence, not formal predictions or current admission identities.
