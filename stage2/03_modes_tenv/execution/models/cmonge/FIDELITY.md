# CMonge tenv adapter fidelity

Source: `pert_model/CMonge/source/`, commit `195f5cd1822e7c7de7c7ee842a7e4082992ccc04`.
Model-specific admission ticket: #50 (parent specification #40, discovery #41).

The adapter calls upstream `ConditionalPerturbationNetwork` unchanged, with its
four 64-unit GELU hidden layers and residual map. The native direct-expression
path has `dim_data=18164`; no autoencoder is inserted or trained separately.
The upstream optional AE preprocessing path is omitted. This is the native
full-dimensional transport architecture, not an HVG or latent reconstruction.

A randomly initialized `Dense(50)` maps binary Morgan-2048 into the existing
continuous chemical-condition interface. Its parameters and all native transport
parameters share one JAX gradient/AdamW update. Actual positive floating dose in
micromolar units enters as `log(dose)` appended to the adapted chemical vector;
this avoids upstream's `int(dose)` truncation. Native two-context boundaries are
`(0,50)` and `(0,51)`, retaining the upstream nested chemistry/chemistry+dose
context convention. Native `dim_cond_map=50` creates context maps of 50 and 1.

Training conditions are equally sampled training drug-dose groups. Each sample
uses its own measured same-plate DMSO mean as source and source + training Δ as
target. This replaces the author's distribution of individual control cells
with the benchmark's approved baseline input. Held-out drugs never contribute
training responses. All expression access is per-group mmap slicing.

The objective calls upstream `fitting_loss` (Sinkhorn divergence, Euclidean cost,
epsilon 1) plus upstream `regularizer` (Monge gap, Euclidean cost, epsilon .01),
with native unit regularizer strength. AdamW uses learning rate .0001, weight
decay .00001 and cosine decay to .01 of initial rate. Formal budget is 10000
iterations, batch size 512; the single OOM retry may halve the batch. Validation
selects the minimum of the same author objective on one randomly sampled held-out
condition every 100 iterations. An independent per-run seeded RNG samples both
the condition and its batch; sampled group IDs are saved in checkpoint history.
A final off-cadence iteration is also validated to support bounded qualifiers.

Qualification uses 3 iterations, batch 4 and at most 3 groups on the complete gene
axis. It can establish train/save/reload/predict and adapter/baseline/dose behavior,
not formal performance or full-budget convergence. Checkpoint evidence reports
parameter change, baseline/dose sensitivity, exact reload agreement and measured
host resources. Real full-axis qualification passed on the frozen tenv release; see admission.json.
The independent review and complete formal nine-run matrix remain separate.

No historical benchmark implementation or upstream source file is modified or
imported. The existing CPU JAX environment is preserved. An isolated matching CUDA12
plugin overlay enables the separately qualified GPU execution; see CUDA.md.
Timings and recipe deviations are explicit rather than claiming the native AE
benchmark was reproduced.

Native sampled-condition validation was restored after resume; the old all-group
profile and its unique OOM retry remain in `pause_record.json`. Fresh native
qualification/profile evidence is recorded in `admission.json`.

Final native batch512 GPU qualification passed 100 training iterations and one
seeded held-out validation event, full finite 1401 × 18164 prediction, exact reload
and joint adapter/baseline/dose behavior. Measured training 72.350 s projects 2.010 h
for the unchanged 10,000-step recipe; formal allowance remains four hours. GPU
allocator/pool peaks and conservative workstation telemetry are in admission.json.
The first CUDA attempt required one assembler-path repair; the CPU-only batch-512 OOM
is preserved separately and does not establish a GPU batch 512 failure. No extra OOM
retry was allowed during the successful GPU batch 512 profile.
