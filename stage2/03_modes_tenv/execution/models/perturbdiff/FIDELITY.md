# PerturbDiff tenv fidelity boundary

Pinned source: `pert_model/PerturbDiff/source/`, commit
`31be89cb8cce5dddef0782bac29d5c748ce55931` (MIT); source manifest remains authoritative.
The adapter calls native `PlModel._compute_loss`, `GaussianDiffusion.training_losses`
and native DDIM sampling. It retains Cross_DiT with 8 MM_DiT blocks, width512,
8 heads, native self-conditioning, ReLU output, categorical dose encoder,
chemical projection and conditional/control branches. Diffusion training has1000
steps. Author `configs/rawdata_diffusion_sampling.yaml` defaults to **100 DDIM
reverse steps**, eta0, guidance1 and clipping, not1000 inference steps.

The existing frozen-chemical-vector lookup is replaced by the common Morgan table,
`Linear(2048,768)` and native `Linear(768,128)` projection. The adapter is jointly
optimized, without chemical pretraining. 768 is the selected existing ChemBERTa
vector interface width; the native loader infers this dimension from its supplied
vector table rather than hard-coding it. Native dose embedding, constant one-celltype
embedding and native gather projection remain trainable. Dose vocabulary includes
all allowed case metadata; test-treated expression is absent. The internal dataset tag `tahoe_tenv` selects the native drug-data router (which
asserts on unknown dataset tags); it does not mean the data came from Tahoe.
Batch covariates are
omitted for this single cell model rather than exposing plate identity as a shortcut.

Inputs use the frozen same-plate DMSO mean and target log expression. Native
`normalize_counts=10` divides both by10 and predictions are multiplied by10 before
subtracting the original baseline. Native `p_drop_control=0` follows the published
scratch command; conditional dropout remains0.5. No upstream checkpoint is loaded.

Training preserves native MSE + energy MMD (factor1, blur0.05), timestep uniform
sampling, AdamW(lr0.0002, betas0.9/0.98, eps1e-8), weight decay0.01 with bias/norm
exemption, gradient norm cap1, warmup200 and cosine floor/plateau0.1. The formal
native scratch command ceiling is200,000 steps; the base trainer default1,560,300 is overridden by the README
`COMMON_TRAIN_RUNTIME trainer.max_steps=200000` used by the scratch Tahoe command;
a four-hour stop before the chosen200,000-step author recipe is an incomplete run,
not a successful completion. No early-stop rule or smaller backbone is invented.

Declared data/runtime deviations: bulk replicate wells form native populations
within each drug-dose group; training resamples three wells per population rather
than native single-cell256. Training batches contain8 populations and use2 CPU
threads. Validation uses up to5 batches (the author's documented runtime limit),
on the same author objective with a fixed noise/resampling seed for checkpoint
comparability, every training epoch (derived per fold and batch size). Test inference buckets groups by actual well count
and retains every well. No evaluator files are read. EMA uses the same README scratch runtime block: decay0.99, updates every10
steps and every epoch. Native `EMAWeightAveraging` swaps averaged parameters for
validation and the selected checkpoint; the jointly trained adapter is included.
The optimizer cosine schedule uses the same200,000-step max_steps override. An isolated Python overlay avoids modifying the shared arc env.

Compact tests deliberately use12 training steps,64 genes and2 reverse sampling
steps. These exercise artifact behavior only and cannot establish full-gene
qualification or formal training feasibility. Real qualification and throughput
must be measured separately before admission.
