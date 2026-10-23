# chemCPA fidelity ledger (ticket 13)

> Derived from `stage2/03_bench_models/temp_env_all_18k/models/chemcpa/FIDELITY.md`; the gene-space rows are re-declared for the 5,000-HVG sensitivity axis (ADR 0028).

| component | author setting | run setting | class | reason | expected effect | ticket |
|---|---|---|---|---|---|---|
| gene space | 977 LINCS -> 2,000 sci-Plex fine-tune | 5,000 rank-sum HVG from scratch | recipe_substitution | fixed benchmark gene universe (ADR 0028) | negative: 2.5x the widest author gene set | 09 |
| corpus scale | ~1.3 M LINCS bulk + sci-Plex | 3,600 wells/fold, 276 train drugs | recipe_substitution | curated CM view is the benchmark | large negative | 09 |
| initialisation | paper models start from the L1000 pretrain | random init | recipe_omission | ticket 09 quarantines 978-landmark artifacts | **largest negative**: the paper's 0.69 -> 0.46 gain is lost | 09 |
| 978-gene pretrain | 20 h on L1000, then fine-tune | skipped | recipe_omission | 978-landmark training is quarantined (ADR 0020/0024) | large negative; do not read a low chemCPA row as "chemCPA is weak" | 09 |
| molecule encoder G | pretrained frozen RDKit2D/chemVAE parquet | frozen RDKit2D descriptors, mean/std and variance filter fitted on train drugs | recipe_substitution | the shipped parquet does not contain this corpus; train-only fit avoids transductive structure statistics | small; cold-drug quality depends on chemical neighbourhood | 09 |
| entrypoint | `train_hydra.py` + Lightning data module | direct `ComPert` + author optimizer/loss loop | adapter_extension | the shipped pipeline cannot run a new dataset (dict val loader, unconditional CUDA collate, DEG bookkeeping) | none if the step is reproduced exactly; it is | 01 / 13 |
| decoder / architecture | dim 32, 256x4 AE, 128x4 drug MLP, 128x3 adversary, 64x3 doser | identical | config_only | paper hyperparameters | none | 09 |
| optimizers | 3 Adam 1.575e-3 / 1.192e-3 / 1.575e-3; StepLR(6, 0.9); grad clip 1.0 on AE+doser | identical | config_only | author defaults | none | 09 |
| adversarial loss | Gaussian NLL - 9.10*drug - 16.17*cov, penalty 0.455, every 3 steps | identical | config_only | author defaults | none | 09 |
| dose | min-max over the dataset; vehicle forced to 1.0 | min-max over training doses; controls at 0 | adapter_extension | no vehicle-treated rows exist; the vehicle hack would put DMSO at max dose | small; removes a known author quirk | 09 |
| epoch selection | validate every `checkpoint_freq=50` | every epoch, early stop patience 50 | config_only | one-shot cold-drug ranking; cheap at this size | none expected | 09 |
| prediction | decode treated state from a basal control | basal = matched plate DMSO mean; subtract it for Δ | config_only | benchmark Δ contract | none | 04 / 09 |
