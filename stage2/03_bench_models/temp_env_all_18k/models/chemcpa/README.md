# chemCPA adapter — run card (ticket 13)

**Status:** five-fold benchmark complete (2026-10-21/22): AUPRC k=200 **0.0608** (95% CI 0.0578–0.0640), rank 10/12, `ideal_fit: partial`. It is `partial` because the 978-gene L1000 pretrain was omitted. The fold-0 smoke below is a wiring proof, not a ranking.

- **Source:** `pert_model/chemCPA/source/` (vendored, MIT). No file under
  `pert_model/` is modified or generated.
- **Environment:** conda `arc-stage2-bench-chemcpa` (Python 3.9, torch
  `2.8.0+cu128`, Lightning 2.2.5, RDKit 2025.09.2), cloned from the setup
  bringup env (already CUDA; no rebuild needed).
- **Config:** [`config.json`](config.json) records the paper's hyperparameters,
  the three-optimizer schedule, and the 18,164-gene from-scratch setting.
- **Code:** `stage2/03_bench_models/temp_env_all_18k/scripts/adapters/chemcpa/run.py`.
- **Fidelity:** [`FIDELITY.md`](FIDELITY.md).

## Command

```bash
# smoke (one epoch, fold 0)
BENCH_SMOKE=1 bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model chemcpa --fold 0 --timeout 5400

# full recipe (ticket 16)
bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model chemcpa --fold 0
```

Artifacts land in `stage2/03_bench_models/temp_env_all_18k/output/chemcpa/fold_0/`.

## Smoke result (fold 0, 1 epoch)

| quantity | value |
|---|---|
| fit wall clock | 7.7 s |
| predict wall clock (valid / test) | 3.1 s / 2.7 s |
| peak RSS | 2.0 GB |
| valid Gaussian NLL | −0.2729 |
| test AUPRC k=200 (headline) | 0.0814 |
| test AUPRC k=100 / k=500 | 0.0594 / 0.1207 |
| test PCC / Spearman | 0.0665 / 0.2314 |
| test R² / RMSE | −7.35 / 0.3392 |
| checkpoint | 14.6 M params, 18,164-gene output |

All 837 test wells produced finite vectors and were averaged to 276 groups.
The model is intentionally un-pretrained (the 978-gene L1000 stage is
quarantined), so the smoke is only a wiring proof.

## Why the adapter drives `ComPert` directly

The shipped Hydra pipeline cannot run a new dataset as shipped: its validation
dataloader returns a dict that Lightning 2.x rejects, the collate moves every
batch to `cuda` unconditionally, and the dataset requires DEG bookkeeping that
does not exist for this corpus. The adapter therefore constructs the authors'
`ComPert` and reproduces their training step: three Adam optimizers, StepLR(6,
0.9), Gaussian NLL, adversarial drug/covariate terms with the published
weights, gradient penalty, and the `global_step % adversary_steps` alternation.

- Molecule encoder G: the author's `rdkit` embedding is RDKit2D descriptors,
  standardised and variance-filtered. Descriptor mean/std and the
  low-variance filter are fitted on **training drugs only** (stricter than the
  author's whole-dataset fit).
- Dose: min-max over training doses; controls sit at 0 (the authors' vehicle
  = 1.0 hack is not copied — this corpus has no vehicle-treated rows).
- Prediction: the matched plate DMSO mean is encoded as the basal state, the
  test drug/dose is added, and the DMSO mean is subtracted to give Δ.

<!-- ticket16-run -->
## Full five-fold run (ticket 16)

| fold | status | wall s | fit s | pred valid s | pred test s | peak GPU MB* | final epoch/step | best epoch/step |
|---|---|---|---|---|---|---|---|---|
| 0 | completed | 77 | 58.7 | 2.8 | 2.6 | 3847 | 201 | 199 |
| 1 | completed | 77 | 58.9 | 2.9 | 2.6 | 3827 | 201 | 197 |
| 2 | completed | 78 | 57.8 | 2.8 | 2.6 | 3847 | 201 | 201 |
| 3 | completed | 75 | 56.8 | 2.8 | 2.6 | 3739 | 201 | 200 |
| 4 | completed | 76 | 57.0 | 2.8 | 2.6 | 3706 | 201 | 199 |

*Total device memory including the ~2.5 GB desktop baseline; the 30 s sampler is `output/gpu_watch.csv`. Machine-readable copy: `results/run_cards.csv`.
