# Compare CarPS-VAE dose inputs on cold-dose-drug folds

**Status** Accepted design, 2026-11-27.

Create Stage 2 Step `10_cpsv_dose` to compare dose inputs for simultaneous
generalization to unseen drug parents and unseen measured doses. Use only the
three frozen `cold_dose_drug` folds on the [tenv Denoise-DEG 5k release](../../cpst/03_data_train/03.7_tenv_denoise-5k/README.md)
and training seeds 17, 42 and 2026. Test the held-parent and held-dose
intersection. Use other parents at the other doses for training and validation.
Exclude conditions that match only one holdout. The held doses are 0.2, 1 and
10 µM. Each dose is tied to its fold and drug population.

Use the same Denoise-DEG 5k axis, treatment conditions and measured same-plate
DMSO controls in every arm. The nine CarPS-VAE arms share the `XΔ-XΔ` backbone,
selected loss and frozen UniMol. `Xmol` is the 512-dimensional frozen chemical
representation at the input to the chemical projection. `Zmol` is the
400-dimensional projected chemical representation supplied to the Gaussian
prediction head.

| Arm | Dose input or reference |
| --- | --- |
| `TranSiGen` | Step 10 model reference using the author `XX-XX` architecture and loss, frozen KPGT and no dose input. Predict treated expression from the measured same-plate control `X1` and score `Δpred = X2pred − X1` against `Δ = X2 − X1`. |
| `no-dose` | CarPS-VAE without dose. |
| `Xmol-scaling` | Multiply all `Xmol` coordinates by a common dose-dependent factor. |
| `Xmol-FiLM` | Apply dose-dependent feature-wise affine modulation to `Xmol`. |
| `Zmol-scaling` | Multiply all `Zmol` coordinates by a common dose-dependent factor. |
| `Zmol-FiLM` | Apply dose-dependent feature-wise affine modulation to `Zmol`. |
| `Xdose-num` | Append 32 copies of the log dose to `Xmol`. |
| `Xdose-MLP` | Append a learned 32-dimensional dose vector to `Xmol`. |
| `Zdose-num` | Append 32 copies of the log dose to the joint control and chemical latent vector. |
| `Zdose-MLP` | Append a learned 32-dimensional dose vector to the same joint latent vector. |

For the concatenation and FiLM arms use `s = log10(dose / 1 µM)` directly.
For both scaling arms multiply by `log10(1 + dose / 1 µM)`. Both dose embedding
MLPs use `1 → 64 → 32`, with hidden ReLU and linear output. Both FiLM arms use
`1 → 64 → 2D` with hidden ReLU and zero-initialized linear output. Split the
output into `D`-dimensional `γ` and `β` and apply `(1 + γ) ⊙ v + β`. The output
width is 1,024 for `Xmol-FiLM` and 800 for `Zmol-FiLM`. Fit all 90 coordinates
on their designated folds.

Select each checkpoint by minimum validation objective for its model. Publish
full-axis held-out results with Δ PCC and Δ RMSE as the primary metrics for the
CarPS-VAE dose comparison. Report the individual folds and pair CarPS-VAE
methods within fold and training seed against `no-dose`. Report the reference
as `TranSiGen` in the figures. The study publishes comparison evidence.
[ADR 0104](0104-carps-vae-no-dose.md) records the downstream default.

## See

- [ADR 0099](0099-carps-vae-alpha-and-chemical-representations.md) for the fixed CarPS-VAE backbone.
- [ADR 0100](0100-carps-vae-unimol.md) for UniMol.
- [ADR 0102](0102-stage2-axis-selection-and-comparison-figures.md) for the selected axis and comparison figures.
