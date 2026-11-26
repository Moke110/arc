# Combine TranSiGen direct Delta readouts with FiLM

**Status:** accepted and revised (2026-11-26).

In Step 06.2, the direct Delta readouts reduced mean Delta RMSE and had
lower mean Delta PCC than the original. FiLM improved both primary means.
Register `stage2/06_tsg_shape/06.3_round3` to compare the two combinations
with their component arms and the original reference.

| Arm | Initially registered configuration |
| --- | --- |
| `0_original` | Round 2 `original` reference. |
| `1_output-delta` | Round 2 `output-delta`, with an independent signed Delta decoder and the Step 05.2 objective. |
| `2_output-delta-no-x2` | Round 2 `output-delta-nox2`, with the treated decoder producing signed Delta. Keep the measured treated encoder and all three KL terms. Use `0.1 K + MSE_X1 + 2 MSE_Delta + 1.75 L1_Delta`. |
| `3_z-interact_FiLM` | Round 2 `z-interact_FiLM`, including identity modulation at initialization and the Step 05.2 objective. |
| `4_1+3` | Arm 1 readout and objective with arm 3's FiLM prediction head. |
| `5_2+3` | Arm 2 readout and objective with arm 3's FiLM prediction head. |

Use the Round 2 protocol, including its three frozen tenv cold-drug folds,
Denoise-DEG 5k, measured same-plate DMSO, frozen KPGT inputs, training-fold dose
transform, optimizer and checkpoint selection settings, and evaluation sampling
seeds. Both direct Delta readouts form predicted treated expression from the
measured control plus predicted Delta. Train both combinations from scratch.

Use training seeds `17`, `42`, `2026`, and `43` on each fold. The fourth seed
is specific to this study and has precedent in ADR 0082. The first attempt
registered 12 seed-42 reuses from arms 0 through 3 and 60 new coordinates.
The three originals retain their sealed Step 06.1 lineage.

Keep the nine Round 2 metrics, with Delta PCC and Delta RMSE as the primary
outcomes. Report scores and signed differences from `0_original` paired within
each `(fold, seed)`. For the primary outcomes, also report `4−1`, `4−3`, `5−2`,
and `5−3`, and the descriptive interaction contrasts `4−1−3+0` and `5−2−3+0`.
Summarize the equally weighted coordinates and show their variability. Training
seeds repeat fitting within the same three held-out folds.

## Revision after the first formal attempt

Eight new coordinates in arms 4 and 5 stopped at a nonfinite training loss.
The first attempt independently replayed 52 successful fits. The failure and
its diagnostics remain in the [blocked receipt](../../stage2_CarPS-VAE/06_tsg_shape/06.3_round3/log/356-blocked.json).

Add `LayerNorm(100, elementwise_affine=False)` after FiLM modulation and before
the Gaussian projection in arms 3, 4 and 5. The conditioning network retains
its identity modulation at initialization before normalization. This revision
uses the same population, optimizer, checkpoint selection, objectives,
training seeds and validation and test sampling seeds.

Train and independently replay all 36 coordinates for arms 3, 4 and 5,
including every seed-42 coordinate. The other 36 coordinates comprise nine
accepted seed-42 fits from the earlier studies and 27 independently replayed
fits from unchanged arms 0 through 2 in the first Round 3 attempt. Verify
their exact artifacts, registered protocol, population and source lineage
before reuse. Score all 72 coordinates together, preserving the paired
contrasts and four-seed summary registered above.

[Round 2 protocol](0096-transigen-shape-round2.md) ·
[Prior use of seed 43](0082-view-variants-tenv-tahoe-replicate-structure.md)
