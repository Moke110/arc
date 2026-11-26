# Compare TranSiGen measured encoder targets and signed Delta readouts

**Status:** accepted (2026-11-26).

In Step 06.3, the independent signed Delta readout had higher mean Delta PCC
and lower mean Delta RMSE than the original. Adding normalized FiLM did not
improve either mean relative to that readout alone. Register
`stage2/06_tsg_shape/06.4_round4` to compare the measured representation behind
the Gaussian alignment and whether the signed Delta decoder is shared. Use the
original non-FiLM stochastic Gaussian prediction head throughout.

Let `X1` be measured same-plate control, `X2` measured treated expression, and
`Δ = X2 − X1`. `K1`, `K2`, and `KΔ` are the control, measured-treated, and
measured-Delta Gaussian KLs against the standard normal where those encoders
exist. `qpred` is the Gaussian predicted from control latent, drug and dose.
`MSE_Δpred` and `L1_Δpred` compare signed predicted and measured Delta.

| Arm | Measured encoders and readout | Objective |
| --- | --- | --- |
| `XX-XX` | Round 3 original. Encode `X1` and `X2`. Decode the predicted latent to `X2_hat` with the measured-treated decoder and derive the signed Delta as `X2_hat − X1`. | `0.1 (K1 + K2 + KL(qpred ‖ q2)) + MSE_X1 + MSE_X2 + MSE_Δpred + 1.75 L1_Δpred` |
| `XX-XXΔ` | Round 3 independent signed Delta arm. Encode and reconstruct `X1` and `X2`. Decode the predicted latent to signed Delta with an independent decoder. | `0.1 (K1 + K2 + KL(qpred ‖ q2)) + MSE_X1 + MSE_X2 + MSE_Δpred + 1.75 L1_Δpred` |
| `X-XΔ` | Encode and reconstruct `X1` only. Decode the predicted latent to signed Delta. Omit the measured-treated encoder and decoder reconstruction path, `K2`, and all KL or prior penalties on the predicted latent. | `0.1 K1 + MSE_X1 + 2 MSE_Δpred + 1.75 L1_Δpred` |
| `XΔ-XΔ` | Encode and reconstruct `X1` and measured signed Delta. Share one signed Delta decoder between the measured Delta latent and the predicted latent. | `0.1 (K1 + KΔ + KL(qpred ‖ qΔ)) + MSE_X1 + MSE_Δrec + MSE_Δpred + 1.75 L1_Δpred` |
| `XΔ-XΔΔ` | Encode and reconstruct `X1` and measured signed Delta. Use a separate signed Delta decoder for the predicted latent. | `0.1 (K1 + KΔ + KL(qpred ‖ qΔ)) + MSE_X1 + MSE_Δrec + MSE_Δpred + 1.75 L1_Δpred` |

`MSE_Δrec` compares the decoded measurement latent with measured signed Delta.
For direct Delta arms, predicted treated expression is `X1 + Δ_hat`. Prediction
in every arm uses measured `X1` and frozen drug and dose features. `X2` and
measured Delta provide training and validation targets.

Keep the [Round 3 protocol](0097-transigen-shape-round3.md) on the frozen tenv
cold-drug folds with Denoise-DEG 5k, same-plate measured DMSO, frozen KPGT,
training-fold dose transform, optimizer, stopping settings and evaluation
sampling seeds. Use training seeds `17`, `42`, and `2026` in each of the three
folds. Verify the 18 accepted fits for `XX-XX` and `XX-XXΔ` through the Round 3
release and their original lineage. Train and independently replay the 27
remaining fits. Score all 45 coordinates with the nine inherited metrics.

Report the paired contrasts `XΔ-XΔΔ − XX-XXΔ` and `XΔ-XΔ − XΔ-XΔΔ` on the
same fold and training seed. Publish matched coordinate values and foldwise
summaries for Delta PCC, Delta RMSE, Delta Spearman, and Delta E-distance.
Four-metric summaries weight all nine fold and seed coordinates equally. The
three folds define the held-out splits, with three repeated fits within each.
