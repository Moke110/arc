# TranSiGen architecture studies

Step 06 used four rounds to test how TranSiGen represents measured expression
and predicts signed perturbation responses. The studies shared the three frozen
tenv cold-drug folds, ordered Denoise-DEG 5k axis, measured same-plate DMSO,
frozen KPGT features and training-fold standardized log dose. They started
from the Step 05.2 objective. The first two rounds used seed 42, and the later
rounds repeated training within the same held-out folds.
[Training release](../03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Starting objective](../05_tsg_loss/05.2_mse_l1/README.md)

1. **06.1 Expression width and depth.** The combined report covers ten
   architectures. The width extension added six fits and reused the 24
   accepted fits. The 4,000-unit encoder had the highest mean Δ PCC (0.433229)
   and lowest mean Δ RMSE (0.248576). Adding a second hidden layer to both
   the encoder and decoder gave the lowest mean Δ E-distance and Δ MMD.
   Round 2 kept the original expression dimensions and reused the three
   original fits as paired references.
   [Round 1 results](06.1_round1/README.md#formal-results) ·
   [Combined acceptance](06.1_round1/log/350-formal.json)
2. **06.2 Prediction architectures.** Seven new architectures tested changes
   to the Δ readout, drug–dose conditioning, encoder sharing and prediction
   mean. FiLM without post-modulation normalization had the highest mean Δ
   PCC (0.439568). The independent signed Δ decoder had the lowest mean Δ
   RMSE (0.241952). Round 3 carried forward FiLM and both signed Δ readouts.
   [Round 2 results](06.2_round2/README.md#formal-results) ·
   [Formal acceptance](06.2_round2/log/353-formal.json)
3. **06.3 Signed Δ readouts and FiLM.** The initial formal attempt stopped
   after the prediction KL overflowed at eight fold–seed coordinates in the
   two unnormalized FiLM combinations (`4_1+3` and `5_2+3`). The accepted
   revision inserted parameter-free `LayerNorm(100)` after modulation and
   before the Gaussian projection in all three FiLM arms. Across six arms,
   three folds and four seeds, 36 normalized fits passed exact independent
   replay and 36 unchanged fits passed reuse verification. Adding normalized
   FiLM to either signed Δ readout lowered mean Δ PCC and raised mean Δ RMSE
   relative to the corresponding readout alone. Round 4 retained the original
   Gaussian prediction head and carried forward the independent signed Δ
   readout as one of its references.
   [Initial failure record](06.3_round3/log/356-blocked.json) ·
   [Round 3 results](06.3_round3/README.md#formal-results) ·
   [Normalized acceptance](06.3_round3/log/356-formal.json)
4. **06.4 Measured encoder targets and decoder sharing.** The five-arm
   comparison used three folds and three seeds, with 27 new fits and 18
   verified reuses. The arm with an independent signed Δ decoder, `XΔ-XΔΔ`,
   had the highest mean Δ PCC (0.429835). The arm with a shared signed Δ
   decoder, `XΔ-XΔ`, had the lowest mean Δ RMSE (0.240752). For the two
   prespecified paired contrasts, the sample SD exceeded the absolute mean
   difference on all four reported metrics.
   [Round 4 results](06.4_round4/README.md#formal-results) ·
   [Formal acceptance](06.4_round4/log/359-formal.json)

Step 07 froze `XΔ-XΔ` as the CarPS-Transcriptome alpha. It encodes measured
control and measured signed Δ, and uses one signed Δ decoder for both
reconstruction and prediction. The alpha verifier loaded the nine accepted
shared-decoder checkpoints and exactly reproduced the validation objective,
validation predictions, and held-out expression and Δ arrays under the fixed
KPGT and standardized log-dose protocol.
[Architecture decision](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) ·
[Alpha checkpoint and replay evidence](../07_cpst_alpha/README.md#checkpoint-and-replay-evidence)

## Rebuild the comparison figures

The four rounds publish learning curves, held-out scores and differences
from the reference within each fold and seed. The weighted signed Jaccard
scores use accepted held-out predictions and matched controls. The other
scores come from the accepted reports. Use the study guide for candidate
rebuilds and round-specific replay routes, and the current figure receipts
for CSV-only verification.
[Study guide](AGENTS.md#rebuild-the-comparison-figures) ·
[Current figure verification](../README.md#verify-the-current-figures)
