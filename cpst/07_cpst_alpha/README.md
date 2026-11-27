# Step 07 CarPS-Transcriptome alpha

Step 07 freezes the `XΔ-XΔ` arm from Step 06.4 as the CarPS-Transcriptome
alpha. The shared signed Δ decoder gave the lowest mean held-out Δ RMSE in
that comparison. The release keeps the ordered tenv Denoise-DEG 5k axis,
measured same-plate DMSO, frozen KPGT features and training-fold standardized
log dose. [Step 06.4 results](../06_tsg_shape/06.4_round4/README.md#formal-results) ·
[ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md)

## Fixed model and public API

Separate 100-dimensional Gaussian encoders model measured control (`X1`) and
measured signed Δ (`X2 − X1`). The prediction head uses the sampled control
latent, a 400-wide projection of the 2,304 KPGT inputs and 32 copies of the
standardized log dose. One signed Δ decoder reconstructs measured Δ and
decodes the predicted latent. Adding predicted Δ to measured `X1` gives
predicted treated expression on the same 5,000 genes. The objective sums
control reconstruction MSE, measured Δ reconstruction MSE, predicted Δ MSE
and 1.75 times predicted Δ L1. It adds 0.1 times the sum of the two KL terms
to the standard normal prior and the KL term from the predicted Gaussian to
the measured Δ Gaussian.

The `cpst_alpha` module exports `AlphaAdapter` and `objective()`. The adapter
takes the pinned TranSiGen source and a fold's dose transform. `build(seed,
device)` initializes the model. `predict(model, inputs, rows)` uses the
matched controls, raw doses and native-width chemical vectors in frozen
condition order to return sampled treated expression in float32. The
default uses the original KPGT initialization. A study can supply
`chemical_width` and optionally a `projection_factory` that returns a CPU
`torch.nn.Linear` of the declared input width and 400 outputs. The adapter
replaces the first chemical layer while preserving the random stream for
latent sampling.
[Public API](output/api/cpst_alpha.py) · [Shared training APIs](../util/README.md#transigen-training-and-evaluation) ·
[TranSiGen source provenance](../../pert_model/TranSiGen/author_log.json)

## Checkpoint and replay evidence

The alpha uses the accepted Step 06.4 `XΔ-XΔ` checkpoints. The release
manifest pins the alpha defaults and runtime sources, the original Step 06.4
code revision and digests, and its formal acceptance receipt.
[Alpha release manifest](output/release.json) ·
[Step 06.4 formal receipt](../06_tsg_shape/06.4_round4/log/359-formal.json)

The alpha verifier loaded the sealed checkpoints for all nine coordinates
across cold-drug folds 0–2 and training seeds 17, 42 and 2026. It reproduced
validation predictions and the validation objective exactly at sampling seed
22617. At sampling seed 100042 it reproduced held-out treated expression and
both predicted and measured signed Δ arrays exactly. The parity receipt
records the registration, checkpoint and prediction hashes for every
coordinate. [Parity receipt](log/362-alpha.json) ·
[Original CUDA resource record](log/362-resources.json)

The verifier uses the materialized Step 06.4 release and its frozen inputs.
The owner guide gives the commands and input requirements for alpha replay.
[Replay verifier](execution/verify_alpha.py) · [Owner guide](AGENTS.md)

## Downstream use

Step 08 used this alpha as the fixed reference for seven chemical
representations, reusing the nine accepted KPGT fits. Its comparison retained
KPGT under the prespecified rule for paired variability. UniMol led the mean
held-out Δ PCC and Δ RMSE, and ADR 0100 adopted it for subsequent studies.
[Step 08 results](../08_cpst_chem-emb/README.md#formal-results) ·
[UniMol adoption](../../docs/adr/0100-carps-vae-unimol.md)

Step 09 selected Denoise-DEG 5k with UniMol and inherited standardized log
dose. Step 10 selected `no-dose` on the cold-dose-drug folds. The selected
CarPS-Transcriptome configuration therefore uses UniMol and `no-dose`,
while the Step 07 release replays the original KPGT and log-dose checkpoints.
[Step 09 results](../09_cpst_axis/README.md#formal-results) ·
[Step 10 results](../10_cpst_dose/README.md#formal-results) ·
[Selected CarPS-Transcriptome configuration](../README.md#selected-carps-transcriptome-configuration)
