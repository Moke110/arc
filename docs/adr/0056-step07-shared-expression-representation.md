# Use a shared frozen expression representation in Step 07

**Status:** accepted (2026-11-10).

Compare matrix compression and reconstruction while holding the expression
representation interface fixed. For each profile, add a fixed 256-channel
gene-position lookup to a shared continuous expression MLP with widths
1 → 64 → 256 and GELU. This produces an 18,164 × 256 transcriptome
representation matrix on the frozen CM gene axis. The shared restorer
subtracts the same gene positions and reads one continuous value per gene
with a shared MLP. Calibrate the generator and restorer jointly on individual
expression values from the joint training split, then freeze the same
weights for every entrant. Report the validation error of this round trip
before compression so its contribution can be distinguished from the
end-to-end reconstruction error. A continuous nonlinear map preserves
within-bin information and can provide distinct channel responses.

| Model | Encoder and decoder |
| --- | --- |
| `PCA_1024` | Fit one gene-axis basis from centered training representations pooled across channels. Project every channel onto the same 1,024 components, then apply the inverse projection and restore the center. |
| `MLP_1024` | Use separate encoder and decoder MLPs that map along the token axis. Each shares its weights across all 256 channels. |
| `AAEc_1024` | Use one gene-to-latent and one latent-to-gene cross-attention layer. |
| `AAEcs_1024` | Add one latent self-attention layer after encoder cross-attention and one before decoder cross-attention. |
| `AAEcs_1024+64` | Use the same attention design, with 1,088 queries and joint encoder self-attention. Pass only the first 1,024 reconstruction tokens to the decoder. |

All five decoders receive 1,024 × 256 reconstruction tokens and return an
18,164 × 256 matrix to the restorer. Joint encoder self-attention allows the
extra 64 queries to learn through their effect on reconstruction tokens.
The ablation measures the effect of allowing this interaction. The training-gene
mean remains the expression-space reference and relative-error denominator.

Use the same effective batch, training order and optimizer update cap for
the neural models, with a common budget fixed after resource smoke tests and
before official fitting. Select checkpoints on tenv validation MSE. Keep
the split, expression scale, and source-specific scoring from
[ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md). Commit the
validation seal before test scoring, independently replay the results, and
publish the new active figure as `models_metric.png`. The figure reuses
the previously reported test plates. Historical reports and seals retain
their replay contracts.

This revises the active comparison in
[ADR 0055](0055-step07-matrix-mlp-reconstruction.md) and defines a new
attention ablation building on
[ADR 0054](0054-step07-separate-reconstruction-and-nuisance-tokens.md).
