# Extend the Step 06.1 encoder-width comparison with sealed prior fits

**Status:** accepted (2026-11-26).

Extend [ADR 0094](0094-transigen-shape-round1.md) with encoder widths 2,000
and 4,000. Both control and treated encoders use the indicated width with
independent parameters. Each added shape uses one 800-unit decoder hidden
layer and expression latent dimension 100. Keep the frozen tenv cold-drug
folds, seed 42, Denoise-DEG 5k, Step 05.2 objective and matched-control
protocol used in the accepted comparison.

| Accepted fit identity | Current coordinate |
| --- | --- |
| `widen_z-100` | `original` |
| `widen_enc` | `widen_enc-6k` |
| `deepen-both` | `deepen_both` |

Fit and independently replay the six new coordinates. Reuse the 24 prior
fits from the accepted Step 06.1 release, including the three `original`
fits that passed exact Step 05.2 parity and serve as the within-fold
reference. The extension accepts sealed prior registrations and replay
evidence from their pinned release after validating compatibility of inputs,
training protocol, and evaluation. Keep their source identities and hashes
in the extension provenance. The canonical names are used for current
registrations and reporting. Historical registrations retain the spellings
that were bound into their original receipts.

Publish the ten-shape report with per-fold scores and learning curves,
equal-fold means and paired differences for all nine metrics, and CSV-only
figure replay. This combines comparable runs while limiting training to the
two newly tested widths.
