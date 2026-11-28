# Compare Tahoe-to-tenv module transfer in CPST Step 11

**Status** Accepted design, 2026-11-28.

Create `11_cpst_transfer_module` to measure how Tahoe pretraining changes CPST
performance on held-out tenv drug parents. Use the selected `XΔ-XΔ` backbone,
objective, ordered Denoise-DEG 5k axis, frozen UniMol and `no-dose` head.

## Training releases

In Step 03, number tenv HVG-DEG 2k, 3k, 4k and 5k as `03.1`–`03.4` and tenv
Denoise-DEG 2k, 3k, 4k and 5k as `03.5`–`03.8`. Preserve the eight accepted
tenv releases and their fold memberships under the new paths. Publish
`03.9_tahoe_denoise-5k` on the same ordered gene axis as tenv Denoise-DEG 5k.
Retire the active `03.2_tahoe_hvg-5k` substep after switching active consumers
and verifying the successor. Keep the sealed historical evidence usable with
its pinned Git revision and materialized inputs.

Build Tahoe directly from the Step 01 integer corpus. Project to the frozen
axis and sum eligible treatment counts within `(plate, cell_line, drug_parent,
dose_um)`. Sum distinct compatible same-plate, same-cell-line DMSO_TF counts
once per control pool. Apply selected-axis CP10K and natural `log1p` separately
to each treatment and its matched control. Compute signed Δ by subtraction.
Keep member identities and exclusions auditable. The plate-level unit preserves
one matched-control interpretation per row and aligns the source and target
training inputs. This replaces the active Tahoe well-level publication boundary.

## Training and comparison

Freeze one Tahoe warm-start split with about 15% of parent–dose combinations
in validation. Assign each combination as a whole across plates and cell lines.
Every validation parent retains another dose in training. Use the Tahoe
validation objective for source checkpoint selection. Tahoe pretraining may
include the parents in any tenv partition. The tenv study uses five new
parent-level cold-drug folds. Keep all doses and plate-level rows of a parent
together. Rotate the five parent groups through validation and test once each,
with the other three groups used for training. Freeze both split ledgers before
the formal fits.

Use seeds 17, 42, 43, 2026 and 2027. For each seed, select one Tahoe checkpoint
and share it across all five tenv folds and all Tahoe-dependent arms.

| Arm | Starting state and fitting data |
| --- | --- |
| `tenv-only` | Random initialization, followed by tenv training. |
| `tahoe-only` | Tahoe checkpoint evaluated directly on tenv. |
| `transfer-x1` | Copy the X1 encoder, its Gaussian mean and log-variance heads, and the X1 decoder. Fine-tune on tenv. |
| `transfer-Δ` | Copy the Δ encoder, its Gaussian mean and log-variance heads, and the shared Δ decoder. Fine-tune on tenv. |
| `transfer-mol` | Copy the learned `Xmol → Zmol` projection, including its BatchNorm state. Fine-tune on tenv. |
| `transfer-pred` | Copy both affine `Zpred` heads, including their weights and biases. Fine-tune on tenv. |
| `transfer-full` | Copy the full model state. Fine-tune on tenv. |

Copy the parameters and buffers of selected modules. Initialize the rest from
the paired `tenv-only` state and match the recipient training random streams
within each fold and seed. Start a fresh Adam optimizer for each tenv fit and
train all modules. Both stages use learning rate `1e-4`, weight decay `1e-5`,
batch size `64`, at most `1000` epochs and early-stopping patience `30`. Select
each fit at its minimum within-domain validation objective.

Score all seven arms on the same held-out tenv conditions with measured
same-plate control. Use full-axis Δ PCC and Δ RMSE as primary metrics and
calibrate Δ MMD on the tenv training partition of each fold. Pair each arm
with `tenv-only` within fold and seed. Report the 25 values per arm, their
equal-weight mean and sample SD, and the observed overlap of held-out tenv
parents with Tahoe training and validation. The matrix comprises five Tahoe
fits, 150 tenv fits and 175 tenv evaluation coordinates. Because CPST shares
the Δ decoder, `transfer-Δ` initializes a path used by both reconstruction and
prediction. The full and single-module arms can capture different interactions.

## Related decisions

- [ADR 0085](0085-stage2-raw-and-training-data.md) for the plate-matched training product.
- [ADR 0099](0099-carps-vae-alpha-and-chemical-representations.md) for the backbone.
- [ADR 0102](0102-stage2-axis-selection-and-comparison-figures.md) for the ordered axis and comparison figures.
- [ADR 0104](0104-carps-vae-no-dose.md) for the selected `no-dose` head.
