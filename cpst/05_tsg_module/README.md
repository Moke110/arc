# Step 05 TranSiGen modules

This study compares five modules on the same tenv Denoise-DEG 5k training
release, measured same-plate DMSO, frozen KPGT and fold-local standardized
log dose. Each module uses the five shared cold-drug folds and training
seeds 17, 42, 43, 2026 and 2027, for 125 new coordinates.
[ADR 0110](../../docs/adr/0110-cpst-transigen-module-rerun.md) ·
[Training and evaluation commands](AGENTS.md)

| Module | Measured encoders | Prediction |
| --- | --- | --- |
| `XX-XX` | Control and treated expression | Treated expression |
| `XX-XXΔ` | Control and treated expression | Independent signed Δ decoder |
| `X-XΔ` | Control only | Signed Δ decoder |
| `XΔ-XΔ` | Control and measured signed Δ | Shared signed Δ decoder |
| `XΔ-XΔΔ` | Control and measured signed Δ | Independent signed Δ decoder |

Training and checkpoint selection minimize the sum of control reconstruction
MSE, signed Δ reconstruction MSE, predicted Δ MSE and 0.1 times the KL sum.
The control-only module uses the control prior KL, control reconstruction and
predicted Δ MSE. The other modules use both prior KL terms and a
predicted-to-measured Gaussian alignment KL. For the two treated-expression
reconstruction modules, subtract the same measured control from reconstructed
and observed treated expression. Their squared reconstruction error therefore
equals the treated-expression squared error algebraically.

The comparison figure subtracts each metric's mean over all valid displayed
`XX-XX` coordinates from every displayed fold-and-seed score. The reference
group has its own centered dots and sample SD. The independent primary
contrasts still subtract scores within the same fold and seed.

## Bounded validation

All five modules passed their two-epoch CUDA fit and exact fresh-process replay
on the complete first fold with seed 17. An independent calculation reproduced
the selected checkpoints' full validation objectives. The bounded report
independently reproduced 35 scores and the paired tables. Two CSV-only renders
reproduced all three figures. The formal register-only run verified 125 new
registered coordinates. [Validation receipt](log/410-gate.json)

## Historical reports

The prior loss studies and the four module rounds remain available at the
revision before the Step 05 cutover. Use their pinned scientific runtime
and saved inputs for historical replays. The existing Step 07–10 model
decisions derive from those historical releases.
[Loss study](https://github.com/Moke110/carps/blob/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst/05_tsg_loss/README.md) ·
[Module rounds and width extension](https://github.com/Moke110/carps/blob/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst/06_tsg_module/README.md) ·
[Local evidence locations](../util/training/verification.md)
