# Compare Tahoe training strategies on a shared tenv cold-drug split

**Status:** proposed (2026-11-08). The comparison scope is set. The split
identity, cross-source drug exclusions and training protocol are open.

`stage2/05_pdiff_transfer/` asks whether Tahoe pretraining improves
PerturbDiff predictions of tenv drug responses. All test responses come from
tenv. Every arm uses the same cold-drug test membership and the comparison
uses six fixed seeds.

The arms are `no_change`, `delta_mean`, `Tahoe` (Tahoe-only training), `CM`
(tenv-only checkpoints from Step 04), `Transfer` (Tahoe pretraining followed
by tenv fine-tuning with separate stage recipes), and `Mix` (pooled Tahoe and
tenv training in a single fit). Deliver validation-selected checkpoints for the
trained models, learning curves, held-out metrics and PR curves, following the
presentation scope of Step 04.

[ADR 0046](0046-perturbdiff-standalone-training.md) supplies only three `CM`
seeds per tenv split. Across the two releases, each tenv test fold also shares
parent-connectivity blocks with Tahoe. These affect how the six paired slots
can be filled and what the cold-drug comparison measures.
