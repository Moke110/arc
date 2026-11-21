# 0080. Use frozen KPGT for the tenv dose-input comparison

## Context

ADR 0078 established a sealed dual-source, axis-local G_5000 dose-input comparison with ECFP4. Step 03 compared drug representations under fixed FiLM on tenv. KPGT produced a parent-equal OOF Δ PCC of 0.28279 and a paired difference of +0.01370 against ECFP4C in that experiment.

## Decision

Create a tenv-only successor to the ADR 0078 experiment. Reuse Step 03's frozen 2304-dimensional float32 KPGT vectors for all seven TranSiGen dose methods and the drug-only KPGT_MLP baseline. Bind the vector artifact and its provenance in the new input manifest. Match all 459 tenv parents on parent InChIKey and canonical SMILES. Stop before fitting and report any invalid, duplicate, missing or mismatched mappings. Train response models from scratch. Retain Mean_cell and keep both baselines outside method selection.

Preserve the tenv axis-local G_5000 inputs, matched controls, exact three-fold memberships for both scenarios, seed 42, sample-equal evaluation, validation-locked cold-drug selection and fixed cold-drug-dose comparison. The matrix has 42 TranSiGen fits, 6 KPGT_MLP fits and 6 Mean_cell evaluations. Use separate input, run and report identities for this generation. Preserve the sealed ADR 0078 artifacts and replay contracts used by other steps.

## Consequences

The successor measures dose-method performance on tenv with KPGT held fixed. A representation effect under the Step 04 protocol would require a matched ECFP4 comparison. The successor publishes comparison evidence and has no downstream response-weight release.
