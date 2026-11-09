# Compare Tahoe pretraining and tenv fine-tuning across the Step 03 models

**Status:** proposed (2026-11-09). The shared experiment boundary is defined. The remaining training and evaluation decisions are under design review.

`stage2/07_models_transfer/` compares PRnet, CPA, chemCPA, CMonge and PerturbDiff, the five learned entrants from [Step 03](0042-tenv-model-benchmark-redesign.md). For each model and seed, compare the Tahoe-pretrained predictor with its tenv-fine-tuned descendant. Use the same frozen [Step 05 cold-cm-drug fold 0](0047-perturbdiff-tahoe-transfer-to-tenv.md) for all models at seeds `17`, `42` and `2026`. Tahoe can contain drug blocks held out in tenv under this exposure policy.

Use the Step 03 Morgan representation with a jointly trained, model-specific linear chemical adapter. Use the [Step 06 continuous log dose](0049-perturbdiff-transfer-optimization-comparisons.md) convention in both stages. For positive doses, apply an affine projection to `log10(d / 0.05 µM)`. Use a separate learned control vector at zero dose. The common scale lets the source dose encoder carry over to tenv and evaluates positive doses directly, including 10 µM.

Fine-tuning uses the Step 05 full-parameter update policy and base learning rate. How that policy applies to model-specific optimizer groups and schedules is pending design review.

Deliver checkpoints under `output/checkpoints/` and the figure `output/models_transfer_metrics.png`. Each metric has 12 bars, with two references and a `pretrain` and `finetune` pair for each model. The metric set, checkpoint selection and aggregation are pending design review.
