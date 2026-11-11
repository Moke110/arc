# Publish and replay Step 06

The AP@100 publication and promotion reports are retired. Historical retrieval
figures and score tables remain only in Git. The shared evaluator now reports
only tenv `AUPRC`, `AUPRC_positive` and `AUPRC_negative` from the frozen predictions.

Run under the Step 06 scheduler and repository memory limits:

```bash
python stage2/evaluation/publish.py --step 6
python stage2/evaluation/evaluate.py replay --step 6
```

Numeric tables and the hash manifest live in `output/retrieval/`; PNG/SVG
result figures live directly in `output/`. The publisher independently scores
twice before replacement. Source exposure is a separate, unchanged figure:

```bash
python stage2/06_pdiff_transfer_optimize/execution/script/exposure.py
```

Training and checkpoint provenance remains in `paired_training_provenance.json`.
The source audit, selected native diffusion-loss checkpoints, prediction arrays
and resource receipts remain available. Historical promotion thresholds must
not be reapplied to the new metrics.
