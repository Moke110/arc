# Positive-dose initialization screen

The intervention restores the three same-seed native tenv positive-dose
embedding rows after loading the Tahoe bridge, preserving the zero-dose row
and remaining recipient state. Its comparator is the original Transfer
validation EMA selected by native diffusion loss.

The original AP@100 scoring seals, promotion gates and score commands are
retired. Their historical records remain only in Git. Existing fitted models
and sealed predictions are preserved; this metric migration does not rerun
training or reinterpret earlier promotion decisions.

All current tenv retrieval reporting uses the shared three-metric evaluator:

```bash
python stage2/evaluation/publish.py --step 6
python stage2/evaluation/evaluate.py replay --step 6
```

See [publication.md](publication.md) for output placement and replay.
