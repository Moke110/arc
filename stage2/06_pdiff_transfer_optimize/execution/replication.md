# Fold 1 paired test replication

The frozen confirmation plan controls whether replication is authorized.
A skipped plan rejects fold 1 fitting artifacts. Existing source preparation,
native fits and independent reloads retain their provenance in
[confirmation.md](confirmation.md); this reporting migration authorizes no
additional training.

Each pair shares the same Tahoe source and mapped bridge. Dose reset changes
only the three native positive-dose rows; learning rate changes only the base
rate. Native loss selection and frozen equal-step EMA remain separate.
Prediction receives only ordered cases, same-plate controls, gene order and
released fingerprints.

The old AP@100 replication scoring, finalization reports and promotion gates
are retired. Historical retrieval scores remain only in Git. Current reporting
uses the shared tenv significance evaluator and its three metrics:

```bash
python stage2/evaluation/publish.py --step 6
python stage2/evaluation/evaluate.py replay --step 6
```

The preserved CPU fixture checks native selection, equal-step EMA, shared Tahoe
bridge, numbered checkpoints and independent target-free prediction replay.
Its recorded run used two threads, 10 GiB memory and 8 GiB swap, took 36.65
seconds and peaked at 4,141,488 KiB RSS.
